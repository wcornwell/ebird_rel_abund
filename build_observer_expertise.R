# Stage 1 of observer-expertise scoring (Kelling et al. 2015; Johnston et al. 2021).
#
# Fits a single global negative-binomial GAM regressing species count per
# complete checklist on effort covariates with a per-observer random
# intercept. The fitted observer random-effect BLUPs are the expertise score
# and are written to ZEROFILL_CACHE/observer_expertise.rds.
#
# Output is then attached to each species' training data by
# attach_observer_expertise() in load_ebird.R, and the per-species GAM gains
# an s(observer_expertise) smooth. Per-species predictions are made at the
# median expertise across training data.
#
# Run once after the per-species zerofill cache has been built. After this
# completes, DELETE ZEROFILL_CACHE/sampling_master.rds so it rebuilds with the
# new observer_expertise column attached.
#
# Citation:
#   Kelling, S. et al. 2015. Can observation skills of citizen scientists be
#     estimated using species accumulation curves? PLoS ONE 10: e0139600.
#   Johnston, A. et al. 2021. Analytical guidelines to increase the value of
#     community science data: An example using eBird data to estimate species
#     distributions. Diversity and Distributions 27: 1265-1277.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(geodata)
  library(data.table)
  library(mgcv)
})

# ── Config (mirror run_batch_nsw.R) ───────────────────────────────────────────
RAW_DATA <- "ebirdabund/raw_data"
EBD <- c(
  file.path(RAW_DATA, "ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"),
  file.path(RAW_DATA, "ebd_AU-ACT_unv_smp_relMar-2026/ebd_AU-ACT_unv_smp_relMar-2026.txt"),
  file.path(RAW_DATA, "ebd_AU-VIC_unv_smp_relMar-2026/ebd_AU-VIC_unv_smp_relMar-2026.txt"),
  file.path(RAW_DATA, "ebd_AU-QLD_unv_smp_relMar-2026/ebd_AU-QLD_unv_smp_relMar-2026.txt"),
  file.path(RAW_DATA, "ebd_AU-SA_unv_smp_relMar-2026/ebd_AU-SA_unv_smp_relMar-2026.txt")
)
SAMP <- c(
  file.path(RAW_DATA, "ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-ACT_unv_smp_relMar-2026/ebd_AU-ACT_unv_smp_relMar-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-VIC_unv_smp_relMar-2026/ebd_AU-VIC_unv_smp_relMar-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-QLD_unv_smp_relMar-2026/ebd_AU-QLD_unv_smp_relMar-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-SA_unv_smp_relMar-2026/ebd_AU-SA_unv_smp_relMar-2026_sampling.txt")
)
CACHE          <- "ebirdabund_cache"
ZEROFILL_CACHE <- "ebirdabund_cache_nsw_buffer"

# Calibration controls — set conservatively so the fit completes overnight.
MIN_CHECKLISTS_PER_OBSERVER <- 20L     # drop rare observers (BLUPs would be heavily shrunk anyway)
HEX_SPACING_KM              <- 5       # spatiotemporal subsampling (same as per-species pipeline)

OUTPUT_FILE <- file.path(ZEROFILL_CACHE, "observer_expertise.rds")
LOG_FILE    <- "observer_expertise_log.txt"
sink(LOG_FILE, split = TRUE)
on.exit(sink(NULL), add = TRUE)

message(sprintf("Stage 1 calibration started: %s", Sys.time()))

# ── Study polygon (NSW + 100 km buffer) ───────────────────────────────────────
message("Loading polygon (NSW + 100 km buffer)...")
aus     <- geodata::gadm("AUS", level = 1, path = CACHE)
nsw     <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(nsw, 3577), 100000),
  4326
)
bbox <- as.numeric(sf::st_bbox(polygon))

# ── Read sampling + spatial clip ──────────────────────────────────────────────
message("Reading sampling files...")
samp <- ebirdabund:::read_sampling(SAMP, bbox)
message(sprintf("  %d complete checklists in bbox.", nrow(samp)))

samp_sf <- sf::st_as_sf(samp, coords = c("longitude", "latitude"),
                        crs = 4326, remove = FALSE)
samp <- sf::st_drop_geometry(sf::st_filter(samp_sf, polygon))
rm(samp_sf); gc(verbose = FALSE)
message(sprintf("  %d after spatial clip to polygon.", nrow(samp)))

# ── Compute n_species per checklist by scanning all EBD files in one pass ─────
message("Scanning EBD files for per-checklist species counts...")
ebd_paths <- vapply(EBD, ebirdabund:::resolve_ebird_path, character(1))
invisible(lapply(ebd_paths, ebirdabund:::validate_ebd_header))
awk_cmd <- paste("awk 'FNR>1 || NR==1'",
                 paste(shQuote(ebd_paths), collapse = " "))
ebd_dt <- data.table::fread(
  cmd = awk_cmd, sep = "\t", quote = "", showProgress = TRUE, na.strings = "",
  select = c(6L, 35L)   # COMMON NAME, SAMPLING EVENT IDENTIFIER
)
data.table::setnames(ebd_dt, c("common_name", "checklist_id"))
message(sprintf("  Read %d EBD rows.", nrow(ebd_dt)))

ebd_dt <- ebd_dt[checklist_id %in% samp$checklist_id]
n_sp_dt <- ebd_dt[, .(n_species = data.table::uniqueN(common_name)),
                  by = checklist_id]
rm(ebd_dt); gc(verbose = FALSE)

samp <- merge(samp, as.data.frame(n_sp_dt), by = "checklist_id", all.x = TRUE)
samp$n_species[is.na(samp$n_species)] <- 0L
message(sprintf("  n_species per checklist: mean=%.2f, median=%d, max=%d, %%zero=%.1f%%",
                mean(samp$n_species), stats::median(samp$n_species),
                max(samp$n_species), 100 * mean(samp$n_species == 0)))

# ── Effort/time filters (mirror clean_ebird) ──────────────────────────────────
n_before <- nrow(samp)
samp <- samp |>
  dplyr::filter(
    duration_minutes >= 10,
    duration_minutes <= 300,
    is.na(effort_distance_km) | effort_distance_km <= 10,
    number_observers <= 10,
    !is.na(ebirdabund:::time_to_decimal(time_observations_started))
  ) |>
  dplyr::mutate(
    observation_date          = as.Date(observation_date),
    day_of_year               = lubridate::yday(observation_date),
    year                      = lubridate::year(observation_date),
    week                      = lubridate::week(observation_date),
    time_observations_started = ebirdabund:::time_to_decimal(time_observations_started),
    effort_distance_km        = dplyr::if_else(
      is.na(effort_distance_km), 0, effort_distance_km),
    protocol_type             = factor(protocol_type)
  )
message(sprintf("Effort/time filters: %d → %d rows", n_before, nrow(samp)))

# ── Observer activity filter ──────────────────────────────────────────────────
obs_counts <- samp |> dplyr::count(observer_id, name = "n_lists")
keep <- obs_counts$observer_id[obs_counts$n_lists >= MIN_CHECKLISTS_PER_OBSERVER]
n_before <- nrow(samp)
samp <- samp[samp$observer_id %in% keep, ]
message(sprintf("Observer filter (>= %d checklists/observer): %d observers, %d → %d rows",
                MIN_CHECKLISTS_PER_OBSERVER, length(keep), n_before, nrow(samp)))

# ── Spatiotemporal subsampling ────────────────────────────────────────────────
samp <- ebirdabund:::subsample_hex(samp, spacing_km = HEX_SPACING_KM)

# ── Factor setup ──────────────────────────────────────────────────────────────
samp$protocol_type <- stats::relevel(samp$protocol_type,
                                     ref = modal_protocol(samp$protocol_type))
samp$observer_id <- factor(samp$observer_id)

n_obs <- nlevels(samp$observer_id)
message(sprintf("\nFinal calibration set: %d rows, %d observers", nrow(samp), n_obs))

# ── Fit Stage 1 calibration GAM ───────────────────────────────────────────────
message(sprintf("\nFitting calibration GAM (this is the expensive step)..."))
t0 <- Sys.time()
fit <- mgcv::bam(
  n_species ~ s(log(duration_minutes), k = 5) +
              s(log1p(effort_distance_km), k = 5) +
              s(time_observations_started, bs = "cc", k = 7) +
              s(day_of_year, bs = "cc", k = 7) +
              protocol_type +
              s(observer_id, bs = "re"),
  data     = samp,
  family   = mgcv::nb(),
  method   = "fREML",
  discrete = TRUE,
  knots    = list(time_observations_started = c(0, 24),
                  day_of_year               = c(0, 365)),
  select   = FALSE
)
t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
dev_expl <- suppressWarnings(summary(fit)$dev.expl)
message(sprintf("Fit complete in %.0f s (%.2f h); dev expl: %.1f%%",
                t_fit, t_fit / 3600, 100 * dev_expl))

# ── Extract observer BLUPs ────────────────────────────────────────────────────
all_coefs <- stats::coef(fit)
re_idx    <- grep("^s\\(observer_id\\)\\.", names(all_coefs))
stopifnot(length(re_idx) == n_obs)

expertise <- data.frame(
  observer_id = levels(samp$observer_id),
  expertise   = unname(all_coefs[re_idx]),
  stringsAsFactors = FALSE
)
message(sprintf("Extracted %d observer BLUPs.", nrow(expertise)))
message(sprintf("  Score: min=%.3f, median=%.3f, max=%.3f, sd=%.3f",
                min(expertise$expertise), stats::median(expertise$expertise),
                max(expertise$expertise), stats::sd(expertise$expertise)))

# ── Save ──────────────────────────────────────────────────────────────────────
attr(expertise, "n_checklists")              <- nrow(samp)
attr(expertise, "n_observers")               <- n_obs
attr(expertise, "min_checklists_per_observer") <- MIN_CHECKLISTS_PER_OBSERVER
attr(expertise, "hex_spacing_km")            <- HEX_SPACING_KM
attr(expertise, "fit_time_seconds")          <- t_fit
attr(expertise, "dev_expl")                  <- dev_expl
attr(expertise, "built_on")                  <- Sys.time()
attr(expertise, "method") <- paste(
  "Per-observer BLUPs from a global negative-binomial GAM regressing species",
  "count per complete checklist on effort covariates with a random observer",
  "intercept. Following Kelling et al. (2015) and Johnston et al. (2021)."
)

saveRDS(expertise, OUTPUT_FILE)
message(sprintf("\nSaved expertise scores: %s", OUTPUT_FILE))
message("\nNEXT STEP: delete sampling_master.rds so the per-species pipeline")
message("rebuilds the master with observer_expertise attached:")
message(sprintf("  rm %s/sampling_master.rds", ZEROFILL_CACHE))
message(sprintf("Stage 1 calibration finished: %s", Sys.time()))
