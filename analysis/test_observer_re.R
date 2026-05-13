# Test: does adding s(observer_id, bs="re") meaningfully change predictions?
#
# Reuses the existing buffered NSW zerofilled cache, but augments one species'
# data with OBSERVER ID read directly from the sampling files (so we don't
# have to rebuild the sampling_master or per-species caches).
#
# Fits the same species twice — without and with the observer random effect —
# at identical subsampling, then predicts on the same grid and compares.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(data.table)
})

# ── Config ────────────────────────────────────────────────────────────────────
SPECIES <- "Superb Fairywren"          # common, well-distributed, in ebirdst
GRID_RES_KM <- 3
RAW_DATA <- "ebirdabund/raw_data"
SAMP <- c(
  file.path(RAW_DATA, "ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-ACT_unv_smp_relMar-2026/ebd_AU-ACT_unv_smp_relMar-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-VIC_unv_smp_relMar-2026/ebd_AU-VIC_unv_smp_relMar-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-QLD_unv_smp_relMar-2026/ebd_AU-QLD_unv_smp_relMar-2026_sampling.txt"),
  file.path(RAW_DATA, "ebd_AU-SA_unv_smp_relMar-2026/ebd_AU-SA_unv_smp_relMar-2026_sampling.txt")
)
CACHE          <- "ebirdabund_cache"
ZEROFILL_CACHE <- "ebirdabund_cache_nsw_buffer"

safe_name_local <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# ── Polygons ──────────────────────────────────────────────────────────────────
message("Loading polygons...")
aus     <- geodata::gadm("AUS", level = 1, path = CACHE)
nsw     <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
# Buffered training polygon (same as run_batch_nsw.R)
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(nsw, 3577), 100000),
  4326
)

# ── Covariate stack (reuses existing tile cache) ──────────────────────────────
message("Preparing covariate stack...")
cov <- prepare_covariates(polygon, cache_dir = CACHE)
cov_cols <- grep(
  "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height)",
  names(cov), value = TRUE
)

# ── Load existing zerofilled and attach observer_id ───────────────────────────
zf_path <- file.path(ZEROFILL_CACHE,
                     sprintf("zerofilled_%s.rds", safe_name_local(SPECIES)))
stopifnot(file.exists(zf_path))
message(sprintf("Loading %s...", zf_path))
zf <- readRDS(zf_path)
message(sprintf("  %d rows (%d detections)", nrow(zf), sum(zf$observation_count != "0")))

message("Reading OBSERVER ID from sampling files...")
obs_lookup <- data.table::rbindlist(lapply(SAMP, function(p) {
  data.table::fread(
    p, sep = "\t", quote = "", showProgress = FALSE, na.strings = "",
    select = c("SAMPLING EVENT IDENTIFIER", "OBSERVER ID")
  )
}))
setnames(obs_lookup, c("checklist_id", "observer_id"))
message(sprintf("  %d sampling events scanned", nrow(obs_lookup)))

zf <- merge(zf, as.data.frame(obs_lookup), by = "checklist_id", all.x = TRUE)
n_na_obs <- sum(is.na(zf$observer_id))
message(sprintf("  %d rows missing observer_id (will be dropped)", n_na_obs))

# ── Apply the same downstream steps as load_ebird's slow path ─────────────────
message("Cleaning, spatial-filtering, extracting covariates, subsampling...")
zf <- ebirdabund:::clean_ebird(zf)
zf_sf <- sf::st_as_sf(zf, coords = c("longitude", "latitude"),
                      crs = 4326, remove = FALSE)
zf <- sf::st_drop_geometry(sf::st_filter(zf_sf, sf::st_transform(polygon, 4326)))
zf <- ebirdabund:::extract_covariates(zf, cov)
zf <- tidyr::drop_na(zf, dplyr::all_of(cov_cols))
ss <- ebirdabund:::subsample_hex(zf, spacing_km = 5)
message(sprintf("Subsampled: %d rows, %d positive",
                nrow(ss), sum(ss$observation_count > 0)))

# ── Fit twice ─────────────────────────────────────────────────────────────────
message("\n==== Fit 1: WITHOUT observer RE ====")
ss_no <- ss[, setdiff(names(ss), "observer_id"), drop = FALSE]
t0 <- Sys.time()
m_no <- ebirdabund:::fit_gam(ss_no)
t_no <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
dev_no <- suppressWarnings(summary(m_no)$dev.expl)
message(sprintf("  fit time: %.1fs, dev expl: %.1f%%", t_no, dev_no * 100))

message("\n==== Fit 2: WITH observer RE ====")
t0 <- Sys.time()
m_with <- ebirdabund:::fit_gam(ss)
t_with <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
dev_with <- suppressWarnings(summary(m_with)$dev.expl)
n_obs <- length(unique(ss$observer_id))
message(sprintf("  fit time: %.1fs, dev expl: %.1f%%, n observers: %d",
                t_with, dev_with * 100, n_obs))

# ── Build prediction surface (use NSW only — that's where we care about REZ) ──
message("\nBuilding prediction surface (NSW)...")
pred_surface <- ebirdabund:::make_prediction_surface(nsw, grid_res_km = GRID_RES_KM)
pred_surface <- ebirdabund:::extract_covariates(pred_surface, cov)
pred_surface <- tidyr::drop_na(pred_surface, dplyr::all_of(cov_cols))

# Use the same peak_doy/peak_time across both models so the comparison is
# purely the effect of the RE, not differences in peak selection.
peak_doy_no   <- NULL  # let predict_abundance pick from m_no
peak_time_no  <- NULL

message("\nPredict: WITHOUT observer RE")
pred_no <- ebirdabund:::predict_abundance(
  m_no, pred_surface, nsw, grid_res_km = GRID_RES_KM
)
message("\nPredict: WITH observer RE (using same peak as above)")
pred_with <- ebirdabund:::predict_abundance(
  m_with, pred_surface, nsw, grid_res_km = GRID_RES_KM,
  peak_doy = pred_no$peak_doy, peak_time = pred_no$peak_time
)

# ── Compare ───────────────────────────────────────────────────────────────────
abd_no   <- terra::values(pred_no$predictions[["abd"]])
abd_with <- terra::values(pred_with$predictions[["abd"]])

stat_line <- function(label, x) {
  q <- stats::quantile(x, c(0.5, 0.9, 0.99), na.rm = TRUE)
  sprintf("  %-25s  mean=%.3f  med=%.3f  q90=%.3f  q99=%.3f  max=%.3f",
          label, mean(x, na.rm = TRUE), q[1], q[2], q[3], max(x, na.rm = TRUE))
}

cat("\n================ SUMMARY ================\n")
cat(sprintf("Species: %s\n", SPECIES))
cat(sprintf("Subsampled n=%d, positive=%d, unique observers=%d\n",
            nrow(ss), sum(ss$observation_count > 0), n_obs))
cat(sprintf("Fit time: no RE=%.1fs, with RE=%.1fs (%.1fx)\n",
            t_no, t_with, t_with / t_no))
cat(sprintf("Dev expl: no RE=%.1f%%, with RE=%.1f%%\n",
            dev_no * 100, dev_with * 100))
cat(sprintf("Peak DOY=%d, peak time=%.2f h (from no-RE model, applied to both)\n",
            pred_no$peak_doy, pred_no$peak_time))
cat("Per-cell abundance over NSW grid:\n")
cat(stat_line("Without observer RE", abd_no), "\n")
cat(stat_line("With observer RE",    abd_with), "\n")
cat(sprintf("Ratio (with / without): mean=%.3f  med=%.3f\n",
            mean(abd_with, na.rm = TRUE) / mean(abd_no, na.rm = TRUE),
            stats::median(abd_with, na.rm = TRUE) / stats::median(abd_no, na.rm = TRUE)))
cat("=========================================\n")

# Save outputs for visual inspection if wanted
saveRDS(list(
  species = SPECIES,
  ss_n = nrow(ss), ss_pos = sum(ss$observation_count > 0),
  n_observers = n_obs,
  fit_time = c(no_re = t_no, with_re = t_with),
  dev_expl = c(no_re = dev_no, with_re = dev_with),
  peak_doy = pred_no$peak_doy, peak_time = pred_no$peak_time,
  abd_summary = list(
    no_re   = summary(abd_no),
    with_re = summary(abd_with)
  )
), "analysis/test_observer_re_result.rds")
message("\nSaved summary to analysis/test_observer_re_result.rds")
