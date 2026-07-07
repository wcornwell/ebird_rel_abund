# rez_seasonality.R  — PROTOTYPE
# Per-REZ (+ statewide) seasonality for a small species set.
#
# For each species it fits ONE negative-binomial GAM that reuses the production
# effort+habitat structure and ADDS an ordered-factor by-REZ cyclic seasonal
# smooth:  s(day_of_year, bs='cc', by = rez).  The reference level ("statewide"
# = rest of state) is carried by the existing global s(day_of_year); each target
# REZ gets a difference smooth, so its own seasonal curve = global + difference.
# Sparse REZs shrink toward the statewide curve (select=TRUE double penalty).
#
# Seasonality metrics per species x region are read off the fitted curve by
# posterior simulation from vcov(fit) (same lower-CI logic as the production
# peak-DOY step, extended to an interval).
#
# Abundance columns (mean_abund, se_abund) are taken from the production
# prediction stacks, exactly as rez_abundance.R computes them.
#
# NOTE: once the columns/approach are locked, the fit + summarise helpers below
# move into ebirdabund/R/seasonality.R and this becomes a thin driver.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf); library(terra); library(mgcv); library(dplyr); library(yaml)
})
set.seed(1)

# ── Config ────────────────────────────────────────────────────────────────────
CONFIG_FILE <- Sys.getenv("EBIRD_CONFIG", "config.yaml")
cfg    <- yaml::read_yaml(CONFIG_FILE)
CACHE  <- cfg$covariate_cache
ZCACHE <- cfg$zerofill_cache
RAW    <- cfg$raw_data

REZ_PATH   <- "RenewableEnergyZones_Spatial/RenewableEnergyZones.shp"
STACK_PATH <- "species_maps/nsw_abundance_stack_3km.tif"
SE_PATH    <- "species_maps/nsw_abundance_se_stack_3km.tif"
REID_PATH  <- "taxonomy/reidbaker_bird_risk_NSW_species.csv"
OUT_DIR    <- "rez_seasonality"
OUT_CSV    <- file.path(OUT_DIR, "seasonality_prototype.csv")

MIN_POS   <- 20L     # per-region positive-detection gate (subsampled data)
N_SIM     <- 1000L   # posterior draws
RATIO_THR <- 2       # "seasonal" if peak >= RATIO_THR x trough on response scale
PROB_THR  <- 0.95    # posterior prob threshold for is_seasonal

# rez factor levels: reference first ("statewide" = rest of state)
STATEWIDE <- "statewide"
TARGET_REZ <- c(central_west = "Central-West Orana",   # "Central West"
                new_england  = "New England",          # "Armidale"
                south_west   = "South West")            # "Hay"
REZ_LEVELS <- c(STATEWIDE, names(TARGET_REZ))

TEST_SPECIES <- c("Pallid Cuckoo", "White-throated Needletail", "Swift Parrot",
                  "Crested Pigeon")   # Crested Pigeon = aseasonal resident control

dir.create(OUT_DIR, showWarnings = FALSE)

# ── Study polygon (for load_ebird bbox + statewide abundance mask) ─────────────
message("Building study polygon...")
aus       <- geodata::gadm(country = cfg$study_polygon$country, level = 1,
                           path = CACHE)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
PP <- cfg$study_polygon$proj_processing
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_sf, PP),
                cfg$study_polygon$buffer_metres), 4326)
# statewide (unbuffered, offshore Lord Howe removed) — for the abundance mask
parts    <- sf::st_cast(sf::st_geometry(region_sf), "POLYGON", warn = FALSE)
xmax     <- vapply(parts, function(g) sf::st_bbox(g)[["xmax"]], numeric(1))
nsw      <- sf::st_transform(sf::st_union(parts[xmax < 155]), 4326)

ebd_files  <- file.path(RAW, cfg$ebd_files)
samp_files <- file.path(RAW, cfg$sampling_files)

# ── REZ polygons + per-checklist REZ membership ───────────────────────────────
rez      <- sf::st_read(REZ_PATH, quiet = TRUE) |> sf::st_transform(4326)
rez_keep <- rez[rez$REZ_Name %in% TARGET_REZ, ]

assign_rez <- function(df) {
  old <- sf::sf_use_s2(); sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(old))
  keep <- !is.na(df$longitude) & !is.na(df$latitude)
  pts <- sf::st_as_sf(df[keep, ], coords = c("longitude", "latitude"),
                      crs = 4326, remove = FALSE)
  j   <- sf::st_join(pts, rez_keep["REZ_Name"], left = TRUE)
  j   <- j[!duplicated(j$checklist_id), ]          # drop overlapping-poly dups
  # map each REZ_Name to its short key; everything else -> statewide
  key_of <- setNames(names(TARGET_REZ), unname(TARGET_REZ))
  lbl    <- key_of[as.character(j$REZ_Name)]
  lbl[is.na(lbl)] <- STATEWIDE
  # join back to df by checklist_id (index-safe)
  lut <- setNames(lbl, j$checklist_id)
  raw <- unname(lut[df$checklist_id]); raw[is.na(raw)] <- STATEWIDE
  df$rez <- ordered(raw, levels = REZ_LEVELS)
  df
}

# ── Taxonomy crosswalk (inline; validates Feature-1 binomial join) ────────────
tax  <- read.csv(cfg$taxonomy_file, stringsAsFactors = FALSE)
reid <- read.csv(REID_PATH, stringsAsFactors = FALSE)
reid$binom <- vapply(strsplit(reid$WLAB_ScientificName, " "),
                     function(x) paste(x[1:2], collapse = " "), character(1))
crosswalk <- function(common) {
  sci <- tax$scientific_name[match(common, tax$common_name)]
  bn  <- paste(strsplit(sci, " ")[[1]][1:2], collapse = " ")
  hit <- reid[reid$binom == bn, ]
  list(
    scientific_name     = sci,
    wlab_taxon_id       = if (nrow(hit)) paste(hit$WLAB_Taxon_ID, collapse = ";") else NA_character_,
    wlab_scientificname = if (nrow(hit)) paste(hit$WLAB_ScientificName, collapse = ";") else NA_character_,
    risk_assessed       = nrow(hit) > 0
  )
}

# ── Model fit: production formula + ordered-factor by-REZ seasonal smooth ──────
fit_seasonality <- function(df) {
  hab_cols <- grep(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height|nightlights|palsar_hv)",
    names(df), value = TRUE)
  hab_cols <- setdiff(hab_cols, "lc_shrubs")
  hab_cols <- hab_cols[vapply(hab_cols, function(c)
    length(unique(stats::na.omit(df[[c]]))) >= 4L, logical(1))]

  base_f <- build_gam_formula(df, hab_cols)            # incl. global s(day_of_year)
  rhs    <- paste(deparse(base_f[[3]]), collapse = " ")
  kdoy   <- safe_k(df$day_of_year, 10L)
  rhs    <- paste0(rhs, " + rez + s(day_of_year, bs = 'cc', k = ", kdoy,
                   ", by = rez)")
  form   <- stats::as.formula(paste("observation_count ~", rhs))

  protos <- levels(df$protocol_type)
  ref    <- if ("Traveling Count" %in% protos) "Traveling Count" else protos[1]
  df$protocol_type <- stats::relevel(df$protocol_type, ref = ref)
  knots  <- list(time_observations_started = c(0, 24), day_of_year = c(0, 365))
  gbic   <- log(nrow(df)) / 2

  try_fit <- function(sel, disc, g) tryCatch(
    mgcv::bam(form, data = df, family = mgcv::nb(), method = "fREML",
              discrete = disc, select = sel, gamma = g, knots = knots),
    error = function(e) NULL)

  fit <- try_fit(TRUE, TRUE, gbic)
  if (is.null(fit)) fit <- try_fit(FALSE, TRUE, 1.4)
  if (is.null(fit)) fit <- try_fit(FALSE, FALSE, 1.4)
  if (is.null(fit)) stop("GAM fitting failed for all attempts.")
  list(fit = fit, hab_cols = hab_cols, ref_proto = ref)
}

# Cyclic contiguous window around the peak where curve >= threshold.
window_doy <- function(above, peak) {
  n <- length(above)
  if (!above[peak]) return(c(NA, NA, FALSE))
  lo <- peak; while (above[(lo - 2) %% n + 1]) { lo <- (lo - 2) %% n + 1; if (lo == peak) break }
  hi <- peak; while (above[hi %% n + 1])       { hi <- hi %% n + 1;       if (hi == peak) break }
  c(lo, hi, lo > hi)   # wraps if start doy > end doy
}

doy_to_date <- function(d) if (is.na(d)) NA_character_ else
  format(as.Date(d - 1, origin = "2023-01-01"), "%d %b")

# Posterior-simulation summary of one region's seasonal curve.
summarise_region <- function(fitobj, df, region, n_pos, n_chk) {
  fit <- fitobj$fit
  hab <- fitobj$hab_cols
  ref <- data.frame(lapply(df[, hab, drop = FALSE],
                           function(x) mean(x, na.rm = TRUE)))
  ref$duration_minutes   <- 60
  ref$effort_distance_km <- 1
  ref$number_observers   <- 1L
  ref$protocol_type      <- factor(fitobj$ref_proto,
                                   levels = levels(df$protocol_type))
  ref$time_observations_started <- 12
  if ("observer_expertise" %in% names(df))
    ref$observer_expertise <- stats::median(df$observer_expertise, na.rm = TRUE)

  sweep <- ref[rep(1L, 365L), , drop = FALSE]
  sweep$day_of_year <- 1:365
  sweep$rez <- ordered(region, levels = REZ_LEVELS)

  Xp    <- predict(fit, newdata = sweep, type = "lpmatrix")
  betas <- mgcv::rmvn(N_SIM, coef(fit), vcov(fit))     # N_SIM x p
  eta   <- Xp %*% t(betas)                             # 365 x N_SIM (link)

  amp_link      <- apply(eta, 2, function(c) max(c) - min(c))
  prob_seasonal <- mean(exp(amp_link) > RATIO_THR)
  eta_med       <- apply(eta, 1, stats::median)
  resp          <- exp(eta_med)
  peak          <- which.max(eta_med)
  suff          <- n_pos >= MIN_POS
  seasonal      <- suff && prob_seasonal >= PROB_THR

  # Bounded seasonality index in [0,1]: 1 - trough/peak on the response scale.
  # ~1 for strict migrants (winter trough ~0), ~0 for a flat resident. Derived
  # from amplitude_link (= log(peak/trough)) so it can't blow up like a raw ratio.
  amp_med <- stats::median(amp_link)
  seas_ix <- 1 - exp(-amp_med)

  win <- if (seasonal) window_doy(resp >= mean(resp), peak) else c(NA, NA, FALSE)

  data.frame(
    rez               = region,
    n_checklists_rez  = n_chk,
    n_positive_rez    = n_pos,
    sufficient        = suff,
    seasonality_index = if (suff) seas_ix else NA_real_,
    amplitude_link    = if (suff) amp_med else NA_real_,
    peak_doy         = if (suff) peak else NA_integer_,
    peak_date        = if (suff) doy_to_date(peak) else NA_character_,
    prob_seasonal    = if (suff) prob_seasonal else NA_real_,
    is_seasonal      = if (suff) seasonal else NA,
    window_start_doy = win[1], window_end_doy = win[2],
    window_start_date = doy_to_date(win[1]), window_end_date = doy_to_date(win[2]),
    window_wraps     = as.logical(win[3]),
    stringsAsFactors = FALSE
  )
}

# ── Abundance from production stacks (rez_abundance.R semantics) ───────────────
abd_stack <- rast(STACK_PATH); se_stack <- rast(SE_PATH)
region_abd <- function(stem, region_vect) {
  la <- abd_stack[[stem]]; ls <- se_stack[[stem]]
  cr  <- mask(crop(la, region_vect), region_vect)
  ins <- mask(crop(setValues(la, 1), region_vect), region_vect)
  nin <- as.numeric(global(ins, "sum", na.rm = TRUE)[[1]])
  se  <- mask(crop(ls, region_vect), region_vect)
  c(mean_abund = as.numeric(global(cr, "sum", na.rm = TRUE)[[1]]) / nin,
    se_abund   = as.numeric(global(se, "mean", na.rm = TRUE)[[1]]))
}
region_vect <- c(
  setNames(list(vect(nsw)), STATEWIDE),
  lapply(names(TARGET_REZ),
         function(k) vect(rez[rez$REZ_Name == TARGET_REZ[[k]], ]))
)
names(region_vect) <- REZ_LEVELS

# ── Main loop ─────────────────────────────────────────────────────────────────
all_rows <- list()
for (sp in TEST_SPECIES) {
  message("\n==== ", sp, " ====")
  stem <- gsub("[^a-z0-9]+", "_", tolower(trimws(sp)))
  cw   <- crosswalk(sp)

  df <- load_ebird(polygon, ebd_files, samp_files, sp, ZCACHE)
  df <- assign_rez(df)
  df <- subsample_hex(df)
  message(sprintf("  %d checklists after subsample; rez counts: %s",
                  nrow(df),
                  paste(names(table(df$rez)), as.integer(table(df$rez)),
                        sep = "=", collapse = ", ")))

  fitobj <- fit_seasonality(df)

  for (region in REZ_LEVELS) {
    idx   <- df$rez == region
    n_chk <- sum(idx)
    n_pos <- sum(df$observation_count[idx] > 0L, na.rm = TRUE)
    s   <- summarise_region(fitobj, df, region, n_pos, n_chk)
    abd <- region_abd(stem, region_vect[[region]])

    all_rows[[length(all_rows) + 1]] <- cbind(
      data.frame(common_name = sp,
                 scientific_name = cw$scientific_name,
                 wlab_taxon_id = cw$wlab_taxon_id,
                 wlab_scientificname = cw$wlab_scientificname,
                 risk_assessed = cw$risk_assessed,
                 stringsAsFactors = FALSE),
      s,
      data.frame(mean_abund = abd[["mean_abund"]],
                 se_abund   = abd[["se_abund"]]),
      data.frame(n_train_total = nrow(df),
                 dev_expl = suppressWarnings(summary(fitobj$fit)$dev.expl),
                 run_date = as.character(Sys.Date()))
    )
  }
}

out <- do.call(rbind, all_rows)
# BH-FDR across all species x REZ rows that were testable
out$q_seasonal <- NA_real_
tst <- !is.na(out$prob_seasonal)
out$q_seasonal[tst] <- p.adjust(1 - out$prob_seasonal[tst], method = "BH")

write.csv(out, OUT_CSV, row.names = FALSE)
message("\nWrote ", OUT_CSV, " (", nrow(out), " rows)")
print(out[, c("common_name", "rez", "n_positive_rez", "sufficient",
              "seasonality_index", "peak_date", "is_seasonal",
              "window_start_date", "window_end_date", "mean_abund")])
