# rez_seasonality_batch.R
# Full-species per-REZ (+ statewide) seasonality run.
#
# Per species: one shared NB GAM (production effort+habitat structure + an
# ordered-factor by-REZ cyclic seasonal smooth); seasonality metrics per
# species x region by posterior simulation (see ebirdabund/R/seasonality.R).
# Abundance (mean_abund, se_abund) and the WLAB crosswalk are joined in the
# MAIN process, so parallel workers never touch terra or the crosswalk.
#
# Robustness: parallel workers (PSOCK) with per-species tryCatch; results
# appended per chunk; resumes on restart (species already in the output CSV are
# skipped); the cluster is rebuilt if a chunk kills its workers.
#
# Usage:  Rscript rez_seasonality_batch.R
#   EBIRD_CONFIG   override config file (default config.yaml)
#   SEASO_NCORES   override worker count (default: config n_cores, capped 5)
#   SEASO_SPECIES  comma-separated species subset (default: full species list)

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf); library(terra); library(mgcv); library(dplyr); library(yaml)
  library(parallel)
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
ALIAS_PATH <- "taxonomy/reidbaker_name_aliases.csv"
OUT_DIR    <- "rez_seasonality"
OUT_CSV    <- file.path(OUT_DIR, "seasonality_all.csv")
RUN_LOG    <- file.path(OUT_DIR, "seasonality_run_log.csv")

MIN_POS   <- 20L
N_SIM     <- 1000L
RATIO_THR <- 2
PROB_THR  <- 0.95

STATEWIDE  <- "statewide"
TARGET_REZ <- c(central_west = "Central-West Orana",
                new_england  = "New England",
                south_west   = "South West")
REZ_LEVELS <- c(STATEWIDE, names(TARGET_REZ))

N_CORES <- as.integer(Sys.getenv("SEASO_NCORES",
                                 min(5L, as.integer(cfg$n_cores))))

dir.create(OUT_DIR, showWarnings = FALSE)

# Final CSV column order (header written once; workers produce the seasonality
# block, main process adds taxonomy + abundance + provenance).
OUT_COLS <- c("common_name", "scientific_name", "wlab_taxon_id",
              "wlab_scientificname", "risk_assessed", "wlab_review",
              "wlab_review_reason", "rez", "n_checklists_rez", "n_positive_rez",
              "sufficient", "seasonality_index", "amplitude_link", "peak_doy",
              "peak_date", "prob_seasonal", "is_seasonal", "window_start_doy",
              "window_end_doy", "window_start_date", "window_end_date",
              "window_wraps", "mean_abund", "se_abund", "range_masked_in_rez",
              "n_train_total", "dev_expl", "run_date")

# ── Study polygon + statewide border ──────────────────────────────────────────
message("Building study polygon...")
aus       <- geodata::gadm(country = cfg$study_polygon$country, level = 1,
                           path = CACHE)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
PP <- cfg$study_polygon$proj_processing
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_sf, PP),
                cfg$study_polygon$buffer_metres), 4326)
parts <- sf::st_cast(sf::st_geometry(region_sf), "POLYGON", warn = FALSE)
xmax  <- vapply(parts, function(g) sf::st_bbox(g)[["xmax"]], numeric(1))
nsw   <- sf::st_transform(sf::st_union(parts[xmax < 155]), 4326)

ebd_files  <- file.path(RAW, cfg$ebd_files)
samp_files <- file.path(RAW, cfg$sampling_files)

# ── REZ polygons + one-time checklist -> REZ membership lookup ─────────────────
message("Deriving REZ membership (once, from cache)...")
rez      <- sf::st_read(REZ_PATH, quiet = TRUE) |> sf::st_transform(4326)
rez_keep <- rez[rez$REZ_Name %in% TARGET_REZ, ]

base_f  <- list.files(ZCACHE, pattern = "zerofilled_.*[.]rds", full.names = TRUE)[1]
base    <- readRDS(base_f)[, c("checklist_id", "latitude", "longitude")]
base    <- base[!is.na(base$latitude) & !is.na(base$longitude), ]
old <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
pts <- sf::st_as_sf(base, coords = c("longitude", "latitude"), crs = 4326)
j   <- sf::st_join(pts, rez_keep["REZ_Name"], left = FALSE)
sf::sf_use_s2(old)
key_of <- setNames(names(TARGET_REZ), unname(TARGET_REZ))
REZ_LUT <- setNames(unname(key_of[as.character(j$REZ_Name)]), j$checklist_id)
REZ_LUT <- REZ_LUT[!duplicated(names(REZ_LUT))]
message(sprintf("  %d checklists mapped to a target REZ.", length(REZ_LUT)))

# ── Taxonomy crosswalk (built once, joined in main process) ────────────────────
tax  <- read.csv(cfg$taxonomy_file, stringsAsFactors = FALSE)
spl  <- read.csv(cfg$species_list_file, stringsAsFactors = FALSE)
reid <- read.csv(REID_PATH, stringsAsFactors = FALSE)
ali  <- read.csv(ALIAS_PATH, stringsAsFactors = FALSE)
CW   <- build_taxonomy_crosswalk(spl$common_name, tax, reid, ali)
rownames(CW) <- CW$common_name

# ── Abundance stacks (read once, main process only) ───────────────────────────
abd_stack <- terra::rast(STACK_PATH)
se_stack  <- terra::rast(SE_PATH)
region_vect <- c(
  setNames(list(terra::vect(nsw)), STATEWIDE),
  setNames(lapply(names(TARGET_REZ),
                  function(k) terra::vect(rez[rez$REZ_Name == TARGET_REZ[[k]], ])),
           names(TARGET_REZ)))

# Region-mean abundance + SE, read off the prediction stack. Returns a length-3
# numeric: mean_abund, se_abund, and range_masked_in_rez (0/1/NA as a flag).
# Three cases:
#   * stem absent from the stack (species not modelled: below the positive-
#     checklist threshold, or excluded) -> NA / NA / NA. No surface exists.
#   * stem present but every cell inside the region is NA (the region falls
#     wholly outside the species' range mask) -> 0 / 0 / 1. Treated as a hard
#     absence so a truly-zero region reads as 0, not NA; the flag records that
#     the zero came from the range mask rather than a modelled ~0.
#   * otherwise -> modelled mean over the region / mean SE / 0. Masked cells
#     within a partially-covered region contribute 0 to the numerator (they are
#     counted in nin), matching the existing whole-region denominator.
region_abd <- function(stem, region) {
  if (!stem %in% names(abd_stack))
    return(c(mean_abund = NA_real_, se_abund = NA_real_, range_masked_in_rez = NA_real_))
  rv <- region_vect[[region]]
  la <- abd_stack[[stem]]; ls <- se_stack[[stem]]
  cr  <- terra::mask(terra::crop(la, rv), rv)
  ins <- terra::mask(terra::crop(terra::setValues(la, 1), rv), rv)
  nin <- as.numeric(terra::global(ins, "sum", na.rm = TRUE)[[1]])
  if (is.na(nin) || nin == 0)              # region doesn't overlap the grid
    return(c(mean_abund = NA_real_, se_abund = NA_real_, range_masked_in_rez = NA_real_))
  n_modelled <- as.numeric(terra::global(!is.na(cr), "sum", na.rm = TRUE)[[1]])
  if (n_modelled == 0)                     # region wholly outside the range mask
    return(c(mean_abund = 0, se_abund = 0, range_masked_in_rez = 1))
  se  <- terra::mask(terra::crop(ls, rv), rv)
  c(mean_abund = as.numeric(terra::global(cr, "sum", na.rm = TRUE)[[1]]) / nin,
    se_abund   = as.numeric(terra::global(se, "mean", na.rm = TRUE)[[1]]),
    range_masked_in_rez = 0)
}

# ── Worker: fit + summarise one species (no terra, no crosswalk) ───────────────
run_species_seasonality <- function(sp) {
  df <- load_ebird(polygon, ebd_files, samp_files, sp, ZCACHE)
  key <- unname(REZ_LUT[df$checklist_id]); key[is.na(key)] <- STATEWIDE
  df$rez <- ordered(key, levels = REZ_LEVELS)
  df <- subsample_hex(df)

  fitobj <- fit_rez_seasonality(df)
  # Deviance explained, with the manual NB fallback used in the production batch
  # (summary(bam)$dev.expl is NaN/NA for some negative-binomial fits).
  dev_expl <- withCallingHandlers(
    tryCatch(summary(fitobj$fit)$dev.expl, error = function(e) NA_real_),
    warning = function(w) invokeRestart("muffleWarning"))
  if (is.null(dev_expl) || is.nan(dev_expl) || is.na(dev_expl)) {
    dev_expl <- tryCatch({
      mod <- fitobj$fit; theta <- exp(mod$family$getTheta())
      y <- mod$y; fv <- mod$fitted.values
      dr <- 2 * (ifelse(y > 0, y * log(y / fv), 0) -
                 (y + theta) * log((y + theta) / (fv + theta)))
      de <- 1 - sum(dr, na.rm = TRUE) / mod$null.deviance
      if (is.nan(de) || is.na(de)) NA_real_ else de
    }, error = function(e) NA_real_)
  }
  n_train  <- nrow(df)

  rows <- do.call(rbind, lapply(REZ_LEVELS, function(region) {
    idx   <- df$rez == region
    summarise_rez_seasonality(
      fitobj, df, region, REZ_LEVELS,
      n_pos = sum(df$observation_count[idx] > 0L, na.rm = TRUE),
      n_chk = sum(idx), n_sim = N_SIM, ratio_thr = RATIO_THR,
      prob_thr = PROB_THR, min_pos = MIN_POS)
  }))
  rows$common_name   <- sp
  rows$dev_expl      <- dev_expl
  rows$n_train_total <- n_train
  rm(fitobj, df); gc(full = TRUE)
  rows
}

# ── Assemble final rows (main process): + taxonomy + abundance + provenance ────
finalise_rows <- function(sp, rows) {
  cw <- CW[sp, ]
  abd <- t(vapply(rows$rez, function(r) region_abd(
    gsub("[^a-z0-9]+", "_", tolower(trimws(sp))), r), numeric(3)))
  out <- data.frame(
    common_name         = sp,
    scientific_name     = cw$scientific_name,
    wlab_taxon_id       = cw$wlab_taxon_id,
    wlab_scientificname = cw$wlab_scientificname,
    risk_assessed       = cw$risk_assessed,
    wlab_review         = cw$wlab_review,
    wlab_review_reason  = cw$wlab_review_reason,
    rows[, setdiff(names(rows), "common_name")],
    mean_abund          = abd[, "mean_abund"],
    se_abund            = abd[, "se_abund"],
    range_masked_in_rez = as.logical(abd[, "range_masked_in_rez"]),
    run_date            = as.character(Sys.Date()),
    stringsAsFactors    = FALSE)
  out[, OUT_COLS]
}

append_rows <- function(df) {
  if (!file.exists(OUT_CSV)) write.csv(df, OUT_CSV, row.names = FALSE)
  else write.table(df, OUT_CSV, sep = ",", col.names = FALSE,
                   row.names = FALSE, append = TRUE, quote = TRUE)
}
log_status <- function(sp, status, msg = NA_character_) {
  row <- data.frame(common_name = sp, status = status,
                    error_message = msg, run_date = as.character(Sys.Date()),
                    stringsAsFactors = FALSE)
  if (!file.exists(RUN_LOG)) write.csv(row, RUN_LOG, row.names = FALSE)
  else write.table(row, RUN_LOG, sep = ",", col.names = FALSE,
                   row.names = FALSE, append = TRUE, quote = TRUE)
}

# ── Species list + resume ──────────────────────────────────────────────────────
species_all <- if (nzchar(Sys.getenv("SEASO_SPECIES")))
  trimws(strsplit(Sys.getenv("SEASO_SPECIES"), ",")[[1]]) else spl$common_name
done <- if (file.exists(OUT_CSV))
  unique(read.csv(OUT_CSV, stringsAsFactors = FALSE)$common_name) else character(0)
species <- setdiff(species_all, done)
message(sprintf("Species: %d total, %d done, %d to run. Workers: %d.",
                length(species_all), length(done), length(species), N_CORES))
if (length(species) == 0L) { message("Nothing to do."); quit(save = "no") }

# ── Parallel cluster ───────────────────────────────────────────────────────────
pkg_src <- normalizePath("ebirdabund", mustWork = FALSE)
export_vars <- c("polygon", "ebd_files", "samp_files", "ZCACHE", "REZ_LUT",
                 "STATEWIDE", "REZ_LEVELS", "N_SIM", "RATIO_THR", "PROB_THR",
                 "MIN_POS", "run_species_seasonality")
setup_workers <- function(cl) {
  parallel::clusterExport(cl, "pkg_src", envir = environment())
  parallel::clusterEvalQ(cl, {
    if (requireNamespace("ebirdabund", quietly = TRUE)) library(ebirdabund)
    else devtools::load_all(pkg_src, quiet = TRUE)
    suppressPackageStartupMessages({ library(sf); library(mgcv); library(dplyr) })
  })
  parallel::clusterExport(cl, export_vars, envir = globalenv())
}

cl <- parallel::makeCluster(N_CORES, timeout = 3600L)
setup_workers(cl)
on.exit(tryCatch(parallel::stopCluster(cl), error = function(e) NULL), add = TRUE)

t0 <- proc.time()[["elapsed"]]
chunks <- split(species, ceiling(seq_along(species) / N_CORES))
n_ok <- 0L; n_err <- 0L

for (ci in seq_along(chunks)) {
  chunk <- chunks[[ci]]
  res <- tryCatch(
    parallel::parLapplyLB(cl, chunk, function(sp)
      tryCatch(run_species_seasonality(sp), error = function(e) e)),
    error = function(e) {
      message(sprintf("  [CLUSTER ERROR] chunk %d (%s); rebuilding.",
                      ci, conditionMessage(e)))
      tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
      cl <<- parallel::makeCluster(N_CORES, timeout = 3600L)
      setup_workers(cl)
      lapply(chunk, function(sp) simpleError(conditionMessage(e)))
    })

  for (i in seq_along(chunk)) {
    sp <- chunk[i]; r <- res[[i]]
    if (inherits(r, "error")) {
      n_err <- n_err + 1L; log_status(sp, "failed", conditionMessage(r))
      message(sprintf("  [FAILED] %s: %s", sp, conditionMessage(r)))
    } else {
      append_rows(finalise_rows(sp, r)); n_ok <- n_ok + 1L
      log_status(sp, "ok")
    }
  }
  el <- proc.time()[["elapsed"]] - t0
  message(sprintf("Chunk %d/%d done | %d ok, %d failed | %.1f min elapsed",
                  ci, length(chunks), n_ok, n_err, el / 60))
}

message(sprintf("\nDone. %d ok, %d failed. Output: %s", n_ok, n_err, OUT_CSV))
