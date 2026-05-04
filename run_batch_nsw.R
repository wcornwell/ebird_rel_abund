suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(data.table)
})

EBD         <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"
SAMP        <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"
CACHE       <- "ebirdabund_cache"
OUTPUT_DIR  <- "species_maps"
LOG_FILE    <- "batch_nsw_log.csv"
BOTW_PATH   <- "botw_species/BOTW_2025.gpkg"
TAXONOMY    <- "nsw_ebird_taxonomy.csv"
GRID_RES_KM <- c(3, 9)

# ── Species list ──────────────────────────────────────────────────────────────
# reporting_rate is already computed against effort-filtered complete checklists
# (see nsw_species_list.R), so this threshold is applied to the correct base.
species_df   <- read.csv("nsw_species_list.csv", stringsAsFactors = FALSE)
species_df   <- species_df[species_df$reporting_rate >= 0.005, ]
species_list <- species_df$common_name
message(sprintf("Species to process: %d (reporting_rate >= 0.5%%)", length(species_list)))

# ── NSW boundary (needed early for bbox in pre-cache block) ───────────────────
message("Getting NSW boundary...")
aus      <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw      <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
nsw_bbox <- as.numeric(sf::st_bbox(sf::st_transform(nsw, 4326)))

# ── Skip already-completed species ────────────────────────────────────────────
safe_name_local <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# ── Migrate log to current schema if it exists ────────────────────────────────
LOG_COLS <- c("common_name", "scientific_name", "status", "exclusion_reason",
              "run_date", "n_checklists", "n_positive", "dev_expl",
              "model_sum", "error_message")
if (file.exists(LOG_FILE)) {
  prior_log <- read.csv(LOG_FILE, stringsAsFactors = FALSE)
  for (col in setdiff(LOG_COLS, names(prior_log))) prior_log[[col]] <- NA
  prior_log <- prior_log[, LOG_COLS]
  write.csv(prior_log, LOG_FILE, row.names = FALSE)
}

done_tif <- sub("[.]tif$", "",
                list.files(file.path(OUTPUT_DIR, "3km"), pattern = "[.]tif$"))
done_log <- character(0)
if (file.exists(LOG_FILE)) {
  prior_log <- read.csv(LOG_FILE, stringsAsFactors = FALSE)
  done_log  <- safe_name_local(
    prior_log$common_name[prior_log$status %in% c("ok", "excluded")]
  )
}
already_done <- union(done_tif, done_log)
species_list <- species_list[!safe_name_local(species_list) %in% already_done]
message(sprintf("Species already completed (skipping): %d",
                length(species_df$common_name) - length(species_list)))

# ── Pre-cache: build missing zerofilled .rds files in one EBD pass ───────────
# Prevents parallel workers competing to fread the large EBD simultaneously.
cached      <- sub("^zerofilled_", "", sub("[.]rds$", "",
               list.files(CACHE, pattern = "zerofilled_.*[.]rds")))
needs_cache <- species_list[!safe_name_local(species_list) %in% cached]

if (length(needs_cache) > 0) {
  message(sprintf(
    "\n── Pre-caching %d species with missing EBD cache (one EBD pass) ──",
    length(needs_cache)
  ))

  # Delegate to package functions so any changes to data-prep logic are picked
  # up automatically. read_sampling() handles column selection, renaming, bbox
  # filtering, and the complete-checklist filter; zero_fill() handles the merge.
  sampling_df <- read_sampling(SAMP, nsw_bbox)

  message("  Reading EBD (full scan for uncached species)...")
  ebd_all <- as.data.frame(fread(
    EBD,
    # Column indices verified against Feb-2026 EBD header:
    #   col 6  = COMMON NAME
    #   col 11 = OBSERVATION COUNT
    #   col 35 = SAMPLING EVENT IDENTIFIER
    select       = c(6L, 11L, 35L),
    sep          = "\t",
    quote        = "",
    showProgress = TRUE,
    na.strings   = ""
  ))
  names(ebd_all) <- c("common_name", "observation_count", "checklist_id")
  ebd_all <- ebd_all[ebd_all[["common_name"]] %in% needs_cache, ]

  message(sprintf("  Building caches for %d species...", length(needs_cache)))
  for (i in seq_along(needs_cache)) {
    sp      <- needs_cache[i]
    cache_f <- file.path(CACHE, sprintf("zerofilled_%s.rds", safe_name_local(sp)))
    if (file.exists(cache_f)) next
    ebd_sp <- ebd_all[ebd_all[["common_name"]] == sp,
                      c("checklist_id", "observation_count"), drop = FALSE]
    saveRDS(zero_fill(sampling_df, ebd_sp), cache_f)
    message(sprintf("    Cached (%d/%d): %s", i, length(needs_cache), sp))
  }
  rm(sampling_df, ebd_all)
  gc()
  message("  Pre-caching complete.\n")
}

# ── Covariates (built once, shared across all workers) ────────────────────────
message("Preparing covariates...")
cov <- prepare_covariates(nsw, cache_dir = CACHE)
message("Covariate layers: ", paste(names(cov), collapse = ", "))

# ── Taxonomy (common name -> scientific name for BOTW range lookup) ───────────
taxonomy <- read.csv(TAXONOMY, stringsAsFactors = FALSE)
message(sprintf("Taxonomy loaded: %d species", nrow(taxonomy)))

# ── Batch run ─────────────────────────────────────────────────────────────────
message(sprintf("\nStarting batch of %d species...", length(species_list)))

if (length(species_list) == 0) {
  message("All species already completed, skipping batch run.")
  log_df  <- if (file.exists(LOG_FILE)) read.csv(LOG_FILE, stringsAsFactors = FALSE) else data.frame()
  t_total <- 0
} else {
  t_start <- proc.time()[["elapsed"]]

  results <- estimate_abundance_batch(
    polygon      = nsw,
    ebird_zip    = EBD,
    sampling_txt = SAMP,
    species_list = species_list,
    taxonomy     = taxonomy,
    cov_stack    = cov,
    cache_dir    = CACHE,
    grid_res_km  = GRID_RES_KM,
    botw_path    = BOTW_PATH,
    output_dir   = OUTPUT_DIR,
    log_file     = LOG_FILE
  )

  t_total <- proc.time()[["elapsed"]] - t_start
}

# ── Build per-resolution raster stacks ───────────────────────────────────────
message("\nBuilding species abundance stacks...")
for (res_km in GRID_RES_KM) {
  res_dir   <- file.path(OUTPUT_DIR, paste0(res_km, "km"))
  tif_files <- sort(list.files(res_dir, pattern = "\\.tif$", full.names = TRUE))

  if (length(tif_files) == 0) {
    message(sprintf("  No .tif files found in %s, skipping stack.", res_dir))
    next
  }

  sp_names   <- sub("\\.tif$", "", basename(tif_files))
  message(sprintf("  %dkm: stacking %d species...", res_km, length(tif_files)))

  # Read each TIF once, then split by layer to avoid opening files twice.
  all_layers <- lapply(tif_files, terra::rast)

  abd_stack <- terra::rast(lapply(all_layers, \(r) r[["abd"]]))
  names(abd_stack) <- sp_names
  abd_path <- file.path(OUTPUT_DIR, sprintf("nsw_abundance_stack_%dkm.tif", res_km))
  terra::writeRaster(abd_stack, abd_path, overwrite = TRUE)
  message(sprintf("  Saved: %s  (%d bands)", abd_path, terra::nlyr(abd_stack)))

  se_stack <- terra::rast(lapply(all_layers, \(r) r[["abd_se"]]))
  names(se_stack) <- sp_names
  se_path <- file.path(OUTPUT_DIR, sprintf("nsw_abundance_se_stack_%dkm.tif", res_km))
  terra::writeRaster(se_stack, se_path, overwrite = TRUE)
  message(sprintf("  Saved: %s  (%d bands)", se_path, terra::nlyr(se_stack)))
}

# ── Summary ───────────────────────────────────────────────────────────────────
log_df   <- read.csv(LOG_FILE, stringsAsFactors = FALSE)
# Summarise only this run's species (prior runs may already be in the log)
run_log  <- log_df[log_df$common_name %in% species_df$common_name[
                     species_df$reporting_rate >= 0.005], ]

hrs      <- floor(t_total / 3600)
mins     <- floor((t_total %% 3600) / 60)
secs     <- round(t_total %% 60)
n_ok     <- sum(run_log$status == "ok")
n_excl   <- sum(run_log$status == "excluded")
n_failed <- sum(run_log$status == "failed")
mean_dev <- mean(run_log$dev_expl[run_log$status == "ok"], na.rm = TRUE)

excl_by_reason <- if (n_excl > 0) {
  excl <- run_log[run_log$status == "excluded", ]
  by_r <- tapply(excl$common_name, excl$exclusion_reason,
                 function(x) paste(x, collapse = ", "))
  paste(names(by_r), by_r, sep = ": ", collapse = "\n  ")
} else "none"

cat(sprintf(
  "\n════════════════════════════════════════\n%d ok  |  %d excluded  |  %d failed  |  %dh %02dm %02ds  |  mean dev.expl = %.3f\nExclusions:\n  %s\nFailed: %s\nLog: %s\n",
  n_ok, n_excl, n_failed, hrs, mins, secs, mean_dev,
  excl_by_reason,
  if (n_failed > 0) paste(run_log$common_name[run_log$status == "failed"], collapse = ", ") else "none",
  LOG_FILE
))
