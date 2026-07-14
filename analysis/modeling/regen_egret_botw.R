# regen_egret_botw.R
#
# Regenerate the Eastern Cattle-Egret abundance surface with BOTW forced
# (Ardea coromanda -> Bubulcus ibis via botw_name_aliases.csv), because its
# ebirdst range clips inland NSW — most visibly the New England REZ, where it
# has 66 positive checklists but was masked to NA. Overwrites the species TIFs
# at every resolution, rebuilds the abundance + SE stacks, and updates the
# range_source in batch_nsw_log.csv. Re-run rejoin_seasonality_abund.R afterward
# to refresh seasonality_all.csv from the rebuilt stack.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf); library(terra); library(geodata); library(yaml)
})

SP     <- "Eastern Cattle-Egret"
cfg    <- yaml::read_yaml("config.yaml")
CACHE  <- cfg$covariate_cache
ZCACHE <- cfg$zerofill_cache
RAW    <- cfg$raw_data
ebd_files  <- file.path(RAW, cfg$ebd_files)
samp_files <- file.path(RAW, cfg$sampling_files)
OUTPUT_DIR <- "species_maps"
GRID_RES   <- sort(unique(as.numeric(unlist(cfg$grid_resolutions_km))))

botw_aliases <- read.csv("botw_name_aliases.csv", stringsAsFactors = FALSE)
tax <- read.csv(cfg$taxonomy_file, stringsAsFactors = FALSE)
sci <- setNames(tax$scientific_name, tax$common_name)[[SP]]
message(sprintf("Species: %s (%s)  -> forcing BOTW", SP, sci))

# ── Study + prediction polygons (mirrors run_batch_nsw.R) ────────────────────
aus       <- geodata::gadm(country = cfg$study_polygon$country, level = 1, path = CACHE)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
PP <- cfg$study_polygon$proj_processing
# Full buffered polygon (incl. Lord Howe) — drives covariate extraction so the
# bbox-keyed cov_stack_v7 / CHMv2 / PALSAR caches HIT. Must match production
# (run_batch_nsw.R), else prepare_covariates re-streams CHMv2 (~35 min).
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_sf, PP), cfg$study_polygon$buffer_metres), 4326)
# Prediction polygon: mainland only (Lord Howe dropped) — spurious-offshore fix.
region_parts <- sf::st_cast(sf::st_geometry(region_sf), "POLYGON", warn = FALSE)
part_xmax    <- vapply(region_parts, function(g) sf::st_bbox(g)[["xmax"]], numeric(1))
region_main  <- sf::st_union(region_parts[part_xmax < 155])
polygon_pred <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_main, PP), cfg$study_polygon$buffer_metres), 4326)
nsw <- sf::st_transform(region_main, 4326)

message("Loading cached covariate stack (full bbox -> cache hit)...")
cov <- prepare_covariates(polygon, cache_dir = CACHE)

# ── Fit once, predict at each resolution with BOTW forced ────────────────────
model_fit <- fit_species_model(
  polygon = polygon_pred, ebird_zip = ebd_files, sampling_txt = samp_files,
  species = SP, cov_stack = cov, cache_dir = ZCACHE, hex_spacing_km = 5)

range_source <- NA_character_
cached_rv <- NULL; cached_doy <- NULL; cached_time <- NULL
for (res_km in GRID_RES) {
  message(sprintf("\nPredicting at %d km...", res_km))
  pred <- predict_species_map(
    model_fit = model_fit, polygon = polygon_pred, species = SP,
    sci_name = sci, grid_res_km = res_km,
    peak_doy = cached_doy, peak_time = cached_time,
    use_range = TRUE, botw_path = cfg$botw_path, botw_aliases = botw_aliases,
    force_botw = TRUE, range_resolution = "27km", border = nsw,
    range_vect = cached_rv)
  if (is.null(cached_rv)   && !is.null(pred$range_vect)) cached_rv   <- pred$range_vect
  if (is.null(cached_doy)) { cached_doy <- pred$peak_doy; cached_time <- pred$peak_time }
  range_source <- pred$range_source

  stem <- file.path(OUTPUT_DIR, paste0(res_km, "km"), safe_name(SP))
  ggplot2::ggsave(paste0(stem, ".png"), pred$plot, width = 10, height = 9, dpi = 150)
  terra::writeRaster(pred$predictions, paste0(stem, ".tif"), overwrite = TRUE)
  message(sprintf("  range_source = %s; wrote %s.tif", range_source, stem))
}
rm(model_fit); gc(full = TRUE)

# ── Rebuild stacks (mirrors run_batch_nsw.R) ─────────────────────────────────
for (res_km in GRID_RES) {
  res_dir   <- file.path(OUTPUT_DIR, paste0(res_km, "km"))
  tif_files <- sort(list.files(res_dir, pattern = "\\.tif$", full.names = TRUE))
  sp_names  <- sub("\\.tif$", "", basename(tif_files))
  message(sprintf("Rebuilding %dkm stack from %d TIFs...", res_km, length(tif_files)))
  all_layers <- lapply(tif_files, terra::rast)
  abd <- terra::rast(lapply(all_layers, \(r) r[["abd"]])); names(abd) <- sp_names
  terra::writeRaster(abd, file.path(OUTPUT_DIR, sprintf("nsw_abundance_stack_%dkm.tif", res_km)),
                     overwrite = TRUE)
  se  <- terra::rast(lapply(all_layers, \(r) r[["abd_se"]])); names(se) <- sp_names
  terra::writeRaster(se, file.path(OUTPUT_DIR, sprintf("nsw_abundance_se_stack_%dkm.tif", res_km)),
                     overwrite = TRUE)
}

# ── Update range_source in batch_nsw_log.csv ─────────────────────────────────
log <- read.csv("batch_nsw_log.csv", stringsAsFactors = FALSE)
idx <- which(log$common_name == SP)
if (length(idx) == 1L) {
  old <- log$range_source[idx]
  log$range_source[idx] <- range_source
  log$run_date[idx]     <- as.character(Sys.Date())
  write.csv(log, "batch_nsw_log.csv", row.names = FALSE)
  message(sprintf("\nbatch_nsw_log: %s range_source %s -> %s", SP, old, range_source))
}
message("\nDone. Now re-run analysis/modeling/rejoin_seasonality_abund.R to refresh the CSV.")
