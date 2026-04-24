suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(data.table)
})

EBD        <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"
SAMP       <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"
CACHE      <- "ebirdabund_cache"
OUTPUT_DIR <- "species_maps"
LOG_FILE   <- "batch_nsw_log.csv"
BOTW_PATH  <- "botw_species/BOTW_2025.gpkg"
TAXONOMY   <- "nsw_ebird_taxonomy.csv"

# ── Species list ──────────────────────────────────────────────────────────────
species_df   <- read.csv("nsw_species_list.csv", stringsAsFactors = FALSE)
species_list <- species_df$common_name
message(sprintf("Species to process: %d", length(species_list)))

# ── Pre-cache: build missing zerofilled .rds files in one EBD pass ───────────
# Prevents parallel workers competing to fread the large EBD simultaneously.
safe_name_local <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
cached <- sub("^zerofilled_", "", sub("[.]rds$", "",
              list.files(CACHE, pattern = "zerofilled_.*[.]rds")))
needs_cache <- species_list[!safe_name_local(species_list) %in% cached]

# Also skip species whose output already exists (already done in a prior run)
already_done <- sub("[.]tif$", "",
                    list.files(file.path(OUTPUT_DIR, "3km"),
                               pattern = "[.]tif$"))
species_list <- species_list[!safe_name_local(species_list) %in% already_done]
message(sprintf("Species already completed (skipping): %d",
                length(species_df$common_name) - length(species_list)))
needs_cache  <- needs_cache[!safe_name_local(needs_cache) %in% already_done]

if (length(needs_cache) > 0) {
  message(sprintf(
    "\n── Pre-caching %d species with missing EBD cache (one EBD pass) ──",
    length(needs_cache)
  ))

  message("  Reading sampling file...")
  samp_all <- as.data.frame(fread(
    SAMP, sep = "\t", quote = "", showProgress = FALSE, na.strings = "",
    select = c("SAMPLING EVENT IDENTIFIER", "ALL SPECIES REPORTED",
               "LATITUDE", "LONGITUDE",
               "OBSERVATION DATE", "TIME OBSERVATIONS STARTED",
               "DURATION MINUTES", "EFFORT DISTANCE KM", "NUMBER OBSERVERS",
               "PROTOCOL NAME")
  ))
  samp_all <- samp_all[samp_all[["ALL SPECIES REPORTED"]] == 1, ]

  message("  Reading EBD (full scan for uncached species)...")
  ebd_all <- as.data.frame(fread(
    EBD, select = c(6L, 11L, 35L),
    sep = "\t", quote = "", showProgress = TRUE, na.strings = ""
  ))
  names(ebd_all) <- c("common_name", "observation_count", "checklist_id")
  ebd_all <- ebd_all[ebd_all[["common_name"]] %in% needs_cache, ]

  nsw_bbox <- as.numeric(sf::st_bbox(sf::st_transform(
    sf::st_as_sf(geodata::gadm("AUS", 1, path = CACHE)[
      geodata::gadm("AUS", 1, path = CACHE)$NAME_1 == "New South Wales", ]),
    4326
  )))

  samp_bbox <- samp_all[
    !is.na(samp_all[["LATITUDE"]]) &
    samp_all[["LATITUDE"]]  >= nsw_bbox[2] & samp_all[["LATITUDE"]]  <= nsw_bbox[4] &
    samp_all[["LONGITUDE"]] >= nsw_bbox[1] & samp_all[["LONGITUDE"]] <= nsw_bbox[3], ]
  names(samp_bbox)[names(samp_bbox) == "SAMPLING EVENT IDENTIFIER"] <- "checklist_id"
  names(samp_bbox)[names(samp_bbox) == "LATITUDE"]  <- "latitude"
  names(samp_bbox)[names(samp_bbox) == "LONGITUDE"] <- "longitude"
  names(samp_bbox)[names(samp_bbox) == "OBSERVATION DATE"] <- "observation_date"
  names(samp_bbox)[names(samp_bbox) == "TIME OBSERVATIONS STARTED"] <- "time_observations_started"
  names(samp_bbox)[names(samp_bbox) == "DURATION MINUTES"] <- "duration_minutes"
  names(samp_bbox)[names(samp_bbox) == "EFFORT DISTANCE KM"] <- "effort_distance_km"
  names(samp_bbox)[names(samp_bbox) == "NUMBER OBSERVERS"] <- "number_observers"
  names(samp_bbox)[names(samp_bbox) == "PROTOCOL NAME"] <- "protocol_type"

  message(sprintf("  Building caches for %d species...", length(needs_cache)))
  for (sp in needs_cache) {
    cache_f <- file.path(CACHE, sprintf("zerofilled_%s.rds", safe_name_local(sp)))
    if (file.exists(cache_f)) next
    ebd_sp  <- ebd_all[ebd_all[["common_name"]] == sp,
                        c("checklist_id", "observation_count"), drop = FALSE]
    zf <- merge(samp_bbox[, setdiff(names(samp_bbox), "ALL SPECIES REPORTED")],
                ebd_sp, by = "checklist_id", all.x = TRUE)
    zf[["observation_count"]][is.na(zf[["observation_count"]])] <- "0"
    zf[["species_observed"]] <- zf[["observation_count"]] != "0"
    saveRDS(zf, cache_f)
    message(sprintf("    Cached: %s", sp))
  }
  rm(samp_all, ebd_all, samp_bbox)
  gc()
  message("  Pre-caching complete.\n")
}

# ── NSW boundary ──────────────────────────────────────────────────────────────
message("Getting NSW boundary...")
aus <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])

# ── Covariates (built once, shared across all workers) ────────────────────────
message("Preparing covariates...")
cov <- prepare_covariates(nsw, cache_dir = CACHE)
message("Covariate layers: ", paste(names(cov), collapse = ", "))

# ── Taxonomy (common name -> scientific name for BOTW range lookup) ───────────
taxonomy <- read.csv(TAXONOMY, stringsAsFactors = FALSE)
message(sprintf("Taxonomy loaded: %d species", nrow(taxonomy)))

# ── Batch run ─────────────────────────────────────────────────────────────────
message(sprintf("\nStarting batch of %d species...", length(species_list)))
t_start <- proc.time()[["elapsed"]]

results <- estimate_abundance_batch(
  polygon      = nsw,
  ebird_zip    = EBD,
  sampling_txt = SAMP,
  species_list = species_list,
  taxonomy     = taxonomy,
  cov_stack    = cov,
  cache_dir    = CACHE,
  grid_res_km  = c(3, 9),
  botw_path    = BOTW_PATH,
  output_dir   = OUTPUT_DIR
)

t_total <- proc.time()[["elapsed"]] - t_start

# ── Write log ─────────────────────────────────────────────────────────────────
ok      <- !vapply(results, inherits, TRUE, "error")
log_df  <- data.frame(
  common_name   = names(results),
  status        = ifelse(ok, "ok", "failed"),
  error_message = vapply(results, function(r) {
    if (inherits(r, "error")) conditionMessage(r) else NA_character_
  }, character(1)),
  n_checklists  = vapply(results, function(r) {
    if (inherits(r, "error")) NA_integer_ else r$n_checklists
  }, integer(1)),
  stringsAsFactors = FALSE
)
# Append to existing log if present (supports resumed runs)
if (file.exists(LOG_FILE)) {
  prior <- read.csv(LOG_FILE, stringsAsFactors = FALSE)
  log_df <- rbind(prior[!prior$common_name %in% log_df$common_name, ], log_df)
}
write.csv(log_df, LOG_FILE, row.names = FALSE)

# ── Build per-resolution raster stacks ───────────────────────────────────────
message("\nBuilding species abundance stacks...")
for (res_km in c(3, 9)) {
  res_dir   <- file.path(OUTPUT_DIR, paste0(res_km, "km"))
  tif_files <- sort(list.files(res_dir, pattern = "\\.tif$", full.names = TRUE))

  if (length(tif_files) == 0) {
    message(sprintf("  No .tif files found in %s, skipping stack.", res_dir))
    next
  }

  sp_names <- sub("\\.tif$", "", basename(tif_files))
  message(sprintf("  %dkm: stacking %d species...", res_km, length(tif_files)))

  abd_layers <- lapply(tif_files, function(f) terra::rast(f)[["abd"]])
  stack      <- terra::rast(abd_layers)
  names(stack) <- sp_names

  stack_path <- file.path(OUTPUT_DIR, sprintf("nsw_abundance_stack_%dkm.tif", res_km))
  terra::writeRaster(stack, stack_path, overwrite = TRUE)
  message(sprintf("  Saved: %s  (%d bands)", stack_path, terra::nlyr(stack)))
}

# ── Summary ───────────────────────────────────────────────────────────────────
hrs  <- floor(t_total / 3600)
mins <- floor((t_total %% 3600) / 60)
secs <- round(t_total %% 60)
cat(sprintf(
  "\n════════════════════════════════════════\n%d/%d species succeeded  |  %dh %02dm %02ds\nFailed: %s\nLog: %s\n",
  sum(ok), length(results), hrs, mins, secs,
  if (any(!ok)) paste(names(results)[!ok], collapse = ", ") else "none",
  LOG_FILE
))
