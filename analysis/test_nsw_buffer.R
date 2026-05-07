suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(data.table)
})

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
OUTPUT_DIR     <- "species_maps_test"
LOG_FILE       <- "test_nsw_buffer_log.csv"
BOTW_PATH      <- "botw_species/BOTW_2025.gpkg"
TAXONOMY       <- "nsw_ebird_taxonomy.csv"
GRID_RES_KM    <- c(3, 9)
SPECIES_LIST   <- c("Gray Currawong", "Australian Logrunner")

dir.create(ZEROFILL_CACHE, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR,     showWarnings = FALSE, recursive = TRUE)

message("Getting NSW boundary (+ 100 km buffer)...")
aus     <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw     <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(nsw, 3577), 100000),
  4326
)
study_bbox <- as.numeric(sf::st_bbox(polygon))

safe_name_local <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# ── Pre-cache ─────────────────────────────────────────────────────────────────
cached      <- sub("^zerofilled_", "", sub("[.]rds$", "",
               list.files(ZEROFILL_CACHE, pattern = "zerofilled_.*[.]rds")))
needs_cache <- SPECIES_LIST[!safe_name_local(SPECIES_LIST) %in% cached]

if (length(needs_cache) > 0) {
  message(sprintf("\n── Pre-caching %d species ──", length(needs_cache)))
  sampling_df <- read_sampling(SAMP, study_bbox)
  valid_ids   <- sampling_df$checklist_id
  message(sprintf("  Scanning %d EBD file(s)...", length(EBD)))
  ebd_all <- read_ebd_observations(EBD,
                                   species_set         = needs_cache,
                                   valid_checklist_ids = valid_ids)
  for (i in seq_along(needs_cache)) {
    sp      <- needs_cache[i]
    cache_f <- file.path(ZEROFILL_CACHE, sprintf("zerofilled_%s.rds", safe_name_local(sp)))
    if (file.exists(cache_f)) next
    ebd_sp  <- ebd_all[ebd_all[["common_name"]] == sp,
                       c("checklist_id", "observation_count"), drop = FALSE]
    saveRDS(zero_fill(sampling_df, ebd_sp), cache_f)
    message(sprintf("    Cached (%d/%d): %s", i, length(needs_cache), sp))
  }
  rm(sampling_df, ebd_all); gc()
  message("  Pre-caching complete.\n")
}

# ── Covariates ────────────────────────────────────────────────────────────────
message("Preparing covariates...")
cov <- prepare_covariates(polygon, cache_dir = CACHE)

# ── Sampling master ───────────────────────────────────────────────────────────
prepare_sampling_master(SAMP, polygon, cov, ZEROFILL_CACHE)

# ── Taxonomy ──────────────────────────────────────────────────────────────────
taxonomy <- read.csv(TAXONOMY, stringsAsFactors = FALSE)

# ── Batch run ─────────────────────────────────────────────────────────────────
message(sprintf("\nRunning batch for: %s", paste(SPECIES_LIST, collapse = ", ")))
results <- estimate_abundance_batch(
  polygon      = polygon,
  ebird_zip    = EBD,
  sampling_txt = SAMP,
  species_list = SPECIES_LIST,
  taxonomy     = taxonomy,
  cov_stack    = cov,
  cache_dir    = ZEROFILL_CACHE,
  grid_res_km  = GRID_RES_KM,
  botw_path    = BOTW_PATH,
  border       = nsw,
  output_dir   = OUTPUT_DIR,
  log_file     = LOG_FILE
)

# ── Summary ───────────────────────────────────────────────────────────────────
log_df <- read.csv(LOG_FILE, stringsAsFactors = FALSE)
print(log_df[, c("common_name", "status", "n_checklists", "n_positive", "dev_expl")])
