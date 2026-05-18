# Concrete uses of test_covariate_variant() — assess alternative tree-height
# layers / aggregations against the current CHMv2 baseline.
#
# Variant 1 (two-stage aggregation): re-aggregates the cached meta_chmv2 raster
# via (~38 m → ~250 m mean → ~1 km p90) instead of the current single-step
# (~38 m → ~1 km p90). This is "option C" from CLAUDE.md's discussion of CHMv2
# stripe artefacts: smoothing within stripes before the p90 step.
#
# Variant 2 (GLAD/Potapov 2019): Landsat + GEDI canopy height at 30 m
# (Potapov et al. 2021, Remote Sensing of Environment). Independent of
# CHMv2's Sentinel-2 lineage, so any orbit-stripe artefacts won't correlate.
# Download: https://glad.umd.edu/dataset/gedi/ (1° tiles, manually mosaicked
# to a NSW-extent GeoTIFF). Set GLAD_TREE_HEIGHT_TIF to point at it.
#
# Variant 3 (ETH-Lang 2020): deliberately skipped here — Google Earth Engine
# is the primary distribution channel and there's no clean HTTP/COG endpoint.
# To add it later, follow the same pattern as the GLAD block.
#
# Each variant is a separate two-way paired CV against the CHMv2 baseline.
# Compare the resulting per_species_metrics.csv files to read the three-way.
#
# Runtime: ~10–20 min per variant on the 19-species test set. Reuses
# sampling_master.rds and zerofilled_*.rds — never invalidates either cache.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(yaml)
  library(sf)
  library(terra)
  library(geodata)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Config + polygon (same as run_batch_nsw.R) ──────────────────────────────
CFG_PATH <- Sys.getenv("EBIRD_CONFIG", "config.yaml")
cfg <- yaml::read_yaml(CFG_PATH)

CACHE          <- cfg$cache_dir %||% "ebirdabund_cache"
ZEROFILL_CACHE <- cfg$zerofill_cache %||% "ebirdabund_cache_nsw_buffer"

message("\n── Building polygon + covariate stack ───────────────────────────")
aus <- geodata::gadm(country = cfg$study_polygon$country, level = 1,
                     path = CACHE)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_sf, cfg$study_polygon$proj_processing),
                cfg$study_polygon$buffer_metres),
  cfg$study_polygon$proj_output
)
cov_stack <- prepare_covariates(polygon, cache_dir = CACHE)

# ── Variant 1: two-stage CHMv2 aggregation ──────────────────────────────────
# Reads the existing cached meta_chmv2_*.tif (no re-download), mean-aggregates
# to ~250 m, then p90-aggregates to the template ~1 km. Returns a single-layer
# SpatRaster named "tree_height" so apply_variant() can swap it in cleanly.

two_stage_chmv2 <- function(bbox, template) {
  cache_f <- file.path(
    CACHE,
    sprintf("meta_chmv2_%.2f_%.2f_%.2f_%.2f.tif",
            bbox[1], bbox[2], bbox[3], bbox[4])
  )
  if (!file.exists(cache_f)) {
    stop("Cached CHMv2 raster not found: ", cache_f,
         "\nRun prepare_covariates() first to materialise it.")
  }
  raw <- terra::rast(cache_f)
  raw_crop <- terra::crop(raw, terra::ext(template))

  # Stage 1: ~38 m → ~250 m via mean (smooths within-stripe heterogeneity).
  fact1 <- max(1L, round(250 / 38))
  m250  <- terra::aggregate(raw_crop, fact = fact1, fun = "mean",
                            na.rm = TRUE)

  # Stage 2: 250 m → template (~1 km) via p90.
  fact2 <- max(1L, round(terra::res(template)[1] / terra::res(m250)[1]))
  p90_1km <- terra::aggregate(
    m250, fact = fact2,
    fun = function(x, ...) stats::quantile(x, 0.9, na.rm = TRUE)
  )
  out <- terra::resample(p90_1km, template[[1]], method = "bilinear")
  names(out) <- "tree_height"
  out
}

# Snap bbox to the same precision used by the existing cache filename so we
# hit the right cached raster. The cov_stack and polygon already agree on
# extent, so just read bbox from the polygon.
polygon_bbox <- as.numeric(sf::st_bbox(sf::st_transform(polygon, 4326)))
# Quick existence check up front for a clearer error than waiting until
# extract time.
expected_chmv2 <- file.path(
  CACHE,
  sprintf("meta_chmv2_%.2f_%.2f_%.2f_%.2f.tif",
          polygon_bbox[1], polygon_bbox[2],
          polygon_bbox[3], polygon_bbox[4])
)
if (!file.exists(expected_chmv2))
  stop("Expected cached CHMv2 not found at: ", expected_chmv2,
       "\nThe two-stage variant cannot run until prepare_covariates() has ",
       "produced it. Trigger that by running the main batch script once.")

message("\n── Running variant: tree_height_two_stage_p90 ──────────────────")
res_two_stage <- test_covariate_variant(
  variant_name   = "tree_height_two_stage_p90",
  variant_source = two_stage_chmv2,
  change         = list(swap = "tree_height"),
  cache_dir      = ZEROFILL_CACHE,
  polygon        = polygon,
  cov_stack      = cov_stack
)

cat("\n── Summary (variant - baseline; positive = better unless RMSE/MAE) ──\n")
print(res_two_stage$summary)


# ── Variant 2: GLAD/Potapov 2019 (Landsat + GEDI, 30 m) ─────────────────────
# Set GLAD_TREE_HEIGHT_TIF to a local mosaic covering the polygon.
# GLAD's GEDI Forest Height product is byte-encoded: 0–100 = height in m,
# values >100 are non-height codes (no-data, water, snow, etc.). We set those
# to NA so they don't contaminate the aggregation.
#
# To keep the comparison apples-to-apples, we apply the same aggregation
# pipeline CHMv2 uses in covariates.R::load_meta_chmv2: p90 over a window
# sized to match the template (~1 km), then resample. This isolates the
# **data source** from any resolution effect.

glad_path <- Sys.getenv("GLAD_TREE_HEIGHT_TIF",
                        "data/glad_potapov_2019_nsw.tif")

glad_source <- function(bbox, template) {
  if (!file.exists(glad_path))
    stop("GLAD raster not found at ", glad_path,
         ". Set GLAD_TREE_HEIGHT_TIF or download from ",
         "https://glad.umd.edu/dataset/gedi/")
  r <- terra::rast(glad_path)
  ext_t <- terra::ext(terra::project(
    terra::vect(terra::ext(template), crs = terra::crs(template)),
    terra::crs(r)
  ))
  r <- terra::crop(r, ext_t)
  r <- terra::ifel(r > 100, NA, r)

  # Same p90 aggregation CHMv2 uses, but starting from GLAD's native ~30 m.
  ratio <- max(1L, round(terra::res(template)[1] / terra::res(r)[1]))
  agg <- terra::aggregate(
    r, fact = ratio,
    fun = function(x, ...) stats::quantile(x, 0.9, na.rm = TRUE)
  )
  out <- terra::resample(agg, template[[1]], method = "bilinear")
  names(out) <- "tree_height"
  out
}

if (file.exists(glad_path)) {
  message("\n── Running variant: tree_height_glad_potapov ────────────────────")
  res_glad <- test_covariate_variant(
    variant_name   = "tree_height_glad_potapov",
    variant_source = glad_source,
    change         = list(swap = "tree_height"),
    cache_dir      = ZEROFILL_CACHE,
    polygon        = polygon,
    cov_stack      = cov_stack
  )
  cat("\n── GLAD summary ──\n"); print(res_glad$summary)
} else {
  message("\n(skipping GLAD variant — raster not found at ", glad_path,
          ".\n Set GLAD_TREE_HEIGHT_TIF or place the mosaic at that path.)")
}


# ── Three-way comparison report (after both variants have run) ──────────────
# Each variant lives in its own subdirectory. Read both per_species_metrics.csv
# files, align on species, and print the deltas side by side.

ts_csv   <- "covariate_tests/tree_height_two_stage_p90/per_species_metrics.csv"
glad_csv <- "covariate_tests/tree_height_glad_potapov/per_species_metrics.csv"

if (file.exists(ts_csv) && file.exists(glad_csv)) {
  ts   <- utils::read.csv(ts_csv,   stringsAsFactors = FALSE)
  glad <- utils::read.csv(glad_csv, stringsAsFactors = FALSE)
  keep <- c("species", "spearman_r_base", "spearman_r_delta", "rmse_delta",
            "holdout_dev_expl_delta")
  ts2   <- ts[,   keep]
  glad2 <- glad[, keep]
  names(ts2)[-1]   <- paste0(names(ts2)[-1],   "_twostage")
  names(glad2)[-1] <- paste0(names(glad2)[-1], "_glad")
  combined <- merge(ts2, glad2, by = "species", all = TRUE)
  combined <- combined[order(-combined$spearman_r_base_twostage), ]
  out_f <- "covariate_tests/tree_height_three_way.csv"
  utils::write.csv(combined, out_f, row.names = FALSE)
  message("\nWrote ", out_f)
  message("\nPer-species deltas vs baseline CHMv2:")
  print(combined, row.names = FALSE)
}
