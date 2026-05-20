# Test PALSAR HV gamma0 (L-band SAR backscatter) as a new habitat covariate,
# using the test_covariate_variant() framework so results sit alongside the
# tree_height and nightlights variant tests in covariate_tests/.
#
# L-band SAR (ALOS-2 PALSAR-2, ~24 cm wavelength) penetrates the canopy and
# responds to woody biomass / volume scattering. Complements CHMv2 (height)
# and ESA WorldCover (cover fractions) — two woodlands at the same height
# and cover can still differ in stem density / biomass, which HV gamma0
# should pick up. Also no Sentinel-2 orbit-stripe artefacts (it's radar).
#
# Reuses sampling_master.rds and zerofilled_*.rds — never invalidates either.
# The PALSAR raster is cached to CACHE/palsar_hv_*.tif (via load_palsar_hv),
# and the per-checklist extract is cached to
# ebirdabund_cache_nsw_buffer/variant_extracts/palsar_hv_2020.rds.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(yaml)
  library(sf)
  library(terra)
  library(geodata)
})

source("analysis/palsar/palsar_helpers.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Config + polygon (same as run_batch_nsw.R) ──────────────────────────────
CFG_PATH <- Sys.getenv("EBIRD_CONFIG", "config.yaml")
cfg <- yaml::read_yaml(CFG_PATH)

CACHE          <- cfg$cache_dir %||% "ebirdabund_cache"
ZEROFILL_CACHE <- cfg$zerofill_cache %||% "ebirdabund_cache_nsw_buffer"
PALSAR_YEAR    <- 2020L

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

# ── PALSAR variant source (function form) ───────────────────────────────────
# Returns a single-layer SpatRaster named "palsar_hv" so apply_variant() can
# add it cleanly. test_covariate_variant() calls this with the polygon bbox
# in EPSG:4326 and the existing cov_stack as template.
palsar_hv_source <- function(bbox, template) {
  ext_buf <- terra::ext(bbox[1], bbox[3], bbox[2], bbox[4])
  r <- load_palsar_hv(
    bb        = bbox,
    ext       = ext_buf,
    template  = template[[1]],
    cache_dir = CACHE,
    year      = PALSAR_YEAR
  )
  names(r) <- "palsar_hv"
  r
}

# ── Run the variant test ────────────────────────────────────────────────────
message("\n── Running variant: palsar_hv_2020 ─────────────────────────────")
res <- test_covariate_variant(
  variant_name   = sprintf("palsar_hv_%d", PALSAR_YEAR),
  variant_source = palsar_hv_source,
  change         = list(add = "palsar_hv"),
  cache_dir      = ZEROFILL_CACHE,
  polygon        = polygon,
  cov_stack      = cov_stack
)

cat("\n── Summary (variant − baseline; positive = better unless RMSE/MAE) ──\n")
print(res$summary)
