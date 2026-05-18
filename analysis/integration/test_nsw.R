suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
})

EBD     <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"
SAMP    <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"
CACHE   <- "ebirdabund_cache"
SPECIES <- "Superb Parrot"

# NSW state boundary
message("Getting NSW boundary...")
aus <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])

# ── Stage 1: prepare covariates once ─────────────────────────────────────────
message("\n=== Stage 1: Prepare covariates ===")
cov <- prepare_covariates(nsw, cache_dir = CACHE)
message("Covariate layers: ", paste(names(cov), collapse = ", "))

# ── Stage 2: estimate abundance ───────────────────────────────────────────────
message("\n=== Stage 2: Estimate abundance ===")
result <- estimate_abundance(
  polygon      = nsw,
  ebird_zip    = EBD,
  sampling_txt = SAMP,
  species      = SPECIES,
  cov_stack    = cov,
  cache_dir    = CACHE,
  grid_res_km  = 3
)

message("Saving outputs...")
ggplot2::ggsave(
  "sparrot_nsw.png",
  result$plot,
  width = 10, height = 9, dpi = 150
)
terra::writeRaster(
  result$predictions, "White-throated Needletail_nsw.tif", overwrite = TRUE
)

cat("\n=== Done ===\n")
cat("Checklists used:", nrow(result$data), "\n")
cat("Deviance explained:",
    round(summary(result$model)$dev.expl * 100, 1), "%\n")
cat("abd range: [",
    round(terra::minmax(result$predictions[["abd"]]), 4), "]\n")
cat("Outputs: superb_fairywren_nsw.png / .tif\n")
