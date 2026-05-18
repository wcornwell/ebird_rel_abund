suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(ggplot2)
})

EBD   <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"
SAMP  <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"
CACHE <- "ebirdabund_cache"

aus <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])

cov <- prepare_covariates(nsw, cache_dir = CACHE)

for (sp in c("White-browed Woodswallow", "Plumed Whistling-Duck",
             "Tawny Frogmouth", "Powerful Owl")) {
  message("\n\n====== ", sp, " ======")
  fit <- fit_species_model(
    polygon      = nsw,
    ebird_zip    = EBD,
    sampling_txt = SAMP,
    species      = sp,
    cov_stack    = cov,
    cache_dir    = CACHE
  )
  safe <- gsub("[^a-z0-9]+", "_", tolower(sp))

  p <- plot_gam_smooths(fit$model, sp)
  ggsave(paste0("test_filter_", safe, ".png"), p, width = 16, height = 12, dpi = 150)
  message("Saved: test_filter_", safe, ".png")

  pred <- predict_species_map(fit, polygon = nsw, species = sp, use_range = FALSE)
  ggsave(paste0("test_filter_map_", safe, ".png"), pred$plot, width = 10, height = 10, dpi = 150)
  message("Saved: test_filter_map_", safe, ".png")
}

message("\nDone.")
