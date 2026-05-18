#!/usr/bin/env Rscript
# Smoke test: does load_palsar_hv() work on a 1°x1° bbox? Pulls a single tile,
# aggregates to a 1 km template, validates the output stats.

suppressPackageStartupMessages({
  library(terra)
})
source("analysis/palsar_helpers.R")

# Test bbox: 1°x1° around the existing CHMv2 test point (Western NSW mallee/scrub).
bb  <- c(141.2, -31.5, 142.2, -30.5)
ext <- terra::ext(bb[1], bb[3], bb[2], bb[4])
tmpl <- terra::rast(ext, resolution = 1000 / 111320, crs = "EPSG:4326")

cache_dir <- "ebirdabund_cache"
dir.create(cache_dir, showWarnings = FALSE)

t0 <- Sys.time()
hv <- load_palsar_hv(bb, ext, tmpl, cache_dir, year = 2020L)
message(sprintf("Done in %.1fs", as.numeric(Sys.time() - t0, units = "secs")))

message(sprintf("Output: %d x %d, res %.5f deg", nrow(hv), ncol(hv), res(hv)[1]))
v <- values(hv); v <- v[!is.na(v)]
message(sprintf(
  "HV gamma0: min %.1f, p10 %.1f, median %.1f, p90 %.1f, max %.1f dB (n=%d, NA=%d)",
  min(v), quantile(v, 0.1), median(v), quantile(v, 0.9),
  max(v), length(v), ncell(hv) - length(v)
))

# Quick PNG to eyeball.
suppressPackageStartupMessages(library(ggplot2)); library(tidyterra)
p <- ggplot() + geom_spatraster(data = hv) +
  scale_fill_viridis_c(name = "HV γ° (dB)", option = "magma",
                       na.value = "transparent") +
  coord_sf(crs = 4326) +
  labs(title = sprintf("PALSAR HV smoke test (%d x %d cells, ~1 km)",
                       nrow(hv), ncol(hv))) +
  theme_minimal(base_size = 11)
ggsave("analysis/palsar_helper_smoke.png", p, width = 7, height = 6, dpi = 130)
message("Wrote analysis/palsar_helper_smoke.png")
