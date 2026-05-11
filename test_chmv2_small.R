suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(terra)
  library(sf)
  library(ggplot2)
})

CACHE <- "ebirdabund_cache"
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

# Small bbox over the Blue Mountains / western Sydney — mix of dense eucalypt
# canopy (Wollemi, Blue Mtns NPs) and cleared land. Easy visual check.
bb  <- c(150, -34, 151, -33)  # ~111 km × 111 km
ext <- terra::ext(bb[1], bb[3], bb[2], bb[4])

# Dummy ~1 km template (matches the resolution prepare_covariates resamples to).
tmpl <- terra::rast(ext, resolution = 1 / 111.32, crs = "EPSG:4326")

message("--- Testing load_meta_chmv2() on a small bbox ---")
message(sprintf("  bbox: lon [%g, %g], lat [%g, %g]", bb[1], bb[3], bb[2], bb[4]))
message(sprintf("  template: %s pixels at ~1 km", paste(dim(tmpl)[1:2], collapse = "x")))

t0 <- Sys.time()
r <- load_meta_chmv2(bb, ext, tmpl, CACHE)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")

message(sprintf("\n--- Done in %.1f sec ---", elapsed))
message(sprintf("  output dim: %s", paste(dim(r)[1:2], collapse = "x")))
message(sprintf("  value range: [%.1f, %.1f]",
                terra::global(r, "min", na.rm = TRUE)[[1]],
                terra::global(r, "max", na.rm = TRUE)[[1]]))
message(sprintf("  fraction NA: %.2f%%",
                100 * terra::global(r, fun = "isNA")[[1]] /
                  terra::ncell(r)))

# Also load the raw 38 m cache for a sharper map
cache_file <- file.path(
  CACHE, sprintf("meta_chmv2_%.2f_%.2f_%.2f_%.2f.tif", bb[1], bb[2], bb[3], bb[4])
)
raw <- terra::rast(cache_file)
message(sprintf("  cache file: %s (%.1f MB)",
                basename(cache_file), file.size(cache_file) / 1e6))
message(sprintf("  cache dim:  %s at %.5f deg (%.0f m)",
                paste(dim(raw)[1:2], collapse = "x"),
                terra::res(raw)[1], terra::res(raw)[1] * 111320))

# Quick map — show the raw 38 m for detail
df <- as.data.frame(raw, xy = TRUE)
names(df) <- c("x", "y", "tree_height")
df <- df[!is.na(df$tree_height), ]

p <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = tree_height)) +
  scale_fill_gradientn(
    colours = c("white", "#cde9c0", "#5ab566", "#1a7337", "#003d12"),
    name = "Tree height (m)",
    limits = c(0, 35), oob = scales::squish
  ) +
  coord_equal(expand = FALSE) +
  labs(title = "Meta/WRI CHMv2 — Blue Mountains / W. Sydney test",
       subtitle = sprintf("OL=4 ≈ 38 m, %d × %d cache pixels",
                          dim(raw)[2], dim(raw)[1])) +
  theme_minimal(base_size = 11)

ggsave("test_chmv2_small.png", p, width = 8, height = 7, dpi = 150)
message("\nSaved test_chmv2_small.png")
