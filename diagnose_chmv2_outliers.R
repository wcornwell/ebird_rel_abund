suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(terra)
  library(ggplot2)
  library(dplyr)
})

CACHE <- "ebirdabund_cache"
bb    <- c(150, -34, 151, -33)
ext   <- terra::ext(bb[1], bb[3], bb[2], bb[4])

# Pull the cached 38 m raster (already 0–100 m clamped). We'll inspect
# values >40 m, which is the 99.9th-percentile and above.
chm <- terra::rast(file.path(
  CACHE, sprintf("meta_chmv2_%.2f_%.2f_%.2f_%.2f.tif", bb[1], bb[2], bb[3], bb[4])
))

# Outlier mask
THRESH <- 40
outliers_mask <- chm > THRESH
n_outliers <- terra::global(outliers_mask, "sum", na.rm = TRUE)[[1]]
n_total    <- terra::global(!is.na(chm), "sum")[[1]]
message(sprintf("Outliers (>%.0f m): %d of %d pixels (%.4f%%)",
                THRESH, n_outliers, n_total, 100 * n_outliers / n_total))

# Pull SRTM elevation for the same bbox and compute slope (degrees)
elev_raw <- geodata::elevation_30s(country = "AUS", path = CACHE)
elev     <- terra::crop(elev_raw, ext)
# Reproject elev to match chm grid for slope-at-pixel join
elev_r   <- terra::resample(elev, chm)
slope    <- terra::terrain(elev_r, v = "slope", unit = "degrees")

# Extract outlier pixels with their elevation + slope
outlier_pts   <- terra::as.points(terra::mask(chm, outliers_mask, maskvalue = 0))
outlier_coord <- terra::crds(outlier_pts)
outlier_chm   <- terra::values(outlier_pts)[, 1]
outlier_elev  <- terra::extract(elev_r, outlier_coord)[, 1]
outlier_slope <- terra::extract(slope,  outlier_coord)[, 1]

outliers_df <- data.frame(
  x = outlier_coord[, 1], y = outlier_coord[, 2],
  height = outlier_chm, elev = outlier_elev, slope = outlier_slope
) |> arrange(desc(height))

message("\n--- Top 20 outlier pixels ---")
print(head(outliers_df, 20))

message("\n--- Slope distribution at outlier locations vs. background ---")
bg_slope <- terra::values(slope)[, 1]
bg_slope <- bg_slope[!is.na(bg_slope)]
message(sprintf("  background slope: median=%.1f°, p95=%.1f°, max=%.1f°",
                median(bg_slope), quantile(bg_slope, 0.95), max(bg_slope)))
message(sprintf("  outlier slope:    median=%.1f°, p95=%.1f°, max=%.1f°",
                median(outliers_df$slope, na.rm = TRUE),
                quantile(outliers_df$slope, 0.95, na.rm = TRUE),
                max(outliers_df$slope, na.rm = TRUE)))
message(sprintf("  fraction outliers on slope >20°: %.1f%%",
                100 * mean(outliers_df$slope > 20, na.rm = TRUE)))

# Build a hillshade for a topo backdrop
slope_rad <- terra::terrain(elev_r, "slope", unit = "radians")
aspect    <- terra::terrain(elev_r, "aspect", unit = "radians")
hillshade <- terra::shade(slope_rad, aspect, angle = 35, direction = 315)

# Render: hillshade base + tree height + red outliers
hill_df <- as.data.frame(hillshade, xy = TRUE)
names(hill_df) <- c("x", "y", "shade")

chm_df <- as.data.frame(chm, xy = TRUE)
names(chm_df) <- c("x", "y", "tree_height")
chm_df <- chm_df[!is.na(chm_df$tree_height), ]

p <- ggplot() +
  geom_raster(data = hill_df, aes(x = x, y = y, fill = shade),
              alpha = 1, show.legend = FALSE) +
  scale_fill_gradient(low = "grey20", high = "white") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = chm_df, aes(x = x, y = y, fill = tree_height),
              alpha = 0.65) +
  scale_fill_gradientn(
    colours = c("white", "#cde9c0", "#5ab566", "#1a7337", "#003d12"),
    name = "Tree height (m)",
    limits = c(0, 35), oob = scales::squish
  ) +
  geom_point(data = outliers_df,
             aes(x = x, y = y, colour = slope),
             size = 1.0, shape = 16) +
  scale_colour_gradient(low = "#ffaa00", high = "red",
                        name = "Slope at\noutlier (°)") +
  coord_equal(expand = FALSE) +
  labs(title = sprintf("CHMv2 outliers >%.0f m, on hillshade backdrop", THRESH),
       subtitle = sprintf("%d outliers of %d pixels (%.4f%%) — colour = local slope",
                          n_outliers, n_total, 100 * n_outliers / n_total)) +
  theme_minimal(base_size = 11)

ggsave("diagnose_chmv2_outliers.png", p, width = 9, height = 8, dpi = 150)
message("\nSaved diagnose_chmv2_outliers.png")
