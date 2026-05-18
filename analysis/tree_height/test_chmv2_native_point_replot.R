#!/usr/bin/env Rscript
# Re-plot the cached CHMv2 tile with a colour scheme that emphasises the 0-1 m
# range (which is where most of this arid-zone window sits).

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(tidyterra)
})

PT_LON <- 141.704623
PT_LAT <- -31.087603
IN_TIF  <- "analysis/chmv2_native_point.tif"
OUT_PNG <- "analysis/chmv2_native_point_lowrange.png"

chm <- rast(IN_TIF)
# Classify into discrete height bins so each bin can be a distinct colour
# (scale_fill_stepsn would interpolate a continuous palette across breaks
# and ignore the bin-specific colours we want).
breaks_classify <- c(-Inf, 0.25, 0.5, 0.75, 1, 2, 5, 10, Inf)
chm_bin <- classify(chm, rcl = cbind(breaks_classify[-length(breaks_classify)],
                                     breaks_classify[-1],
                                     seq_along(breaks_classify[-1])))

# Hand-built palette: 0 = white, 0.25/0.5/0.75/1 m get distinct, *light*
# colours (so they pop against the background and from each other), then
# greens/yellows/oranges/reds up to 20 m. Discrete bins so adjacent low-end
# heights are obviously different categories rather than nearly-the-same hue.
bin_labels <- c("0–0.25", "0.25–0.5", "0.5–0.75", "0.75–1",
                "1–2", "2–5", "5–10", ">10")
pal <- c(
  "0–0.25"  = "#ffffff",  # white (true zero / no canopy)
  "0.25–0.5"  = "#fee391",  # pale yellow
  "0.5–0.75"  = "#fe9929",  # orange
  "0.75–1"  = "#cc4c02",  # dark red
  "1–2"  = "#a1d99b",  # light green
  "2–5"  = "#41ab5d",  # green
  "5–10" = "#225ea8",  # blue
  ">10"  = "#54278f"   # purple
)
# Make the raster a categorical SpatRaster so geom_spatraster maps it via
# scale_fill_manual rather than as a continuous gradient.
levels(chm_bin) <- data.frame(ID = seq_along(bin_labels), label = bin_labels)

p <- ggplot() +
  geom_spatraster(data = chm_bin) +
  scale_fill_manual(
    name = "Canopy height (m)",
    values = pal,
    breaks = bin_labels,
    drop = FALSE,
    na.value = "transparent"
  ) +
  geom_point(aes(x = PT_LON, y = PT_LAT), colour = "black", size = 2.5, shape = 4) +
  coord_sf(crs = 4326) +
  labs(
    title = sprintf("Meta CHMv2 (DINOv3) native resolution @ (%.4f, %.4f)",
                    PT_LAT, PT_LON),
    subtitle = "0–1 m: white→yellow→orange→red (each 0.25 m a distinct hue); 1 m+: green→blue→purple",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = "grey92", colour = NA),
        panel.grid = element_line(colour = "grey80"))

ggsave(OUT_PNG, p, width = 7.5, height = 7, dpi = 150)
message(sprintf("Wrote %s", OUT_PNG))
