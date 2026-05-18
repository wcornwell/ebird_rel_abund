#!/usr/bin/env Rscript
# Log-scale histogram of CHMv2 native pixel heights for the 1x1 km window.

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
})

IN_TIF  <- "analysis/chmv2_native_point.tif"
OUT_PNG <- "analysis/chmv2_native_point_hist.png"

chm <- rast(IN_TIF)
v   <- values(chm)
v   <- v[!is.na(v)]
n_total <- length(v)
n_zero  <- sum(v == 0)
n_pos   <- sum(v > 0)

message(sprintf("Pixels: %d total, %d == 0 (%.1f%%), %d > 0 (%.1f%%)",
                n_total, n_zero, 100 * n_zero / n_total,
                n_pos,  100 * n_pos  / n_total))
message(sprintf("Non-zero: min %.3f, median %.2f, p90 %.2f, p99 %.2f, max %.2f",
                min(v[v > 0]), median(v[v > 0]),
                quantile(v[v > 0], 0.9), quantile(v[v > 0], 0.99), max(v)))

# Two-panel layout: a single 0-bar (linear count axis would dwarf everything;
# we use log10(count) y-axis throughout so the long tail is visible too) and
# the >0 distribution on a log x-axis.
df_zero <- data.frame(panel = "= 0 m", height = 0, count = n_zero)
df_pos  <- data.frame(height = v[v > 0])

# Log-spaced breaks across the observed range.
brks_x <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20)
brks_minor <- exp(seq(log(0.05), log(20), length.out = 50))

p_zero <- ggplot(df_zero, aes(x = factor(panel), y = count)) +
  geom_col(width = 0.6, fill = "grey40") +
  scale_y_log10(labels = scales::comma) +
  labs(x = NULL, y = "Pixel count (log)", title = "Exactly 0 m") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 10))

p_pos <- ggplot(df_pos, aes(x = height)) +
  geom_histogram(breaks = brks_minor, fill = "steelblue", colour = "white",
                 linewidth = 0.1) +
  scale_x_log10(breaks = brks_x, labels = brks_x) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Canopy height (m, log scale)",
       y = "Pixel count (log)",
       title = sprintf("> 0 m  (n = %s pixels)", format(n_pos, big.mark = ","))) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 10))

# Combine via patchwork if available, else save the >0 panel alone.
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p <- (p_zero + p_pos) +
    plot_layout(widths = c(1, 5)) +
    plot_annotation(
      title = "Meta CHMv2 (DINOv3) native pixel heights",
      subtitle = sprintf("1 km x 1 km window at (-31.0876, 141.7046); %s pixels",
                         format(n_total, big.mark = ","))
    )
} else {
  p <- p_pos +
    labs(title = "Meta CHMv2 (DINOv3) heights (> 0 m), 1 km x 1 km window")
}

ggsave(OUT_PNG, p, width = 9, height = 4.5, dpi = 150)
message(sprintf("Wrote %s", OUT_PNG))
