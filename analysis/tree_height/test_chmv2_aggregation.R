suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
})

CACHE <- "ebirdabund_cache"
bb    <- c(150, -34, 151, -33)
ext   <- terra::ext(bb[1], bb[3], bb[2], bb[4])

cache_file <- file.path(
  CACHE, sprintf("meta_chmv2_%.2f_%.2f_%.2f_%.2f.tif", bb[1], bb[2], bb[3], bb[4])
)
raw <- terra::rast(cache_file)
message(sprintf("Cache: %s pixels at %.0f m  range [%.1f, %.1f]",
                paste(dim(raw)[1:2], collapse = "x"),
                terra::res(raw)[1] * 111320,
                terra::global(raw, "min", na.rm = TRUE)[[1]],
                terra::global(raw, "max", na.rm = TRUE)[[1]]))

# Master ~1 km template — what the model actually sees.
tmpl <- terra::rast(ext, resolution = 1 / 111.32, crs = "EPSG:4326")
ratio <- round(terra::res(tmpl)[1] / terra::res(raw)[1])
message(sprintf("Template: %s at ~1 km  (aggregation factor ~%d)",
                paste(dim(tmpl)[1:2], collapse = "x"), ratio))

# --- Three aggregation methods ---
message("\nComputing mean, max, p90 aggregations...")

t0 <- Sys.time()
agg_mean <- terra::resample(raw, tmpl, method = "average")
message(sprintf("  mean: %.1f sec",
                as.numeric(Sys.time() - t0, units = "secs")))

t0 <- Sys.time()
agg_max <- terra::resample(raw, tmpl, method = "max")
message(sprintf("  max:  %.1f sec",
                as.numeric(Sys.time() - t0, units = "secs")))

# p90 requires a custom function — terra::resample has no q90 method.
# Aggregate at integer factor (~26), then resample to template grid.
t0 <- Sys.time()
agg_p90_raw <- terra::aggregate(
  raw, fact = ratio,
  fun = function(x, ...) quantile(x, 0.9, na.rm = TRUE)
)
agg_p90 <- terra::resample(agg_p90_raw, tmpl, method = "bilinear")
message(sprintf("  p90:  %.1f sec",
                as.numeric(Sys.time() - t0, units = "secs")))

# --- Stats per method ---
summary_one <- function(r, name) {
  v <- terra::values(r); v <- v[!is.na(v)]
  data.frame(
    method = name,
    n_cells = length(v),
    n_na    = sum(is.na(terra::values(r))),
    min     = min(v),
    median  = median(v),
    mean    = mean(v),
    p90     = quantile(v, 0.90),
    p99     = quantile(v, 0.99),
    max     = max(v)
  )
}

stats <- rbind(
  summary_one(agg_mean, "mean"),
  summary_one(agg_max,  "max"),
  summary_one(agg_p90,  "p90")
)
message("\n--- Per-cell stats across the 1 km template ---")
print(stats, row.names = FALSE, digits = 3)

# --- Pairwise differences ---
message("\n--- Pairwise differences (max-mean, p90-mean, max-p90) ---")
diff_mm  <- agg_max  - agg_mean
diff_pm  <- agg_p90  - agg_mean
diff_mp  <- agg_max  - agg_p90
for (nm in c("max-mean", "p90-mean", "max-p90")) {
  d <- switch(nm, "max-mean" = diff_mm, "p90-mean" = diff_pm, "max-p90" = diff_mp)
  v <- terra::values(d); v <- v[!is.na(v)]
  message(sprintf("  %-9s : median=%5.2f  mean=%5.2f  p95=%5.2f  max=%5.2f",
                  nm, median(v), mean(v), quantile(v, 0.95), max(v)))
}

# --- 3-panel map ---
to_df <- function(r, name) {
  d <- as.data.frame(r, xy = TRUE); names(d) <- c("x", "y", "h")
  d$method <- name; d
}
plot_df <- rbind(
  to_df(agg_mean, "mean"),
  to_df(agg_max,  "max"),
  to_df(agg_p90,  "p90")
)
plot_df$method <- factor(plot_df$method, levels = c("mean", "p90", "max"))

p <- ggplot(plot_df, aes(x = x, y = y, fill = h)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("white", "#cde9c0", "#5ab566", "#1a7337", "#003d12"),
    name = "Tree\nheight\n(m)",
    limits = c(0, 50), oob = scales::squish
  ) +
  facet_wrap(~ method, ncol = 3) +
  coord_equal(expand = FALSE) +
  labs(title = "CHMv2 aggregation comparison — 38 m → 1 km",
       subtitle = "Same colour scale; note how max preserves single-pixel peaks") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggsave("test_chmv2_aggregation.png", p, width = 13, height = 5, dpi = 150)
message("\nSaved test_chmv2_aggregation.png")

# --- Diff map: where do max & mean disagree most? ---
diff_df <- as.data.frame(diff_mm, xy = TRUE)
names(diff_df) <- c("x", "y", "diff")

p2 <- ggplot(diff_df, aes(x = x, y = y, fill = diff)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "max − mean\n(m)"
  ) +
  coord_equal(expand = FALSE) +
  labs(title = "Where do max and mean diverge?",
       subtitle = "High values = isolated tall canopy patches; low = uniform forest or open") +
  theme_minimal(base_size = 11)

ggsave("test_chmv2_aggregation_diff.png", p2, width = 7, height = 6, dpi = 150)
message("Saved test_chmv2_aggregation_diff.png")
