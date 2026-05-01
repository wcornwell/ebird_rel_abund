suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(ggplot2)
  library(dplyr)
})

STACK_PATH <- "species_maps_old/nsw_abundance_stack_3km.tif"
SE_PATH    <- "species_maps_old/nsw_abundance_se_stack_3km.tif"
REZ_PATH   <- "RenewableEnergyZones_Spatial/RenewableEnergyZones.shp"
OUT_DIR    <- "rez_plots"
TOP_N      <- 50

dir.create(OUT_DIR, showWarnings = FALSE)

message("Loading data...")
stack    <- rast(STACK_PATH)
se_stack <- rast(SE_PATH)
rez      <- st_read(REZ_PATH, quiet = TRUE) |>
  st_transform(crs(stack))

# Pretty-print species names from snake_case layer names
pretty_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  # capitalise first letter only
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

message(sprintf("Computing mean abundance for %d species across %d REZs...",
                nlyr(stack), nrow(rez)))

for (i in seq_len(nrow(rez))) {
  rez_name <- rez$REZ_Name[i]
  message(sprintf("  %s", rez_name))

  rez_vect <- vect(rez[i, ])

  cropped <- crop(stack, rez_vect)
  masked  <- mask(cropped, rez_vect)
  means   <- global(masked, fun = "mean", na.rm = TRUE)[[1]]
  names(means) <- names(stack)

  se_cropped <- crop(se_stack, rez_vect)
  se_masked  <- mask(se_cropped, rez_vect)
  mean_se    <- global(se_masked, fun = "mean", na.rm = TRUE)[[1]]
  names(mean_se) <- names(stack)

  top_idx <- order(means, decreasing = TRUE)[
    seq_len(min(TOP_N, sum(!is.na(means))))
  ]
  df <- data.frame(
    species   = pretty_name(names(means)[top_idx]),
    abundance = means[top_idx],
    se_link   = mean_se[top_idx],
    stringsAsFactors = FALSE
  )
  df$species <- factor(df$species, levels = rev(df$species))
  # Asymmetric 95% CI on the log scale, back-transformed
  z <- 1.96
  df$ci_lo <- df$abundance * exp(-z * df$se_link)
  df$ci_hi <- df$abundance * exp( z * df$se_link)

  p <- ggplot(df, aes(x = abundance, y = species)) +
    geom_col(fill = "#2c7bb6") +
    geom_errorbarh(
      aes(xmin = ci_lo, xmax = ci_hi),
      height = 0.4, colour = "grey30", linewidth = 0.4
    ) +
    labs(
      title    = sprintf(
        "Top %d species by mean relative abundance\n%s REZ",
        TOP_N, rez_name
      ),
      x        = "Mean relative abundance (95% CI, log scale)",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank()
    )

  safe_name <- gsub("[^a-z0-9]+", "_", tolower(rez_name))
  out_path  <- file.path(OUT_DIR, sprintf("top%d_%s.png", TOP_N, safe_name))
  ggsave(out_path, p, width = 10, height = 14, dpi = 150)
  message(sprintf("    Saved: %s", out_path))
}

message("\nDone. Plots written to: ", OUT_DIR)
