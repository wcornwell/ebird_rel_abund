# Convert HH:MM:SS time string to decimal hours (vectorised)
time_to_decimal <- function(x) {
  vapply(x, function(t) {
    if (is.na(t) || t == "") return(NA_real_)
    parts <- as.numeric(strsplit(t, ":", fixed = TRUE)[[1]])
    parts[1] + parts[2] / 60 + (if (length(parts) >= 3) parts[3] / 3600 else 0)
  }, numeric(1), USE.NAMES = FALSE)
}

# Make a filesystem-safe name from a species string
safe_name <- function(x) {
  gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
}

# Simple ggplot abundance map
plot_abundance <- function(r_pred, polygon, species) {
  abd_df <- as.data.frame(r_pred[["abd"]], xy = TRUE, na.rm = TRUE)
  names(abd_df) <- c("x", "y", "abd")

  abd_df <- abd_df[abd_df$abd > 0.01, ]

  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = abd_df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$abd)
    ) +
    ggplot2::geom_sf(data = sf::st_transform(polygon, 4326),
                     fill = NA, colour = "grey40", linewidth = 0.4) +
    ggplot2::scale_fill_viridis_c(
      option    = "plasma",
      direction = -1,
      name      = "Expected\nbird count"
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(
      title    = paste("Relative Abundance:", species),
      subtitle = "Expected birds detected on a standard 1-hr, 1-km travelling count",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
