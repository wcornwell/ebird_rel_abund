utils::globalVariables(".data")

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

# Try to load a range polygon from ebirdst.
# Returns an sf object (possibly multi-row, one per season) or NULL.
load_range_ebirdst <- function(species, resolution) {
  if (!requireNamespace("ebirdst", quietly = TRUE)) return(NULL)

  species_code <- tryCatch(
    ebirdst::get_species(species),
    error = function(e) {
      warning("ebirdst: could not find species code for '", species,
              "': ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(species_code)) return(NULL)

  tryCatch(
    ebirdst::load_ranges(
      species_code, resolution = resolution, smoothed = TRUE
    ),
    error = function(e) {
      warning("ebirdst: could not load range for '", species_code,
              "': ", conditionMessage(e),
              "\n  Fix: ebirdst::ebirdst_download_status('", species_code,
              "', download_ranges = TRUE)")
      NULL
    }
  )
}

# Load a species range from the BirdLife/HBW BOTW geopackage.
#
# botw_path  : path to BOTW_2025.gpkg (or equivalent)
# sci_name   : scientific name matching the BOTW sci_name field
#
# Filters to extant/probably-extant (presence 1-2), established populations
# (origin 1-3, 6; excludes vagrants=4 and uncertain=5), and non-passage
# seasons (seasonal 1-3; excludes passage=4 and uncertain=5).
#
# Returns an sf polygon or NULL if not found.
load_range_botw <- function(sci_name, botw_path) {
  if (is.null(sci_name) || is.na(sci_name) || !file.exists(botw_path))
    return(NULL)

  query <- sprintf(
    "SELECT geom FROM all_species
     WHERE sci_name = '%s'
       AND presence IN (1, 2)
       AND origin   IN (1, 2, 3, 6)
       AND seasonal IN (1, 2, 3)",
    gsub("'", "''", sci_name)   # escape any apostrophes in name
  )

  tryCatch(
    {
      r <- sf::st_read(botw_path, query = query, quiet = TRUE)
      if (nrow(r) == 0) NULL else r
    },
    error = function(e) {
      warning("BOTW: error loading range for '", sci_name, "': ",
              conditionMessage(e))
      NULL
    }
  )
}

# GAM partial-effect smooth plots for one fitted model.
# Returns a patchwork ggplot; caller is responsible for saving.
plot_gam_smooths <- function(mod, species) {
  train      <- mod$model
  term_names <- vapply(mod$smooth, function(s) s$term, character(1L))

  plots <- lapply(term_names, function(vname) {
    if (!vname %in% names(train)) return(NULL)
    x_raw <- train[[vname]]
    xvals <- seq(min(x_raw, na.rm = TRUE), max(x_raw, na.rm = TRUE),
                 length.out = 200L)
    nd <- train[rep(1L, 200L), ]
    nd[[vname]] <- xvals
    for (v in setdiff(term_names, vname)) {
      if (v %in% names(nd) && is.numeric(nd[[v]]))
        nd[[v]] <- mean(train[[v]], na.rm = TRUE)
    }
    p       <- mgcv::predict.gam(mod, newdata = nd,
                                 type = "terms", se.fit = TRUE)
    col_idx <- grep(vname, colnames(p$fit), fixed = TRUE)[1L]
    if (is.na(col_idx)) return(NULL)
    df <- data.frame(
      x  = xvals,
      y  = p$fit[, col_idx],
      lo = p$fit[, col_idx] - 2 * p$se.fit[, col_idx],
      hi = p$fit[, col_idx] + 2 * p$se.fit[, col_idx]
    )
    ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                           alpha = 0.2, fill = "steelblue") +
      ggplot2::geom_line(ggplot2::aes(y = .data$y),
                         colour = "steelblue", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "grey50") +
      ggplot2::labs(x = vname, y = NULL) +
      ggplot2::theme_minimal(base_size = 9)
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0L) return(NULL)

  nc  <- 4L
  nr  <- ceiling(length(plots) / nc)
  patchwork::wrap_plots(plots, ncol = nc) +
    patchwork::plot_annotation(
      title    = paste(species, "— GAM partial effects"),
      subtitle = sprintf(
        "Deviance explained: %.1f%%  |  n = %d checklists",
        summary(mod)$dev.expl * 100, nrow(train)
      )
    )
}

# Simple ggplot abundance map
plot_abundance <- function(r_pred, polygon, species) {
  abd_df <- as.data.frame(r_pred[["abd"]], xy = TRUE, na.rm = TRUE)
  names(abd_df) <- c("x", "y", "abd")

  thresh <- stats::quantile(abd_df$abd, probs = 0.05, na.rm = TRUE)
  abd_df <- abd_df[abd_df$abd > thresh, ]

  poly_wgs84 <- sf::st_transform(polygon, 4326)

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = poly_wgs84, fill = "grey93", colour = NA) +
    ggplot2::geom_raster(
      data = abd_df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$abd)
    ) +
    ggplot2::geom_sf(data = poly_wgs84,
                     fill = NA, colour = "grey40", linewidth = 0.4) +
    ggplot2::scale_fill_gradient(
      low  = "grey90",
      high = "#3B0F70",
      name = "Expected\nbird count"
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(
      title    = paste("Relative Abundance:", species),
      subtitle = paste(
        "Expected birds detected on a standard",
        "1-hr, 1-km travelling count"
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

# Plot species relative abundance within a circular area around a point.
#
# stack      : SpatRaster (multi-band, one band per species) or path to one
# lon / lat  : centre of circle in decimal degrees (WGS84)
# radius_km  : radius of extraction circle in km (default 10)
# top_n      : number of top species to show (default 100)
# species_names : optional named character vector mapping safe-names to display
#                 names; if NULL, safe-names are title-cased for display
#
# Returns a ggplot2 object. For 100 species, save at roughly 22 x 10 inches.
plot_circle_abundance <- function(stack, lon = 149.0041, lat = -32.4553,
                                  radius_km = 10, top_n = 100,
                                  species_names = NULL) {
  if (is.character(stack)) stack <- terra::rast(stack)

  pt     <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  circle <- sf::st_transform(
    sf::st_buffer(sf::st_transform(pt, 32755), dist = radius_km * 1e3),
    4326
  )

  ex <- terra::extract(stack, terra::vect(circle), fun = "mean", na.rm = TRUE)

  abd_df <- data.frame(
    species  = names(stack),
    mean_abd = as.numeric(ex[1L, -1L])
  )
  abd_df <- abd_df[!is.na(abd_df$mean_abd) & abd_df$mean_abd > 0, ]
  abd_df <- abd_df[order(-abd_df$mean_abd), ]
  abd_df <- utils::head(abd_df, top_n)

  if (!is.null(species_names)) {
    labels <- species_names[abd_df$species]
    fallback <- is.na(labels)
    labels[fallback] <- tools::toTitleCase(
      gsub("_", " ", abd_df$species[fallback])
    )
    abd_df$label <- labels
  } else {
    abd_df$label <- tools::toTitleCase(gsub("_", " ", abd_df$species))
  }

  abd_df$label <- factor(abd_df$label, levels = rev(abd_df$label))

  ggplot2::ggplot(abd_df, ggplot2::aes(x = .data$mean_abd, y = .data$label)) +
    ggplot2::geom_col(fill = "#3B0F70", alpha = 0.85, width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$mean_abd)),
      hjust = -0.15, size = 2.5, colour = "grey30"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.18))
    ) +
    ggplot2::labs(
      title    = sprintf(
        "Bird abundance within %g km of %.4f°S, %.4f°E",
        radius_km, abs(lat), lon
      ),
      subtitle = paste(
        "Mean expected birds on a standard 1-hr, 1-km travelling count",
        sprintf("• top %d species by abundance", nrow(abd_df))
      ),
      x = "Mean relative abundance",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_text(size = 7.5),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 8, colour = "grey40")
    )
}
