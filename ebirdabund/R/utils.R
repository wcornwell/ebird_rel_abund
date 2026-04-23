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

# Simple ggplot abundance map
plot_abundance <- function(r_pred, polygon, species) {
  abd_df <- as.data.frame(r_pred[["abd"]], xy = TRUE, na.rm = TRUE)
  names(abd_df) <- c("x", "y", "abd")

  thresh <- max(abd_df$abd, na.rm = TRUE) * 0.01
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
