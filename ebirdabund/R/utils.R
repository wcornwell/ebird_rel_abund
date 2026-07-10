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

# Most frequent value of a protocol_type vector/factor — used as the GAM
# reference level so predictions are standardised on the dominant protocol.
# Picking an arbitrary (e.g. alphabetically-first) level instead can select a
# rare protocol with too few detections to identify a stable baseline
# coefficient, causing the intercept to blow up toward -Inf (quasi-complete
# separation) and the whole predicted surface to collapse toward zero.
modal_protocol <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
}

# Buffer a polygon by `buffer_km` kilometres using a metric projection.
#
# polygon       : sf or sfc (any CRS).
# buffer_km     : single non-negative numeric, buffer distance in km.
# land_polygon  : optional sf polygon to intersect with after buffering
#                 (e.g., AU mainland) so the fit region doesn't extend
#                 indefinitely over open ocean.
# proj_crs      : metric CRS used for the buffer step. Default is
#                 EPSG:3577 (GDA94 / Australian Albers); override for other
#                 regions.
#
# Returns an sf object in WGS84 (EPSG:4326).
make_fit_polygon <- function(polygon,
                             buffer_km,
                             land_polygon = NULL,
                             proj_crs     = "EPSG:3577") {
  if (!inherits(polygon, c("sf", "sfc"))) {
    stop("`polygon` must be an sf or sfc object.")
  }
  if (!is.numeric(buffer_km) || length(buffer_km) != 1L || is.na(buffer_km) ||
      buffer_km < 0) {
    stop("`buffer_km` must be a single non-negative numeric value.")
  }

  buffered <- sf::st_transform(polygon, proj_crs)
  buffered <- sf::st_buffer(buffered, dist = buffer_km * 1000)
  buffered <- sf::st_transform(buffered, 4326)

  if (!is.null(land_polygon)) {
    if (!inherits(land_polygon, c("sf", "sfc"))) {
      stop("`land_polygon` must be an sf or sfc object.")
    }
    land_wgs <- sf::st_transform(land_polygon, 4326)
    # BirdLife / GADM polygons often have near-duplicate vertices that the S2
    # spherical engine rejects. Disable S2 so GEOS handles the repair.
    old_s2 <- sf::sf_use_s2()
    sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(old_s2), add = TRUE)
    buffered <- sf::st_intersection(
      sf::st_make_valid(buffered),
      sf::st_make_valid(land_wgs)
    )
  }

  buffered
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
# aliases    : optional data frame with columns ebird_sci_name + botw_sci_name
#              for species where eBird and BirdLife disagree on taxonomy
#              (genus reshuffles, gender-agreement renames, splits).
#              On a hit, the BOTW query is rewritten to use botw_sci_name.
#              An empty/NA botw_sci_name short-circuits the lookup (used to
#              record "no BOTW equivalent exists"; avoids a guaranteed-zero
#              SQL round-trip).
#
# Filters to extant/probably-extant (presence 1-2) and established populations
# (origin 1-3, 6; excludes vagrants=4 and uncertain=5). All seasonal codes are
# kept (1-5: resident/breeding/non-breeding/passage/uncertain): nomadic and
# passage species (e.g. Red-necked Avocet, jaegers, sand-plovers) carry their
# coastal/eastern occurrence in passage=4 / uncertain=5 polygons, and excluding
# those clipped their range mid-NSW or left them unmasked entirely.
#
# Returns an sf polygon or NULL if not found.
load_range_botw <- function(sci_name, botw_path, aliases = NULL) {
  if (is.null(sci_name) || is.na(sci_name) || !file.exists(botw_path))
    return(NULL)

  query_name <- sci_name
  if (!is.null(aliases) && nrow(aliases) > 0L) {
    hit <- aliases[!is.na(aliases$ebird_sci_name) &
                     aliases$ebird_sci_name == sci_name, , drop = FALSE]
    if (nrow(hit) > 0L) {
      alias_to <- hit$botw_sci_name[1L]
      if (is.na(alias_to) || !nzchar(alias_to)) {
        message(sprintf(
          "  BOTW: alias map records no BOTW equivalent for '%s' — skipping.",
          sci_name))
        return(NULL)
      }
      message(sprintf("  BOTW: aliasing '%s' -> '%s' (per botw_name_aliases).",
                      sci_name, alias_to))
      query_name <- alias_to
    }
  }

  query <- sprintf(
    "SELECT geom FROM all_species
     WHERE sci_name = '%s'
       AND presence IN (1, 2)
       AND origin   IN (1, 2, 3, 6)
       AND seasonal IN (1, 2, 3, 4, 5)",
    gsub("'", "''", query_name)   # escape any apostrophes in name
  )

  tryCatch(
    {
      r <- sf::st_read(botw_path, query = query, quiet = TRUE)
      if (nrow(r) == 0) NULL else r
    },
    error = function(e) {
      warning("BOTW: error loading range for '", query_name, "': ",
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

# Simple ggplot abundance map.
# border: optional sf polygon drawn as a distinct border line (e.g. NSW state
#   boundary when polygon is a buffered study area).
plot_abundance <- function(r_pred, polygon, species, border = NULL) {
  abd_df <- as.data.frame(r_pred[["abd"]], xy = TRUE, na.rm = TRUE)
  names(abd_df) <- c("x", "y", "abd")

  thresh <- stats::quantile(abd_df$abd, probs = 0.05, na.rm = TRUE)
  abd_df <- abd_df[abd_df$abd > thresh, ]

  poly_wgs84 <- sf::st_transform(polygon, 4326)

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = poly_wgs84, fill = "grey93", colour = NA) +
    ggplot2::geom_raster(
      data = abd_df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$abd)
    ) +
    ggplot2::geom_sf(data = poly_wgs84,
                     fill = NA, colour = "grey40", linewidth = 0.4)

  if (!is.null(border)) {
    p <- p + ggplot2::geom_sf(
      data      = sf::st_transform(border, 4326),
      fill      = NA,
      colour    = "black",
      linewidth = 0.7
    )
  }

  p +
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

# Pairwise scatterplot matrix comparing five species-level abundance metrics
# for all species present in the stack within a given polygon.
#
# The five metrics are:
#   1. model_sum  — sum of the modelled abundance raster over the polygon
#   2. freq_all   — detection rate across all cleaned checklists in the polygon
#   3. abd_all    — mean observation count (incl. zeros) across all checklists
#   4. freq_sub   — detection rate after spatiotemporal subsampling
#   5. abd_sub    — mean observation count after spatiotemporal subsampling
#
# stack          : multi-band SpatRaster (one band per species, safe-name layers)
#                  or a path to one
# polygon        : sf object defining the study area
# cache_dir      : directory containing zerofilled_{species}.rds files
# species_names  : optional named character vector mapping safe-name → display name
# top_label_n    : number of species (by model_sum) to label on every panel
# hex_spacing_km : hex-cell diameter for subsampling (matches model fitting default)
#
# Returns a patchwork ggplot. Save at ~22 × 14 inches for 15 labelled species.
plot_abundance_comparison <- function(stack,
                                      polygon,
                                      cache_dir      = "ebirdabund_cache",
                                      species_names  = NULL,
                                      top_label_n    = 15,
                                      hex_spacing_km = 5) {

  if (is.character(stack)) stack <- terra::rast(stack)

  poly_wgs84 <- sf::st_transform(polygon, 4326)
  poly_vect  <- terra::vect(poly_wgs84)

  # ── (1) Modelled sum in polygon per species ─────────────────────────────────
  model_sum <- as.numeric(
    terra::global(terra::mask(stack, poly_vect), "sum", na.rm = TRUE)[[1L]]
  )
  sp_safe <- names(stack)

  # ── Build hex grid once so dgconstruct() isn't called 500+ times ────────────
  invisible(utils::capture.output(
    suppressMessages(dgg_obj <- dggridR::dgconstruct(spacing = hex_spacing_km))
  ))

  # ── Precompute valid checklist IDs with one polygon intersection ────────────
  # Cache files are bbox-filtered; do the precise polygon test once on the first
  # available species, then reuse the resulting checklist_id set for all others.
  message(sprintf("Loading raw checklist stats for %d species...", length(sp_safe)))
  valid_ids <- NULL
  for (sp_init in sp_safe) {
    cache_init <- file.path(cache_dir, sprintf("zerofilled_%s.rds", sp_init))
    if (!file.exists(cache_init)) next
    zf_init <- tryCatch(readRDS(cache_init), error = function(e) NULL)
    if (is.null(zf_init) || nrow(zf_init) == 0L) next
    zf_init <- tryCatch(clean_ebird(zf_init), error = function(e) NULL)
    if (is.null(zf_init) || nrow(zf_init) == 0L) next
    zf_sf_init <- sf::st_as_sf(zf_init, coords = c("longitude", "latitude"),
                                crs = 4326L, remove = FALSE)
    valid_ids <- sf::st_drop_geometry(
      sf::st_filter(zf_sf_init, poly_wgs84)
    )$checklist_id
    break
  }
  if (is.null(valid_ids)) stop("Could not determine valid checklist IDs from cache.")
  message(sprintf("  %d checklists inside polygon.", length(valid_ids)))

  raw_list <- lapply(seq_along(sp_safe), function(i) {
    sp      <- sp_safe[i]
    cache_f <- file.path(cache_dir, sprintf("zerofilled_%s.rds", sp))
    if (!file.exists(cache_f)) return(NULL)

    zf <- tryCatch(readRDS(cache_f), error = function(e) NULL)
    if (is.null(zf) || nrow(zf) == 0L) return(NULL)

    zf <- tryCatch(clean_ebird(zf), error = function(e) NULL)
    if (is.null(zf) || nrow(zf) == 0L) return(NULL)

    zf <- zf[zf$checklist_id %in% valid_ids, ]
    if (nrow(zf) == 0L) return(NULL)

    freq_all <- mean(zf$species_observed, na.rm = TRUE)
    abd_all  <- mean(zf$observation_count, na.rm = TRUE)

    # Subsample using the shared dgg object (no dgconstruct call per species)
    zf_ss <- tryCatch({
      zf |>
        dplyr::mutate(
          hex_cell = dggridR::dgGEO_to_SEQNUM(dgg_obj,
                                               .data$longitude,
                                               .data$latitude)$seqnum
        ) |>
        dplyr::group_by(.data$year, .data$week, .data$hex_cell) |>
        dplyr::slice_sample(n = 1L) |>
        dplyr::ungroup() |>
        dplyr::select(-"hex_cell")
    }, error = function(e) NULL)

    if (!is.null(zf_ss) && nrow(zf_ss) > 0L) {
      freq_sub <- mean(zf_ss$species_observed, na.rm = TRUE)
      abd_sub  <- mean(zf_ss$observation_count, na.rm = TRUE)
    } else {
      freq_sub <- NA_real_
      abd_sub  <- NA_real_
    }

    data.frame(
      species   = sp,
      model_sum = model_sum[i],
      freq_all  = freq_all,
      abd_all   = abd_all,
      freq_sub  = freq_sub,
      abd_sub   = abd_sub,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, Filter(Negate(is.null), raw_list))
  if (is.null(df) || nrow(df) == 0L) stop("No species data could be computed.")
  # Require model_sum to exceed eps so that species absent from the polygon
  # (tiny raster residuals) don't pile up at log10(eps) = -5 in the plot.
  eps <- 1e-5
  df <- df[!is.na(df$model_sum) & df$model_sum > eps, ]

  # ── Display names ───────────────────────────────────────────────────────────
  if (!is.null(species_names)) {
    labs     <- species_names[df$species]
    fallback <- is.na(labs)
    labs[fallback] <- tools::toTitleCase(gsub("_", " ", df$species[fallback]))
    df$display <- labs
  } else {
    df$display <- tools::toTitleCase(gsub("_", " ", df$species))
  }

  # Notable = top N by modelled sum; their labels appear on every panel
  notable_spp <- df$species[order(-df$model_sum)][seq_len(min(top_label_n, nrow(df)))]
  df$notable  <- df$species %in% notable_spp
  df$label    <- ifelse(df$notable, df$display, "")

  # ── Log10 transform all variables (all span orders of magnitude) ────────────
  # model_sum is already > eps (guaranteed above); no epsilon needed.
  # Rates/counts can be legitimately zero so still need eps.
  plot_df <- df
  plot_df$model_sum <- log10(df$model_sum)
  plot_df$freq_all  <- log10(df$freq_all  + eps)
  plot_df$abd_all   <- log10(df$abd_all   + eps)
  plot_df$freq_sub  <- log10(df$freq_sub  + eps)
  plot_df$abd_sub   <- log10(df$abd_sub   + eps)

  var_meta <- list(
    model_sum = list(label = "Modelled total\n(log₁₀ stack sum)"),
    freq_all  = list(label = "Detection rate\n(log₁₀, all checklists)"),
    abd_all   = list(label = "Mean count\n(log₁₀, all checklists)"),
    freq_sub  = list(label = "Detection rate\n(log₁₀, subsampled)"),
    abd_sub   = list(label = "Mean count\n(log₁₀, subsampled)")
  )
  vars <- names(var_meta)

  use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  # All C(5,2) = 10 pairs, iterating in natural combn order
  pairs <- utils::combn(vars, 2L, simplify = FALSE)

  panels <- lapply(pairs, function(p) {
    xv     <- p[1L]
    yv     <- p[2L]
    sub_df <- plot_df[!is.na(plot_df[[xv]]) & !is.na(plot_df[[yv]]), ]

    g <- ggplot2::ggplot(sub_df,
           ggplot2::aes(x = .data[[xv]], y = .data[[yv]])) +
      ggplot2::geom_point(
        ggplot2::aes(colour = .data$notable),
        alpha = 0.55, size = 1.5
      ) +
      ggplot2::scale_colour_manual(
        values = c("FALSE" = "grey70", "TRUE" = "#3B0F70"),
        guide  = "none"
      ) +
      ggplot2::labs(
        x = var_meta[[xv]]$label,
        y = var_meta[[yv]]$label
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        axis.title   = ggplot2::element_text(size = 7.5),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey85",
                                             linewidth = 0.4),
        plot.margin  = ggplot2::margin(4, 6, 4, 6)
      )

    if (use_ggrepel) {
      g <- g + ggrepel::geom_text_repel(
        ggplot2::aes(label = .data$label),
        size               = 2.2,
        max.overlaps       = 20L,
        segment.colour     = "grey60",
        segment.size       = 0.25,
        min.segment.length = 0.1,
        force              = 2,
        seed               = 42L,
        show.legend        = FALSE
      )
    } else {
      g <- g + ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        size = 2.2, check_overlap = TRUE, vjust = -0.5
      )
    }
    g
  })

  patchwork::wrap_plots(panels, ncol = 4L) +
    patchwork::plot_annotation(
      title    = "Abundance metrics — all pairwise species comparisons",
      subtitle = sprintf(
        "%d species  •  top %d by modelled sum highlighted and labelled",
        nrow(df), sum(df$notable)
      )
    )
}
