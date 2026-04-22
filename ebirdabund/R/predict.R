# Scan 300 evenly-spaced start times (0-24 h) and return the time
# whose lower confidence limit is maximised -- a conservative peak.
find_peak_time <- function(model, template_row) {
  times   <- seq(0, 24, length.out = 300)
  df_scan <- template_row[rep(1L, 300L), ]
  df_scan$time_observations_started <- times

  preds <- mgcv::predict.gam(
    model, newdata = df_scan, type = "link", se.fit = TRUE
  )
  times[which.max(preds$fit - preds$se.fit)]
}

# Convenience: get the reference protocol_type used during model fitting.
ref_protocol <- function(model) {
  lvls <- levels(model$model$protocol_type)
  if ("Traveling Count" %in% lvls) "Traveling Count" else lvls[1]
}

# Predict relative abundance across the prediction surface.
#
# model        : fitted mgcv::gam from fit_gam()
# pred_surface : data.frame with longitude, latitude + covariate columns
# polygon      : sf polygon (for masking output raster)
# grid_res_km  : resolution used to build pred_surface (sets raster cell size)
# peak_time    : optional numeric override (decimal hours); estimated if NULL
#
# Returns a terra::SpatRaster with two layers: abd and abd_se.
predict_abundance <- function(model, pred_surface, polygon,
                              grid_res_km = 1, peak_time = NULL) {
  train      <- model$model
  proto      <- ref_protocol(model)
  median_doy <- round(stats::median(train$day_of_year, na.rm = TRUE))

  # Build a single template row for peak-time estimation
  template <- pred_surface[1L, , drop = FALSE]
  template$day_of_year          <- median_doy
  template$duration_minutes     <- 60
  template$effort_distance_km   <- 1
  template$number_observers     <- 1
  template$protocol_type        <- factor(
    proto, levels = levels(train$protocol_type)
  )

  if (is.null(peak_time)) {
    message("Estimating peak detection time...")
    peak_time <- find_peak_time(model, template)
    message(sprintf(
      "Peak detection time: %.2f h (%.0f:%02.0f)",
      peak_time, floor(peak_time), (peak_time %% 1) * 60
    ))
  }

  # Add standard-effort columns to full prediction surface
  pred_data <- pred_surface |>
    dplyr::mutate(
      day_of_year               = median_doy,
      time_observations_started = peak_time,
      duration_minutes          = 60,
      effort_distance_km        = 1,
      number_observers          = 1,
      protocol_type             = factor(
        proto, levels = levels(train$protocol_type)
      )
    )

  message(sprintf("Predicting at %d grid points...", nrow(pred_data)))
  preds <- mgcv::predict.gam(
    model, newdata = pred_data, type = "link", se.fit = TRUE
  )

  # as.numeric() strips the row-name attributes mgcv attaches to the vector,
  # which would otherwise confuse terra's type inference.
  pred_data$abd    <- as.numeric(exp(preds$fit))
  pred_data$abd_se <- as.numeric(exp(preds$fit) * preds$se.fit)

  r_abd <- terra::rast(
    data.frame(x = pred_data$longitude, y = pred_data$latitude, z = pred_data$abd),
    type = "xyz", crs = "EPSG:4326"
  )
  r_se <- terra::rast(
    data.frame(x = pred_data$longitude, y = pred_data$latitude, z = pred_data$abd_se),
    type = "xyz", crs = "EPSG:4326"
  )

  r_out        <- c(r_abd, r_se)
  names(r_out) <- c("abd", "abd_se")

  poly_vect <- terra::vect(sf::st_transform(polygon, 4326))
  terra::mask(r_out, poly_vect)
}
