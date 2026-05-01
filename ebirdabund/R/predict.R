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

  # Use the circular mean DOY from detections so that species peaking
  # around Dec-Jan (austral summer) are not mis-assigned to mid-winter.
  # Simple median fails for records spanning the year boundary (DOY 340
  # and DOY 15 have a linear median of ~178, but a circular mean of ~1).
  detected <- train[train$observation_count > 0, ]
  ref_rows <- if (nrow(detected) >= 10L) detected else train
  doy_circ_mean <- function(doy) {
    theta <- doy * 2 * pi / 365
    C     <- mean(cos(theta), na.rm = TRUE)
    S     <- mean(sin(theta), na.rm = TRUE)
    raw   <- atan2(S, C) * 365 / (2 * pi)
    round(((raw %% 365) + 365) %% 365)
  }
  median_doy <- doy_circ_mean(ref_rows$day_of_year)
  if (median_doy == 0L) median_doy <- 365L
  message(sprintf("Standard DOY: %d (circular mean of detections)", median_doy))

  if (is.null(peak_time)) {
    time_circ_mean <- function(t) {
      theta <- t * 2 * pi / 24
      C <- mean(cos(theta), na.rm = TRUE)
      S <- mean(sin(theta), na.rm = TRUE)
      ((atan2(S, C) * 24 / (2 * pi)) %% 24 + 24) %% 24
    }
    peak_time <- time_circ_mean(ref_rows$time_observations_started)
    message(sprintf("Peak observation time: %.2f (circular mean of detections)", peak_time))
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

  # Cap on the log scale before exponentiating to prevent overflow when the GAM
  # extrapolates into covariate regions with no training data. The ceiling is
  # the log of the 90th percentile of observed non-zero counts.
  obs_cap   <- stats::quantile(
    train$observation_count[train$observation_count > 0],
    probs = 0.9, na.rm = TRUE
  )
  log_cap   <- log(max(obs_cap, 1))
  fit_capped       <- pmin(as.numeric(preds$fit), log_cap)
  pred_data$abd    <- exp(fit_capped)
  # Store SE on the link (log) scale so callers can form correct asymmetric
  # CIs: [abd * exp(-z * abd_se), abd * exp(+z * abd_se)]
  pred_data$abd_se <- as.numeric(preds$se.fit)

  # Build raster from an explicit bounding-box template so that concave
  # polygon sections don't leave NA gap rows in the output.
  res_deg   <- grid_res_km / 111.32
  poly_vect <- terra::vect(sf::st_transform(polygon, 4326))
  r_template <- terra::rast(
    terra::ext(poly_vect),
    resolution = res_deg,
    crs        = "EPSG:4326"
  )
  cells <- terra::cellFromXY(
    r_template,
    cbind(pred_data$longitude, pred_data$latitude)
  )

  r_abd <- r_se <- r_template
  terra::values(r_abd) <- NA_real_
  terra::values(r_se)  <- NA_real_
  r_abd[cells] <- pred_data$abd
  r_se[cells]  <- pred_data$abd_se

  r_out        <- c(r_abd, r_se)
  names(r_out) <- c("abd", "abd_se")

  terra::mask(r_out, poly_vect)
}
