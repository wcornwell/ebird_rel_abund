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
# peak_doy     : integer DOY override; estimated from model if NULL
# peak_time    : decimal-hour override; estimated from model if NULL
#
# DOY and time are found via the Cornell eBird Best Practices approach: sweep
# each variable over its range at mean habitat + standard effort, then take the
# argmax of the lower 95% CI (fit - 1.96 * se). This reflects when the species
# is most detectable according to the model, not when observers happen to go out.
#
# Returns a list: list(predictions = SpatRaster, peak_doy = int, peak_time = dbl)
predict_abundance <- function(model, pred_surface, polygon,
                              grid_res_km = 1, peak_doy = NULL, peak_time = NULL) {
  train <- model$model
  proto <- ref_protocol(model)

  if (is.null(peak_doy) || is.null(peak_time)) {
    # Build reference row from the prediction surface — it has raw column names
    # (e.g. precip_annual, not log(precip_annual)) that predict.gam expects.
    # The model frame (model$model) stores transformed names so cannot be used.
    non_hab <- c("day_of_year", "time_observations_started", "duration_minutes",
                 "effort_distance_km", "number_observers", "protocol_type",
                 "longitude", "latitude")
    hab_cols <- setdiff(names(pred_surface), non_hab)

    ref_base <- as.data.frame(lapply(
      pred_surface[, hab_cols, drop = FALSE],
      function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1L]
    ))
    ref_base$duration_minutes   <- 60
    ref_base$effort_distance_km <- 1
    ref_base$number_observers   <- 1L
    ref_base$protocol_type      <- factor(proto, levels = levels(train$protocol_type))

    lcl <- function(df) {
      p <- mgcv::predict.gam(model, newdata = df, type = "link", se.fit = TRUE)
      p$fit - 1.96 * p$se.fit
    }

    if (is.null(peak_doy)) {
      seq_doy        <- seq(1L, 365L, by = 1L)
      doy_df         <- ref_base[rep(1L, length(seq_doy)), ]
      doy_df$day_of_year               <- seq_doy
      doy_df$time_observations_started <- 12   # hold time fixed while sweeping DOY
      peak_doy <- seq_doy[which.max(lcl(doy_df))]
      message(sprintf("Standard DOY: %d (model-based peak)", peak_doy))
    }

    if (is.null(peak_time)) {
      seq_tod        <- seq(0, 24, length.out = 300)
      tod_df         <- ref_base[rep(1L, length(seq_tod)), ]
      tod_df$day_of_year               <- peak_doy
      tod_df$time_observations_started <- seq_tod
      peak_time <- seq_tod[which.max(lcl(tod_df))]
      message(sprintf("Peak observation time: %.2f (model-based peak)", peak_time))
    }
  }

  # Add standard-effort columns to full prediction surface
  pred_data <- pred_surface |>
    dplyr::mutate(
      day_of_year               = peak_doy,
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

  list(
    predictions = terra::mask(r_out, poly_vect),
    peak_doy    = peak_doy,
    peak_time   = peak_time
  )
}
