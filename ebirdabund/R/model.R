# Choose k for a smooth term: at most default_k, but never more than
# (n_unique - 1). Minimum of 3 so the spline is identifiable.
safe_k <- function(x, default_k) {
  n_uniq <- length(unique(stats::na.omit(x)))
  as.integer(min(default_k, max(3L, n_uniq - 1L)))
}

# Build a GAM formula with data-driven k for each smooth term.
build_gam_formula <- function(df, hab_cols) {
  effort_terms <- c(
    sprintf("s(day_of_year,                k = %d)", safe_k(df$day_of_year,                5L)),
    sprintf("s(time_observations_started,  bs = 'cc', k = %d)",
                                                      safe_k(df$time_observations_started, 7L)),
    sprintf("s(duration_minutes,           k = %d)", safe_k(df$duration_minutes,           5L)),
    sprintf("s(effort_distance_km,         k = %d)", safe_k(df$effort_distance_km,         5L)),
    sprintf("s(number_observers,           k = %d)", safe_k(df$number_observers,           5L)),
    "protocol_type"
  )

  hab_terms <- vapply(hab_cols, function(col) {
    sprintf("s(%s, k = %d)", col, safe_k(df[[col]], 5L))
  }, character(1))

  stats::as.formula(
    paste("observation_count ~",
          paste(c(effort_terms, hab_terms), collapse = " + "))
  )
}

# Fit a negative-binomial GAM to training data.
# Returns the fitted mgcv::gam object.
fit_gam <- function(df) {
  hab_cols <- grep("^(lc_|elevation)", names(df), value = TRUE)

  if (length(hab_cols) == 0) {
    stop("No habitat covariate columns found (expected names starting with 'lc_' or 'elevation').")
  }

  # Drop any hab col with < 4 unique values (can't build even a k=3 spline)
  hab_cols <- hab_cols[vapply(hab_cols, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]

  formula <- build_gam_formula(df, hab_cols)

  # Set "Traveling Count" as reference level when present, else most common
  protocols <- levels(df$protocol_type)
  ref_proto <- if ("Traveling Count" %in% protocols) "Traveling Count" else protocols[1]
  df$protocol_type <- stats::relevel(df$protocol_type, ref = ref_proto)

  time_knots <- list(time_observations_started = c(0, 24))

  message("Fitting negative-binomial GAM (", nrow(df), " checklists)...")

  mgcv::bam(
    formula,
    data     = df,
    family   = mgcv::nb(),
    method   = "fREML",
    discrete = TRUE,
    knots    = time_knots
  )
}
