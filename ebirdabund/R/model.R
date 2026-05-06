# Choose k for a smooth term: at most default_k, but never more than
# (n_unique - 1). Minimum of 3 so the spline is identifiable.
safe_k <- function(x, default_k) {
  n_uniq <- length(unique(stats::na.omit(x)))
  as.integer(min(default_k, max(3L, n_uniq - 1L)))
}
# Build a GAM formula with data-driven k for each smooth term.
build_gam_formula <- function(df, hab_cols) {
  sk <- function(var, k) safe_k(df[[var]], k)
  effort_terms <- c(
    sprintf("s(day_of_year, bs = 'cc', k = %d)",
            sk("day_of_year", 10L)),
    sprintf("s(time_observations_started, bs = 'cc', k = %d)",
            sk("time_observations_started", 10L)),
    sprintf("s(log(duration_minutes), k = %d)",
            sk("duration_minutes", 4L)),
    sprintf("s(log1p(effort_distance_km), k = %d)",
            sk("effort_distance_km", 4L)),
    sprintf("s(number_observers, k = %d)",
            sk("number_observers", 4L)),
    "protocol_type"
  )

  # Right-skewed covariates spanning wide ranges benefit from log scaling so
  # knots are distributed more evenly across the data.
  log1p_cols <- c("pop_density", "tree_height")  # include zeros
  log_cols   <- "precip_annual"                  # strictly positive
  k6_cols    <- c("elevation", "precip_annual", "temp_annual")

  hab_terms <- vapply(hab_cols, function(col) {
    var <- if (col %in% log1p_cols) sprintf("log1p(%s)", col)
           else if (col %in% log_cols) sprintf("log(%s)", col)
           else col
    k   <- if (col %in% k6_cols) 6L else 4L
    sprintf("s(%s, k = %d)", var, safe_k(df[[col]], k))
  }, character(1))

  stats::as.formula(
    paste("observation_count ~",
          paste(c(effort_terms, hab_terms), collapse = " + "))
  )
}

# Fit a negative-binomial GAM to training data.
# Returns the fitted mgcv::gam object.
fit_gam <- function(df) {
  hab_cols <- grep(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height)", names(df), value = TRUE
  )
  hab_cols <- setdiff(hab_cols, "lc_shrubs")

  if (length(hab_cols) == 0) {
    stop("No habitat covariate columns found (expected lc_*, ",
         "elevation, precip_*, temp_*, or pop_*).")
  }

  # Drop any hab col with < 4 unique values (can't build even a k=3 spline)
  hab_cols <- hab_cols[vapply(hab_cols, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]

  formula <- build_gam_formula(df, hab_cols)

  # Set "Traveling Count" as reference level when present, else most common
  protocols <- levels(df$protocol_type)
  ref_proto <- if ("Traveling Count" %in% protocols) "Traveling Count" else
    protocols[1]
  df$protocol_type <- stats::relevel(df$protocol_type, ref = ref_proto)

  time_knots <- list(time_observations_started = c(0, 24),
                     day_of_year               = c(0, 365))

  message("Fitting negative-binomial GAM (", nrow(df), " checklists)...")

  fit_bam <- function(select) {
    converged <- TRUE
    mod <- withCallingHandlers(
      tryCatch(
        mgcv::bam(
          formula,
          data     = df,
          family   = mgcv::nb(),
          method   = "fREML",
          discrete = TRUE,
          knots    = time_knots,
          gamma    = 1.4,
          select   = select
        ),
        error = function(e) {
          converged <<- FALSE
          NULL
        }
      ),
      warning = function(w) {
        if (grepl("did not converge", conditionMessage(w), fixed = TRUE)) {
          converged <<- FALSE
          invokeRestart("muffleWarning")
        }
      }
    )
    list(mod = mod, converged = converged)
  }

  res <- fit_bam(select = TRUE)

  if (!res$converged || is.null(res$mod)) {
    message("  Attempt 1 (select=TRUE) did not converge; ",
            "retrying without term selection.")
    res <- fit_bam(select = FALSE)
  }

  if (!res$converged || is.null(res$mod)) {
    message("  Attempt 2 (select=FALSE) did not converge; ",
            "retrying with discrete=FALSE.")
    res$mod <- mgcv::bam(
      formula,
      data     = df,
      family   = mgcv::nb(),
      method   = "fREML",
      discrete = FALSE,
      knots    = time_knots,
      gamma    = 1.4,
      nthreads = 1L
    )
  }

  mod <- res$mod
  print_predictor_importance(mod)
  mod
}

# Print a ranked table of predictor importance from a fitted GAM.
print_predictor_importance <- function(mod) {
  gam_sum  <- suppressWarnings(summary(mod))
  gam_anov <- suppressWarnings(anova(mod))

  # Smooth terms: EDF and approximate F from summary
  stbl <- as.data.frame(gam_sum$s.table)
  rows <- data.frame(
    term    = rownames(stbl),
    edf     = stbl$edf,
    f_stat  = stbl[, "F"],
    p_value = stbl[, "p-value"],
    stringsAsFactors = FALSE
  )

  # Parametric terms: overall F from anova (e.g. protocol_type)
  if (!is.null(gam_anov$pTerms.table) && nrow(gam_anov$pTerms.table) > 0) {
    ptbl <- as.data.frame(gam_anov$pTerms.table)
    prows <- data.frame(
      term    = rownames(ptbl),
      edf     = ptbl$df,
      f_stat  = ptbl$F,
      p_value = ptbl[, "p-value"],
      stringsAsFactors = FALSE
    )
    rows <- rbind(rows, prows)
  }

  rows <- rows[order(-rows$f_stat), ]

  message("\nPredictor importance (sorted by F-statistic):")
  message(sprintf("  %-42s  %5s  %9s  %9s", "Term", "EDF", "F", "p-value"))
  message("  ", strrep("-", 70))
  for (i in seq_len(nrow(rows))) {
    message(sprintf("  %-42s  %5.2f  %9.2f  %9s",
      rows$term[i],
      rows$edf[i],
      rows$f_stat[i],
      format.pval(rows$p_value[i], digits = 2, eps = 1e-4)
    ))
  }
  message("")
}
