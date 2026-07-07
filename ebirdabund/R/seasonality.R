# Per-REZ seasonality: a shared negative-binomial GAM per species that reuses
# the production effort+habitat structure and adds an ordered-factor by-REZ
# cyclic seasonal smooth  s(day_of_year, bs='cc', by = rez).  The reference
# level (statewide) is carried by the global s(day_of_year); each target REZ
# gets a difference smooth, so its own seasonal curve = global + difference.
# Sparse REZs shrink toward the statewide curve (select=TRUE double penalty).
#
# Seasonality metrics per species x region are read off the fitted curve by
# posterior simulation from vcov(fit). Orchestration lives in
# rez_seasonality_batch.R; abundance columns and the WLAB crosswalk are joined
# there (they need terra / the crosswalk, not the fit).

# Fit the shared seasonality GAM. `df` must already carry an ordered factor
# `rez` (reference level first) plus the usual effort/habitat columns.
#' @export
fit_rez_seasonality <- function(df) {
  hab_cols <- grep(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height|nightlights|palsar_hv)",
    names(df), value = TRUE)
  hab_cols <- setdiff(hab_cols, "lc_shrubs")
  hab_cols <- hab_cols[vapply(hab_cols, function(c)
    length(unique(stats::na.omit(df[[c]]))) >= 4L, logical(1))]

  kdoy <- safe_k(df$day_of_year, 10L)
  build_form <- function(simple) {
    bf  <- build_gam_formula(df, hab_cols, simple = simple)  # incl. global s(day_of_year)
    rhs <- paste0(paste(deparse(bf[[3]]), collapse = " "),
                  " + rez + s(day_of_year, bs = 'cc', k = ", kdoy, ", by = rez)")
    stats::as.formula(paste("observation_count ~", rhs))
  }
  full_form   <- build_form(FALSE)
  simple_form <- build_form(TRUE)

  protos <- levels(df$protocol_type)
  ref    <- if ("Traveling Count" %in% protos) "Traveling Count" else protos[1]
  df$protocol_type <- stats::relevel(df$protocol_type, ref = ref)
  knots  <- list(time_observations_started = c(0, 24), day_of_year = c(0, 365))
  gbic   <- log(nrow(df)) / 2

  # Stay on discrete=TRUE for every attempt. The non-discrete bam fallback is
  # ruinously slow on ~334k rows (30+ min per species) and was the cause of the
  # apparent "hangs"; a species that cannot be fit discretely fails fast instead.
  try_fit <- function(form_, sel, g) tryCatch(
    mgcv::bam(form_, data = df, family = mgcv::nb(), method = "fREML",
              discrete = TRUE, select = sel, gamma = g, knots = knots),
    error = function(e) NULL)

  fit <- try_fit(full_form,   TRUE,  gbic)
  if (is.null(fit)) fit <- try_fit(full_form,   TRUE,  1.4)
  if (is.null(fit)) fit <- try_fit(full_form,   FALSE, 1.4)
  if (is.null(fit)) fit <- try_fit(simple_form, FALSE, 1.4)
  if (is.null(fit)) stop("seasonality GAM failed (discrete attempts).")
  list(fit = fit, hab_cols = hab_cols, ref_proto = ref)
}

# Cyclic contiguous window around the peak where the curve is >= threshold.
# Returns c(start_doy, end_doy, wraps) with wraps=TRUE if start > end.
#' @export
seasonality_window <- function(above, peak) {
  n <- length(above)
  if (!above[peak]) return(c(NA, NA, FALSE))
  lo <- peak; while (above[(lo - 2) %% n + 1]) { lo <- (lo - 2) %% n + 1; if (lo == peak) break }
  hi <- peak; while (above[hi %% n + 1])       { hi <- hi %% n + 1;       if (hi == peak) break }
  c(lo, hi, lo > hi)
}

#' @export
doy_to_date <- function(d) if (is.na(d)) NA_character_ else
  format(as.Date(d - 1, origin = "2023-01-01"), "%d %b")

# Posterior-simulation summary of one region's seasonal curve. Returns a
# one-row data.frame of seasonality metrics (no abundance/crosswalk columns).
#' @export
summarise_rez_seasonality <- function(fitobj, df, region, rez_levels,
                                      n_pos, n_chk, n_sim = 1000L,
                                      ratio_thr = 2, prob_thr = 0.95,
                                      min_pos = 20L) {
  fit <- fitobj$fit
  hab <- fitobj$hab_cols
  ref <- data.frame(lapply(df[, hab, drop = FALSE],
                           function(x) mean(x, na.rm = TRUE)))
  ref$duration_minutes   <- 60
  ref$effort_distance_km <- 1
  ref$number_observers   <- 1L
  ref$protocol_type      <- factor(fitobj$ref_proto,
                                   levels = levels(df$protocol_type))
  ref$time_observations_started <- 12
  if ("observer_expertise" %in% names(df))
    ref$observer_expertise <- stats::median(df$observer_expertise, na.rm = TRUE)

  sweep <- ref[rep(1L, 365L), , drop = FALSE]
  sweep$day_of_year <- 1:365
  sweep$rez <- ordered(region, levels = rez_levels)

  Xp    <- stats::predict(fit, newdata = sweep, type = "lpmatrix")
  betas <- mgcv::rmvn(n_sim, stats::coef(fit), stats::vcov(fit))
  eta   <- Xp %*% t(betas)                             # 365 x n_sim (link)

  amp_link      <- apply(eta, 2, function(c) max(c) - min(c))
  prob_seasonal <- mean(exp(amp_link) > ratio_thr)
  eta_med       <- apply(eta, 1, stats::median)
  resp          <- exp(eta_med)
  peak          <- which.max(eta_med)
  suff          <- n_pos >= min_pos
  seasonal      <- suff && prob_seasonal >= prob_thr

  # Bounded seasonality index in [0,1]: 1 - trough/peak on the response scale.
  amp_med <- stats::median(amp_link)
  seas_ix <- 1 - exp(-amp_med)

  win <- if (seasonal) seasonality_window(resp >= mean(resp), peak) else c(NA, NA, FALSE)

  data.frame(
    rez               = region,
    n_checklists_rez  = n_chk,
    n_positive_rez    = n_pos,
    sufficient        = suff,
    seasonality_index = if (suff) seas_ix else NA_real_,
    amplitude_link    = if (suff) amp_med else NA_real_,
    peak_doy          = if (suff) peak else NA_integer_,
    peak_date         = if (suff) doy_to_date(peak) else NA_character_,
    prob_seasonal     = if (suff) prob_seasonal else NA_real_,
    is_seasonal       = if (suff) seasonal else NA,
    window_start_doy  = win[1], window_end_doy = win[2],
    window_start_date = doy_to_date(win[1]), window_end_date = doy_to_date(win[2]),
    window_wraps      = as.logical(win[3]),
    stringsAsFactors  = FALSE
  )
}
