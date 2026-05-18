#' Cross-validate a GAM on subsampled eBird data
#'
#' Runs k-fold cross-validation on the subsampled data frame returned by
#' [fit_species_model()] and reports held-out prediction skill.  Pass a custom
#' `formula` to evaluate an alternative covariate set; leave it `NULL` to use
#' the same formula as [fit_gam()] (effort + all habitat covariates).
#'
#' @param df Data frame from `fit_species_model()$data` — subsampled, with
#'   covariate columns already attached.
#' @param formula Optional model formula.  If `NULL`, built automatically via
#'   [build_gam_formula()] (identical to [fit_gam()]).
#' @param k Number of CV folds.  Default `5`.
#' @param seed Random seed for reproducible fold assignment.  Default `42`.
#'   Ignored when `fold_ids` is supplied.
#' @param fold_ids Optional integer vector of length `nrow(df)` assigning each
#'   row to a fold (values in `1..k`). When supplied, `seed` and the internal
#'   random assignment are skipped — enabling paired CV across two formulas on
#'   the same data (e.g. baseline vs. covariate variant).
#'
#' @return A list with:
#' \describe{
#'   \item{fold_metrics}{Data frame with one row per fold: `spearman_r`,
#'     `pearson_r`, `mae`, `rmse`, `holdout_dev_expl`.}
#'   \item{summary}{Named numeric vector: mean of each metric across folds.}
#' }
#'
#' @examples
#' \dontrun{
#' fit   <- fit_species_model(polygon, ebd, samp, "Superb Fairywren", cov)
#' cv    <- evaluate_model_cv(fit$data)
#' cv$summary
#'
#' # Compare effort-only vs. full model
#' cmp <- compare_covariate_models(fit$data)
#' print(cmp)
#' }
#'
#' @export
evaluate_model_cv <- function(df, formula = NULL, k = 5L, seed = 42L,
                              fold_ids = NULL) {
  k <- as.integer(k)
  if (k < 2L || k > nrow(df))
    stop("`k` must be between 2 and nrow(df).")

  if (is.null(formula)) {
    hab_cols <- grep(
      "^(lc_|elevation|precip_|temp_|pop_|water_)", names(df), value = TRUE
    )
    hab_cols <- hab_cols[vapply(hab_cols, function(col) {
      length(unique(stats::na.omit(df[[col]]))) >= 4L
    }, logical(1))]
    formula <- build_gam_formula(df, hab_cols)
  }

  if (is.null(fold_ids)) {
    set.seed(seed)
    fold_id <- sample(rep_len(seq_len(k), nrow(df)))
  } else {
    if (length(fold_ids) != nrow(df))
      stop("`fold_ids` must have length nrow(df).")
    if (!all(fold_ids %in% seq_len(k)))
      stop("`fold_ids` values must be integers in 1..k.")
    fold_id <- as.integer(fold_ids)
  }

  # Protocol reference level
  protocols <- levels(df$protocol_type)
  ref_proto <- if ("Traveling Count" %in% protocols) "Traveling Count" else
    protocols[1L]

  time_knots <- list(time_observations_started = c(0, 24))

  fold_metrics <- lapply(seq_len(k), function(i) {
    train <- df[fold_id != i, ]
    test  <- df[fold_id == i, ]

    train$protocol_type <- stats::relevel(
      droplevels(train$protocol_type), ref = ref_proto
    )
    # Align test levels to train to avoid prediction errors on unseen levels
    test$protocol_type <- factor(
      as.character(test$protocol_type),
      levels = levels(train$protocol_type)
    )

    mod <- tryCatch(
      mgcv::bam(
        formula,
        data     = train,
        family   = mgcv::nb(),
        method   = "fREML",
        discrete = TRUE,
        knots    = time_knots,
        gamma    = 1.4,
        select   = TRUE
      ),
      error = function(e) {
        warning("Fold ", i, " failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(mod)) return(NULL)

    pred <- tryCatch(
      stats::predict(mod, newdata = test, type = "response"),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)

    obs <- test$observation_count

    # Holdout deviance explained: 1 - D(model) / D(null)
    # Null = mean of training counts; use NB log-likelihood ratio
    mu_null <- mean(train$observation_count)
    theta   <- mod$family$getTheta(TRUE)   # estimated NB dispersion
    nb_ll   <- function(y, mu) {
      sum(stats::dnbinom(y, mu = mu, size = theta, log = TRUE))
    }
    ll_mod  <- nb_ll(obs, pred)
    ll_null <- nb_ll(obs, rep(mu_null, length(obs)))
    # Saturated model: mu = obs. For obs = 0, P(0 | NB(0, θ)) = 1 so ll = 0;
    # only non-zero observations contribute. Using pmax(obs, 0.5) for zeros
    # gives mu >> true optimum and makes ll_sat < ll_null for sparse species,
    # flipping the denominator sign and blowing up the ratio.
    ll_sat  <- if (any(obs > 0))
      sum(stats::dnbinom(obs[obs > 0], mu = obs[obs > 0], size = theta, log = TRUE))
    else 0
    dev_expl <- (ll_mod - ll_null) / (ll_sat - ll_null)

    data.frame(
      fold          = i,
      n_test        = length(obs),
      spearman_r    = stats::cor(pred, obs, method = "spearman"),
      pearson_r     = stats::cor(pred, obs),
      mae           = mean(abs(pred - obs)),
      rmse          = sqrt(mean((pred - obs)^2)),
      holdout_dev_expl = dev_expl,
      stringsAsFactors = FALSE
    )
  })

  fold_metrics <- do.call(rbind, Filter(Negate(is.null), fold_metrics))

  metric_cols <- c("spearman_r", "pearson_r", "mae", "rmse", "holdout_dev_expl")
  summary_vec <- colMeans(fold_metrics[, metric_cols], na.rm = TRUE)

  list(fold_metrics = fold_metrics, summary = summary_vec)
}


#' Compare effort-only vs. full habitat model via cross-validation
#'
#' Convenience wrapper around [evaluate_model_cv()] that fits two formulas on
#' the same data and returns a side-by-side comparison:
#' \itemize{
#'   \item \strong{effort_only} — day-of-year, time, duration, distance,
#'     observers, protocol (no habitat covariates).
#'   \item \strong{full} — effort terms + all habitat covariate columns
#'     (names starting with `lc_`, `elevation`, `precip_`, or `temp_`).
#' }
#'
#' Additional named formulas can be passed via `...` to evaluate custom
#' covariate sets.
#'
#' @param df Data frame from `fit_species_model()$data`.
#' @param k Number of CV folds.  Default `5`.
#' @param seed Random seed.  Default `42`.
#' @param ... Named formulas for additional model variants to compare.
#'
#' @return A data frame with one row per model variant and columns for each
#'   cross-validated metric.
#'
#' @export
compare_covariate_models <- function(df, k = 5L, seed = 42L, ...) {

  effort_formula <- build_effort_formula(df)
  full_formula   <- {
    hab_cols <- grep(
      "^(lc_|elevation|precip_|temp_|pop_|water_)", names(df), value = TRUE
    )
    hab_cols <- hab_cols[vapply(hab_cols, function(col) {
      length(unique(stats::na.omit(df[[col]]))) >= 4L
    }, logical(1))]
    build_gam_formula(df, hab_cols)
  }

  model_list <- c(
    list(effort_only = effort_formula, full = full_formula),
    list(...)
  )

  rows <- lapply(names(model_list), function(nm) {
    message(sprintf("  Evaluating model: %s", nm))
    cv  <- evaluate_model_cv(df,
                             formula = model_list[[nm]], k = k, seed = seed)
    cbind(data.frame(model = nm, stringsAsFactors = FALSE),
          as.data.frame(t(cv$summary)))
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  message("\n── Covariate model comparison ─────────────────────────────────────")
  message(sprintf("  %-20s  %8s  %8s  %6s  %6s  %8s",
                  "Model", "Spearman", "Pearson", "MAE", "RMSE", "DevExpl"))
  message("  ", strrep("-", 68))
  for (i in seq_len(nrow(result))) {
    message(sprintf("  %-20s  %8.3f  %8.3f  %6.3f  %6.3f  %8.3f",
                    result$model[i],
                    result$spearman_r[i],
                    result$pearson_r[i],
                    result$mae[i],
                    result$rmse[i],
                    result$holdout_dev_expl[i]))
  }

  invisible(result)
}


# Build an effort-only GAM formula (no habitat covariates).
build_effort_formula <- function(df) {
  sk <- function(var, k) safe_k(df[[var]], k)
  terms <- c(
    sprintf("s(day_of_year, k = %d)",
            sk("day_of_year", 5L)),
    sprintf("s(time_observations_started, bs = 'cc', k = %d)",
            sk("time_observations_started", 4L)),
    sprintf("s(log(duration_minutes), k = %d)",
            sk("duration_minutes", 4L)),
    sprintf("s(log1p(effort_distance_km), k = %d)",
            sk("effort_distance_km", 4L)),
    sprintf("s(number_observers, k = %d)",
            sk("number_observers", 4L)),
    "protocol_type"
  )
  stats::as.formula(
    paste("observation_count ~", paste(terms, collapse = " + "))
  )
}
