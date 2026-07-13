# test_tweedie_vs_nb.R
# Compare Tweedie vs negative-binomial GAM for two species.
# Primary metric: 5-fold CV Spearman r (distribution-agnostic).
# Secondary: AIC on full data (note: NB is discrete, Tweedie continuous —
# direct AIC comparison is approximate but widely used in ecology).

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
})

# tw() looks up ldTweedie by name in the calling frame; expose it globally.
ldTweedie <- mgcv:::ldTweedie

EBD   <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"
SAMP  <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"
CACHE <- "ebirdabund_cache"

SPECIES <- c("White-browed Woodswallow", "Plumed Whistling-Duck")

aus <- geodata::gadm("AUS", 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])

message("Loading covariate stack...")
cov <- prepare_covariates(nsw, cache_dir = CACHE)

# 5-fold CV for a Tweedie GAM, returning Spearman r and RMSE per fold.
cv_tweedie <- function(df, formula, k = 5L, seed = 42L) {
  time_knots <- list(time_observations_started = c(0, 24))
  ref_proto  <- modal_protocol(df$protocol_type)

  set.seed(seed)
  fold_id <- sample(rep_len(seq_len(k), nrow(df)))

  fold_metrics <- lapply(seq_len(k), function(i) {
    train <- df[fold_id != i, ]
    test  <- df[fold_id == i, ]
    train$protocol_type <- stats::relevel(droplevels(train$protocol_type), ref = ref_proto)
    test$protocol_type  <- factor(as.character(test$protocol_type),
                                  levels = levels(train$protocol_type))

    mod <- tryCatch(
      mgcv::bam(formula, data = train, family = mgcv::tw(),
                method = "fREML", discrete = TRUE,
                knots = time_knots, gamma = 1.4, select = TRUE),
      error = function(e) { warning("Fold ", i, " failed: ", conditionMessage(e)); NULL }
    )
    if (is.null(mod)) return(NULL)

    pred <- tryCatch(predict(mod, newdata = test, type = "response"),
                     error = function(e) NULL)
    if (is.null(pred)) return(NULL)

    obs <- as.numeric(test$observation_count)
    data.frame(
      fold       = i,
      spearman_r = cor(pred, obs, method = "spearman"),
      rmse       = sqrt(mean((pred - obs)^2)),
      stringsAsFactors = FALSE
    )
  })

  fold_metrics <- do.call(rbind, Filter(Negate(is.null), fold_metrics))
  colMeans(fold_metrics[, c("spearman_r", "rmse")], na.rm = TRUE)
}

results <- list()

for (sp in SPECIES) {
  message("\n══════════════════════════════════════════")
  message("Species: ", sp)
  message("══════════════════════════════════════════")

  # Use fit_species_model to get the prepared subsampled data, then re-fit
  # both families explicitly so we can benchmark them head-to-head.
  fit <- fit_species_model(
    polygon      = nsw,
    ebird_zip    = EBD,
    sampling_txt = SAMP,
    species      = sp,
    cov_stack    = cov,
    cache_dir    = CACHE
  )

  df <- fit$data
  df$observation_count <- as.numeric(df$observation_count)

  hab_cols <- grep("^(lc_|elevation|precip_|temp_|pop_|water_)", names(df), value = TRUE)
  hab_cols <- hab_cols[vapply(hab_cols, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]
  formula <- build_gam_formula(df, hab_cols)

  time_knots <- list(time_observations_started = c(0, 24))
  ref_proto  <- modal_protocol(df$protocol_type)
  df$protocol_type <- stats::relevel(df$protocol_type, ref = ref_proto)

  message("\nFitting NB GAM on full data...")
  t_nb <- system.time(
    nb_model <- mgcv::bam(formula, data = df, family = mgcv::nb(),
                          method = "fREML", discrete = TRUE,
                          knots = time_knots, gamma = 1.4, select = TRUE)
  )
  message(sprintf("  NB fit time: %.1f s", t_nb["elapsed"]))

  message("\nFitting Tweedie GAM on full data...")
  t_tw <- system.time(
    tw_model <- tryCatch(
      mgcv::bam(formula, data = df, family = mgcv::tw(),
                method = "fREML", discrete = TRUE,
                knots = time_knots, gamma = 1.4, select = TRUE),
      error = function(e) { message("  Tweedie failed: ", conditionMessage(e)); NULL }
    )
  )

  if (is.null(tw_model)) {
    message("Skipping ", sp, " — Tweedie did not converge.")
    next
  }
  message(sprintf("  Tweedie fit time: %.1f s", t_tw["elapsed"]))

  tw_p   <- tw_model$family$getTheta(TRUE)
  aic_nb <- AIC(nb_model)
  aic_tw <- AIC(tw_model)
  dev_nb <- summary(nb_model)$dev.expl
  dev_tw <- summary(tw_model)$dev.expl

  message(sprintf("\n  %-22s  AIC: %10.2f  dev.expl: %.1f%%  time: %.1fs",
                  "Negative Binomial", aic_nb, dev_nb * 100, t_nb["elapsed"]))
  message(sprintf("  %-22s  AIC: %10.2f  dev.expl: %.1f%%  time: %.1fs  (p=%.3f)",
                  "Tweedie", aic_tw, dev_tw * 100, t_tw["elapsed"], tw_p))
  message(sprintf("  ΔAIC (NB - Tweedie): %.2f  [positive = Tweedie better]", aic_nb - aic_tw))

  message("\nRunning 5-fold CV (NB)...")
  cv_nb <- evaluate_model_cv(df, k = 5L)

  message("Running 5-fold CV (Tweedie)...")
  cv_tw <- cv_tweedie(df, formula, k = 5L)

  message(sprintf("  CV Spearman r  — NB: %.4f  |  Tweedie: %.4f",
                  cv_nb$summary["spearman_r"], cv_tw["spearman_r"]))
  message(sprintf("  CV RMSE        — NB: %.4f  |  Tweedie: %.4f",
                  cv_nb$summary["rmse"], cv_tw["rmse"]))

  results[[sp]] <- list(
    n              = nrow(df),
    aic_nb         = aic_nb,
    aic_tw         = aic_tw,
    delta_aic      = aic_nb - aic_tw,
    dev_nb         = dev_nb,
    dev_tw         = dev_tw,
    tw_p           = tw_p,
    t_nb_s         = t_nb["elapsed"],
    t_tw_s         = t_tw["elapsed"],
    cv_spearman_nb = cv_nb$summary["spearman_r"],
    cv_spearman_tw = cv_tw["spearman_r"],
    cv_rmse_nb     = cv_nb$summary["rmse"],
    cv_rmse_tw     = cv_tw["rmse"]
  )
}

message("\n\n══════════════════════════════════════════════════════════════")
message("SUMMARY: Tweedie vs Negative Binomial")
message("══════════════════════════════════════════════════════════════")
message("Note: AIC is not strictly comparable across NB (discrete) and")
message("Tweedie (continuous) families; CV Spearman r is the fairer metric.\n")

for (sp in names(results)) {
  r <- results[[sp]]
  message(sprintf("  %s  (n=%d)", sp, r$n))
  message(sprintf("    AIC:         NB=%.1f  Tweedie=%.1f  ΔAIC=%.1f  [%s favoured]",
                  r$aic_nb, r$aic_tw, r$delta_aic,
                  if (r$delta_aic > 0) "Tweedie" else "NB"))
  message(sprintf("    dev.expl:    NB=%.1f%%  Tweedie=%.1f%%",
                  r$dev_nb * 100, r$dev_tw * 100))
  message(sprintf("    CV Spearman: NB=%.4f  Tweedie=%.4f  [%s favoured]",
                  r$cv_spearman_nb, r$cv_spearman_tw,
                  if (r$cv_spearman_tw > r$cv_spearman_nb) "Tweedie" else "NB"))
  message(sprintf("    CV RMSE:     NB=%.4f  Tweedie=%.4f",
                  r$cv_rmse_nb, r$cv_rmse_tw))
  message(sprintf("    Tweedie p parameter = %.3f", r$tw_p))
  message(sprintf("    Fit time:    NB=%.1fs  Tweedie=%.1fs\n", r$t_nb_s, r$t_tw_s))
}
