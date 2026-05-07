# test_log_transforms.R
# Compare linear vs log-transformed predictors in the GAM.
# Changes being tested:
#   effort: log(duration_minutes), log1p(effort_distance_km)
#   habitat: log(precip_annual), log1p(pop_density)
#   k bumped to 6 for elevation, precip_annual, temp_annual
#
# Primary metric: 5-fold CV Spearman r (distribution-agnostic rank correlation).
# Secondary: held-out deviance explained, RMSE.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
})

EBD   <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"
SAMP  <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"
CACHE <- "ebirdabund_cache"

# Span a range of reporting rates to check robustness across common/uncommon spp
SPECIES <- c(
  "Superb Fairywren",       # very common
  "White-browed Woodswallow", # common, migratory — the species that prompted this
  "Jacky Winter",           # moderately common
  "Painted Honeyeater"      # less common, habitat specialist
)

aus <- geodata::gadm("AUS", 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])

message("Loading covariate stack...")
cov <- prepare_covariates(nsw, cache_dir = CACHE)

# ── Old (linear) formula — reproduces the pre-change behaviour ────────────────
build_gam_formula_linear <- function(df, hab_cols) {
  sk <- function(var, k) safe_k(df[[var]], k)
  effort_terms <- c(
    sprintf("s(day_of_year, bs = 'cc', k = %d)",         sk("day_of_year", 10L)),
    sprintf("s(time_observations_started, bs = 'cc', k = %d)",
            sk("time_observations_started", 10L)),
    sprintf("s(duration_minutes, k = %d)",                sk("duration_minutes", 4L)),
    sprintf("s(effort_distance_km, k = %d)",              sk("effort_distance_km", 4L)),
    sprintf("s(number_observers, k = %d)",                sk("number_observers", 4L)),
    "protocol_type"
  )
  hab_terms <- vapply(hab_cols, function(col) {
    sprintf("s(%s, k = %d)", col, safe_k(df[[col]], 4L))
  }, character(1))
  stats::as.formula(
    paste("observation_count ~", paste(c(effort_terms, hab_terms), collapse = " + "))
  )
}

results <- list()

for (sp in SPECIES) {
  message("\n══════════════════════════════════════════")
  message("Species: ", sp)
  message("══════════════════════════════════════════")

  fit <- tryCatch(
    fit_species_model(
      polygon      = nsw,
      ebird_zip    = EBD,
      sampling_txt = SAMP,
      species      = sp,
      cov_stack    = cov,
      cache_dir    = CACHE
    ),
    error = function(e) { message("  fit_species_model failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) next

  df <- fit$data

  hab_cols <- grep("^(lc_|elevation|precip_|temp_|pop_|water_)", names(df), value = TRUE)
  hab_cols <- hab_cols[vapply(hab_cols, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]

  formula_linear <- build_gam_formula_linear(df, hab_cols)
  formula_log    <- build_gam_formula(df, hab_cols)   # new version with transforms

  message(sprintf("  n checklists: %d  |  n detections: %d",
                  nrow(df), sum(df$observation_count > 0)))
  message("  Linear formula:  ", deparse(formula_linear))
  message("  Log formula:     ", deparse(formula_log))

  message("\nRunning 5-fold CV (linear)...")
  cv_linear <- tryCatch(
    evaluate_model_cv(df, formula = formula_linear, k = 5L),
    error = function(e) { message("  CV failed: ", conditionMessage(e)); NULL }
  )

  message("Running 5-fold CV (log-transformed)...")
  cv_log <- tryCatch(
    evaluate_model_cv(df, formula = formula_log, k = 5L),
    error = function(e) { message("  CV failed: ", conditionMessage(e)); NULL }
  )

  if (is.null(cv_linear) || is.null(cv_log)) next

  message(sprintf("  CV Spearman r  — linear: %.4f  |  log: %.4f  [Δ = %+.4f]",
                  cv_linear$summary["spearman_r"],
                  cv_log$summary["spearman_r"],
                  cv_log$summary["spearman_r"] - cv_linear$summary["spearman_r"]))
  message(sprintf("  CV dev.expl    — linear: %.4f  |  log: %.4f  [Δ = %+.4f]",
                  cv_linear$summary["holdout_dev_expl"],
                  cv_log$summary["holdout_dev_expl"],
                  cv_log$summary["holdout_dev_expl"] - cv_linear$summary["holdout_dev_expl"]))
  message(sprintf("  CV RMSE        — linear: %.4f  |  log: %.4f  [Δ = %+.4f]",
                  cv_linear$summary["rmse"],
                  cv_log$summary["rmse"],
                  cv_log$summary["rmse"] - cv_linear$summary["rmse"]))

  results[[sp]] <- list(
    n                    = nrow(df),
    n_det                = sum(df$observation_count > 0),
    spearman_linear      = cv_linear$summary["spearman_r"],
    spearman_log         = cv_log$summary["spearman_r"],
    dev_expl_linear      = cv_linear$summary["holdout_dev_expl"],
    dev_expl_log         = cv_log$summary["holdout_dev_expl"],
    rmse_linear          = cv_linear$summary["rmse"],
    rmse_log             = cv_log$summary["rmse"]
  )
}

message("\n\n══════════════════════════════════════════════════════════════════")
message("SUMMARY: Linear vs log-transformed predictors (5-fold CV)")
message("══════════════════════════════════════════════════════════════════")
message(sprintf("  %-28s  %6s  %8s  %8s  %8s  %8s  %8s",
                "Species", "n", "Sp.r lin", "Sp.r log", "Δ Sp.r",
                "Δ devExpl", "Δ RMSE"))
message("  ", strrep("-", 88))

for (sp in names(results)) {
  r <- results[[sp]]
  message(sprintf(
    "  %-28s  %6d  %8.4f  %8.4f  %+8.4f  %+9.4f  %+8.4f  [%s]",
    sp, r$n,
    r$spearman_linear, r$spearman_log,
    r$spearman_log - r$spearman_linear,
    r$dev_expl_log - r$dev_expl_linear,
    r$rmse_log - r$rmse_linear,
    if (r$spearman_log > r$spearman_linear) "log better" else "linear better"
  ))
}
