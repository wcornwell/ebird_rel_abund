# test_gamma_sweep.R
#
# Sweep mgcv::bam's `gamma` complexity multiplier across {1.4, 2, 3, 4.5, 6.4}
# and evaluate held-out predictive skill on the covariate test species set.
#
# Why: gamma=1.4 is the standard mild anti-overfit setting (Wood 2017).
# At n=334k checklists, log(n)/2 ≈ 6.4 corresponds to BIC-like behaviour —
# strongly preferring simpler models. We want to know empirically where, if
# anywhere, holdout deviance starts degrading.
#
# Uses paired k-fold CV: identical fold assignments across all gamma values
# per species, so deltas are matched-pair comparisons. Reuses
# sampling_master.rds and zerofilled_*.rds caches — never rewrites them.
#
# Runtime: ~30-60 min serial for 19 species × 5 gammas × 5 folds = 475 fits.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
})

# ── Config ──────────────────────────────────────────────────────────────────

GAMMAS         <- c(1.4, 2.0, 3.0, 4.5, 6.4)
SPECIES_CSV    <- "analysis/covariate_test_species.csv"
CACHE_DIR      <- "ebirdabund_cache_nsw_buffer"
OUTPUT_DIR     <- "gamma_sweep"
K_FOLDS        <- 5L
SEED           <- 42L
HEX_SPACING_KM <- 5
MAX_COUNT      <- 200L

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load shared inputs ──────────────────────────────────────────────────────

master_f <- file.path(CACHE_DIR, "sampling_master.rds")
if (!file.exists(master_f))
  stop("sampling_master.rds not found in ", CACHE_DIR,
       " — run run_batch_nsw.R first to build it.")

message("Loading sampling master ...")
master <- readRDS(master_f)

test_set <- utils::read.csv(SPECIES_CSV, stringsAsFactors = FALSE)
species  <- test_set$common_name

message(sprintf("Species: %d   Gammas: %s   Folds: %d",
                length(species),
                paste(GAMMAS, collapse = ", "),
                K_FOLDS))

# ── Build per-species training data + fold IDs once ─────────────────────────
# Mirrors paired_cv_one_species() in test_covariate.R so the inputs are
# identical to what variant testing uses.

build_species_data <- function(sp) {
  zf_f <- file.path(CACHE_DIR, sprintf("zerofilled_%s.rds", safe_name(sp)))
  if (!file.exists(zf_f))
    stop("zerofilled cache not found: ", zf_f)

  obs <- readRDS(zf_f)[, c("checklist_id", "observation_count"), drop = FALSE]
  zf  <- merge(master, obs, by = "checklist_id", all.x = TRUE)
  zf$observation_count[is.na(zf$observation_count)] <- "0"
  zf$species_observed <- zf$observation_count != "0"

  df <- clean_ebird(zf, max_count = MAX_COUNT)

  set.seed(SEED)
  df <- subsample_hex(df, spacing_km = HEX_SPACING_KM)

  n_pos <- sum(df$observation_count > 0L)
  if (n_pos < 50L)
    stop(sprintf("Insufficient detections after subsampling (%d, need >= 50).",
                 n_pos))

  set.seed(SEED + 1L)
  fold_ids <- sample(rep_len(seq_len(K_FOLDS), nrow(df)))

  hab_cols <- detect_hab_cols(df)
  formula  <- build_gam_formula(df, hab_cols)

  list(df = df, fold_ids = fold_ids, formula = formula, n_pos = n_pos)
}

# ── Sweep ───────────────────────────────────────────────────────────────────

rows     <- list()
excluded <- character(0)
t_start  <- Sys.time()

for (sp in species) {
  message(sprintf("\n── [%s] ─────────────────────────────────────────", sp))

  prepared <- tryCatch(build_species_data(sp),
                       error = function(e) {
                         message("  Skipping: ", conditionMessage(e))
                         excluded <<- c(excluded, sp)
                         NULL
                       })
  if (is.null(prepared)) next

  for (g in GAMMAS) {
    message(sprintf("  gamma = %.2f ...", g))
    cv <- tryCatch(
      evaluate_model_cv(prepared$df,
                        formula  = prepared$formula,
                        k        = K_FOLDS,
                        fold_ids = prepared$fold_ids,
                        gamma    = g),
      error = function(e) {
        message("    Failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(cv)) next

    s <- cv$summary
    rows[[length(rows) + 1L]] <- data.frame(
      species          = sp,
      n_pos            = prepared$n_pos,
      gamma            = g,
      spearman_r       = unname(s["spearman_r"]),
      pearson_r        = unname(s["pearson_r"]),
      mae              = unname(s["mae"]),
      rmse             = unname(s["rmse"]),
      holdout_dev_expl = unname(s["holdout_dev_expl"]),
      stringsAsFactors = FALSE
    )
  }
}

if (length(rows) == 0L) stop("No species produced metrics.")

per_species <- do.call(rbind, rows)

per_species_f <- file.path(OUTPUT_DIR, "per_species_metrics.csv")
utils::write.csv(per_species, per_species_f, row.names = FALSE)
message("\nWrote ", per_species_f)

# ── Summary: mean metrics by gamma, plus paired deltas vs baseline (1.4) ────

metric_cols <- c("spearman_r", "pearson_r", "mae", "rmse", "holdout_dev_expl")

summary_by_gamma <- aggregate(
  per_species[, metric_cols],
  by = list(gamma = per_species$gamma),
  FUN = mean, na.rm = TRUE
)

baseline_gamma <- 1.4
base <- per_species[per_species$gamma == baseline_gamma, , drop = FALSE]

delta_rows <- lapply(setdiff(GAMMAS, baseline_gamma), function(g) {
  sub <- per_species[per_species$gamma == g, , drop = FALSE]
  m   <- merge(base[, c("species", metric_cols)],
               sub[,  c("species", metric_cols)],
               by = "species", suffixes = c("_base", "_g"))
  out <- data.frame(gamma = g, stringsAsFactors = FALSE)
  for (mt in metric_cols) {
    d <- m[[paste0(mt, "_g")]] - m[[paste0(mt, "_base")]]
    higher_better <- mt %in% c("spearman_r", "pearson_r", "holdout_dev_expl")
    out[[paste0(mt, "_mean_delta")]] <- mean(d, na.rm = TRUE)
    out[[paste0(mt, "_wins")]] <- if (higher_better)
      sum(d > 0, na.rm = TRUE) else sum(d < 0, na.rm = TRUE)
    out[[paste0(mt, "_losses")]] <- if (higher_better)
      sum(d < 0, na.rm = TRUE) else sum(d > 0, na.rm = TRUE)
  }
  out
})
deltas <- do.call(rbind, delta_rows)

summary_f <- file.path(OUTPUT_DIR, "summary.txt")
con <- file(summary_f, "w")
writeLines(c(
  sprintf("Gamma sweep — %d species, %d folds, seed %d",
          length(unique(per_species$species)), K_FOLDS, SEED),
  sprintf("Started:  %s", format(t_start)),
  sprintf("Finished: %s", format(Sys.time())),
  if (length(excluded))
    sprintf("Excluded: %d (%s)", length(excluded),
            paste(excluded, collapse = ", "))
  else "Excluded: 0",
  "",
  "── Mean metrics by gamma ───────────────────────────────────────────",
  sprintf("  %6s  %8s  %8s  %7s  %7s  %8s",
          "gamma", "Spearman", "Pearson", "MAE", "RMSE", "DevExpl")
), con)
for (i in seq_len(nrow(summary_by_gamma))) {
  r <- summary_by_gamma[i, ]
  writeLines(sprintf("  %6.2f  %8.4f  %8.4f  %7.3f  %7.3f  %8.4f",
                     r$gamma, r$spearman_r, r$pearson_r,
                     r$mae, r$rmse, r$holdout_dev_expl), con)
}
writeLines(c("",
  sprintf("── Paired deltas vs gamma = %.1f (mean across species) ──────────",
          baseline_gamma),
  sprintf("  %6s  %12s  %12s  %12s  %12s  %12s",
          "gamma", "dSpearman", "dPearson", "dMAE", "dRMSE", "dDevExpl")
), con)
for (i in seq_len(nrow(deltas))) {
  r <- deltas[i, ]
  writeLines(sprintf("  %6.2f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f",
                     r$gamma,
                     r$spearman_r_mean_delta,
                     r$pearson_r_mean_delta,
                     r$mae_mean_delta,
                     r$rmse_mean_delta,
                     r$holdout_dev_expl_mean_delta), con)
}
writeLines(c("",
  "── Win/loss counts vs baseline (species moving in better direction) ─",
  sprintf("  %6s  %12s  %12s  %12s  %12s  %12s",
          "gamma", "Spearman", "Pearson", "MAE", "RMSE", "DevExpl")
), con)
for (i in seq_len(nrow(deltas))) {
  r <- deltas[i, ]
  fmt <- function(w, l) sprintf("%d / %d", w, l)
  writeLines(sprintf("  %6.2f  %12s  %12s  %12s  %12s  %12s",
                     r$gamma,
                     fmt(r$spearman_r_wins, r$spearman_r_losses),
                     fmt(r$pearson_r_wins,  r$pearson_r_losses),
                     fmt(r$mae_wins,        r$mae_losses),
                     fmt(r$rmse_wins,       r$rmse_losses),
                     fmt(r$holdout_dev_expl_wins,
                         r$holdout_dev_expl_losses)), con)
}
writeLines(c("",
  "Convention: 'win' = gamma_g moved in the better direction vs gamma=1.4",
  "(higher Spearman/Pearson/DevExpl, lower MAE/RMSE)."
), con)
close(con)

message("Wrote ", summary_f)
message("\nDone. Inspect ", summary_f, " for the summary table.")
