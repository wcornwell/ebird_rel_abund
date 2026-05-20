# Paired CV: does lowering the minimum-duration filter from 10 min to 5 min
# add useful training data, or inflate detection noise?
#
# Comparison is not row-level paired (the two filters yield different
# eligibility sets) but uses a COMMON evaluation set: the held-out fold from
# the long-duration pipeline (df_long_sub). Training:
#   baseline → train on (long minus test fold)
#   variant  → train on (long minus test fold)  ∪  (short-only checklists)
# Both predict on the same test rows from df_long_sub. Skill differences are
# therefore attributable to whether the extra 5–9 min checklists help or hurt
# predictions on the standard (≥10 min) distribution.
#
# Requires clean_ebird() to accept `min_duration_minutes` (parameterised
# 2026-05-18). Reuses ebirdabund_cache_nsw_buffer/{sampling_master, zerofilled_*}.
#
# Writes:
#   analysis/effort_terms/min_duration_results.csv
#   analysis/effort_terms/min_duration_summary.txt

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
})

cache_dir   <- "ebirdabund_cache_nsw_buffer"
species_csv <- "analysis/covariate_test_species.csv"
out_dir     <- "analysis/effort_terms"
k_folds     <- 5L
seed_sub    <- 42L
seed_fold   <- 43L
hex_km      <- 5
max_count   <- 200L
MIN_BASE    <- 10L
MIN_VAR     <- 5L

safe_name <- function(x) {
  s <- gsub("[^A-Za-z0-9]+", "_", tolower(x))
  gsub("^_|_$", "", s)
}

detect_hab_cols <- function(df) {
  hab <- grep("^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height)",
              names(df), value = TRUE)
  hab <- setdiff(hab, "lc_shrubs")
  hab[vapply(hab, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]
}

# Mirror of evaluate_model_cv's single-fold logic, but with a custom train set
# (different from the standard "all rows except the test fold").
fit_and_score <- function(train, test, formula) {
  protocols <- levels(train$protocol_type)
  ref_proto <- if ("Traveling Count" %in% protocols) "Traveling Count" else
    protocols[1L]
  train$protocol_type <- stats::relevel(
    droplevels(train$protocol_type), ref = ref_proto)
  test$protocol_type <- factor(as.character(test$protocol_type),
                               levels = levels(train$protocol_type))
  time_knots <- list(time_observations_started = c(0, 24),
                     day_of_year               = c(0, 365))
  mod <- tryCatch(
    mgcv::bam(formula, data = train, family = mgcv::nb(),
              method = "fREML", discrete = TRUE,
              knots = time_knots, gamma = 1.4, select = TRUE),
    error = function(e) { message("    fit failed: ",
                                  conditionMessage(e)); NULL })
  if (is.null(mod)) return(NULL)
  pred <- tryCatch(
    stats::predict(mod, newdata = test, type = "response"),
    error = function(e) NULL)
  if (is.null(pred)) return(NULL)
  obs <- test$observation_count
  mu_null <- mean(train$observation_count)
  theta   <- mod$family$getTheta(TRUE)
  nb_ll   <- function(y, mu) sum(stats::dnbinom(y, mu = mu, size = theta,
                                                log = TRUE))
  ll_mod  <- nb_ll(obs, pred)
  ll_null <- nb_ll(obs, rep(mu_null, length(obs)))
  ll_sat  <- if (any(obs > 0))
    sum(stats::dnbinom(obs[obs > 0], mu = obs[obs > 0], size = theta,
                       log = TRUE)) else 0
  dev_expl <- (ll_mod - ll_null) / (ll_sat - ll_null)
  list(spearman_r = stats::cor(pred, obs, method = "spearman"),
       pearson_r  = stats::cor(pred, obs),
       mae        = mean(abs(pred - obs)),
       rmse       = sqrt(mean((pred - obs)^2)),
       holdout_dev_expl = dev_expl,
       n_train    = nrow(train),
       n_test     = nrow(test))
}

test_set <- utils::read.csv(species_csv, stringsAsFactors = FALSE)
species  <- test_set$common_name

master <- readRDS(file.path(cache_dir, "sampling_master.rds"))
message(sprintf("Loaded sampling_master: %d checklists", nrow(master)))

results  <- list()
excluded <- character(0)
t0_all   <- Sys.time()

for (sp in species) {
  message(sprintf("\n── [%s] %s", sp, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()

  zf_f <- file.path(cache_dir, sprintf("zerofilled_%s.rds", safe_name(sp)))
  if (!file.exists(zf_f)) {
    message("  Skipping: zerofilled cache not found"); excluded <- c(excluded, sp); next
  }
  obs <- readRDS(zf_f)[, c("checklist_id", "observation_count"), drop = FALSE]
  zf  <- merge(master, obs, by = "checklist_id", all.x = TRUE)
  zf$observation_count[is.na(zf$observation_count)] <- "0"
  zf$species_observed <- zf$observation_count != "0"

  df_long <- ebirdabund:::clean_ebird(zf, max_count = max_count,
                                      min_duration_minutes = MIN_BASE)
  df_full <- ebirdabund:::clean_ebird(zf, max_count = max_count,
                                      min_duration_minutes = MIN_VAR)

  set.seed(seed_sub); df_long_sub <- ebirdabund:::subsample_hex(df_long,
                                                                spacing_km = hex_km)
  set.seed(seed_sub); df_full_sub <- ebirdabund:::subsample_hex(df_full,
                                                                spacing_km = hex_km)

  n_pos_long <- sum(df_long_sub$observation_count > 0L)
  n_pos_full <- sum(df_full_sub$observation_count > 0L)
  if (n_pos_long < 50L) {
    message(sprintf("  Skipping: long-pipeline n_pos = %d < 50", n_pos_long))
    excluded <- c(excluded, sp); next
  }
  message(sprintf(
    "  rows: long=%d (n_pos=%d), full=%d (n_pos=%d), extra short rows=%d",
    nrow(df_long_sub), n_pos_long,
    nrow(df_full_sub), n_pos_full,
    nrow(df_full_sub) - nrow(df_long_sub)))

  # Fold assignment on the COMMON test universe (df_long_sub)
  set.seed(seed_fold)
  fold_ids <- sample(rep_len(seq_len(k_folds), nrow(df_long_sub)))

  hab_cols  <- detect_hab_cols(df_long_sub)
  formula_b <- ebirdabund:::build_gam_formula(df_long_sub, hab_cols)
  formula_v <- ebirdabund:::build_gam_formula(df_full_sub, hab_cols)

  fold_metrics <- list()
  for (i in seq_len(k_folds)) {
    test_ids   <- df_long_sub$checklist_id[fold_ids == i]
    test_rows  <- df_long_sub[fold_ids == i, , drop = FALSE]
    train_base <- df_long_sub[fold_ids != i, , drop = FALSE]
    train_var  <- df_full_sub[!df_full_sub$checklist_id %in% test_ids, ,
                              drop = FALSE]

    m_b <- fit_and_score(train_base, test_rows, formula_b)
    m_v <- fit_and_score(train_var,  test_rows, formula_v)
    if (is.null(m_b) || is.null(m_v)) {
      message(sprintf("    fold %d: skipped (fit/predict failed)", i)); next
    }
    fold_metrics[[i]] <- data.frame(
      fold = i,
      spearman_base = m_b$spearman_r, spearman_var = m_v$spearman_r,
      pearson_base  = m_b$pearson_r,  pearson_var  = m_v$pearson_r,
      mae_base      = m_b$mae,        mae_var      = m_v$mae,
      rmse_base     = m_b$rmse,       rmse_var     = m_v$rmse,
      dev_base      = m_b$holdout_dev_expl, dev_var = m_v$holdout_dev_expl,
      n_train_base  = m_b$n_train,    n_train_var  = m_v$n_train,
      n_test        = m_b$n_test
    )
  }
  fold_df <- do.call(rbind, fold_metrics)
  if (is.null(fold_df) || nrow(fold_df) == 0L) {
    message("  All folds failed; skipping species."); excluded <- c(excluded, sp); next
  }

  agg <- colMeans(fold_df[, -1, drop = FALSE], na.rm = TRUE)
  row <- data.frame(
    species = sp, n_pos_long = n_pos_long, n_pos_full = n_pos_full,
    n_extra_short = nrow(df_full_sub) - nrow(df_long_sub),
    spearman_base = agg["spearman_base"], spearman_var = agg["spearman_var"],
    spearman_delta = agg["spearman_var"] - agg["spearman_base"],
    pearson_base = agg["pearson_base"], pearson_var = agg["pearson_var"],
    pearson_delta = agg["pearson_var"] - agg["pearson_base"],
    mae_base = agg["mae_base"], mae_var = agg["mae_var"],
    mae_delta = agg["mae_var"] - agg["mae_base"],
    rmse_base = agg["rmse_base"], rmse_var = agg["rmse_var"],
    rmse_delta = agg["rmse_var"] - agg["rmse_base"],
    dev_base = agg["dev_base"], dev_var = agg["dev_var"],
    dev_delta = agg["dev_var"] - agg["dev_base"],
    n_train_base = agg["n_train_base"], n_train_var = agg["n_train_var"],
    row.names = NULL, stringsAsFactors = FALSE
  )
  results[[sp]] <- row

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf(
    "  done in %5.1fs |  Δspearman=%+0.4f  ΔRMSE=%+0.4f  Δdev=%+0.4f  (+%d train rows)",
    dt, row$spearman_delta, row$rmse_delta, row$dev_delta,
    round(row$n_train_var - row$n_train_base)))
}

per_species <- do.call(rbind, results)
rownames(per_species) <- NULL

metric_dir <- list(spearman_delta = TRUE, pearson_delta = TRUE,
                   mae_delta = FALSE, rmse_delta = FALSE, dev_delta = TRUE)
summary_rows <- lapply(names(metric_dir), function(m) {
  deltas <- per_species[[m]]
  higher_better <- metric_dir[[m]]
  alt <- if (higher_better) "greater" else "less"
  p_val <- tryCatch(
    stats::wilcox.test(deltas, alternative = alt, exact = FALSE)$p.value,
    error = function(e) NA_real_)
  wins   <- if (higher_better) sum(deltas > 0, na.rm = TRUE)
            else sum(deltas < 0, na.rm = TRUE)
  losses <- if (higher_better) sum(deltas < 0, na.rm = TRUE)
            else sum(deltas > 0, na.rm = TRUE)
  data.frame(metric = m, mean_delta = mean(deltas, na.rm = TRUE),
             median_delta = stats::median(deltas, na.rm = TRUE),
             wilcox_p = p_val, wins = wins,
             ties = sum(deltas == 0, na.rm = TRUE), losses = losses,
             n = sum(!is.na(deltas)), stringsAsFactors = FALSE)
})
summary_df <- do.call(rbind, summary_rows)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
csv_f <- file.path(out_dir, "min_duration_results.csv")
txt_f <- file.path(out_dir, "min_duration_summary.txt")
utils::write.csv(per_species, csv_f, row.names = FALSE)

con <- file(txt_f, "w")
writeLines(c(
  sprintf("Min-duration test: baseline = %d min, variant = %d min",
          MIN_BASE, MIN_VAR),
  "Evaluation: common held-out test set from the long-duration pipeline;",
  "variant trains on long + short, baseline trains on long-only.",
  sprintf("Species evaluated: %d   (excluded: %d)", nrow(per_species),
          length(excluded)),
  if (length(excluded))
    sprintf("Excluded: %s", paste(excluded, collapse = ", ")) else "",
  "",
  sprintf("Total wall time: %.1f min",
          as.numeric(difftime(Sys.time(), t0_all, units = "mins"))),
  "",
  "Per-metric summary (variant - baseline, in the better direction):",
  sprintf("  %-18s  %10s  %10s  %10s  %5s  %5s  %5s",
          "metric", "mean_delta", "median_dlt", "wilcox_p",
          "win", "tie", "loss")
), con)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  writeLines(sprintf("  %-18s  %10.4f  %10.4f  %10.4f  %5d  %5d  %5d",
                     r$metric, r$mean_delta, r$median_delta, r$wilcox_p,
                     r$wins, r$ties, r$losses), con)
}
close(con)

message("\nWrote: ", csv_f)
message("Wrote: ", txt_f)
message(sprintf("Total wall time: %.1f min",
                as.numeric(difftime(Sys.time(), t0_all, units = "mins"))))
