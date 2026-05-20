# Paired k-fold CV: baseline GAM vs baseline + s(log1p(speed_kph), k = 4)
#
# Speed (km/h) = effort_distance_km / (duration_minutes / 60). For stationary
# protocols (distance = 0) speed is 0, which log1p handles naturally; the
# protocol_type factor already absorbs the stationary-vs-traveling intercept.
#
# eBird Status & Trends 2023 uses speed as a separate effort variable
# alongside distance and duration — the speed smooth can pick up curvature
# (long fast drives vs. long slow walks) that the additive distance/duration
# smooths cannot.
#
# Reuses ebirdabund_cache_nsw_buffer/{sampling_master.rds, zerofilled_*.rds}.
# Writes:
#   analysis/effort_terms/speed_effort_results.csv
#   analysis/effort_terms/speed_effort_summary.txt

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

test_set <- utils::read.csv(species_csv, stringsAsFactors = FALSE)
species  <- test_set$common_name

master <- readRDS(file.path(cache_dir, "sampling_master.rds"))
message(sprintf("Loaded sampling_master: %d checklists", nrow(master)))

results  <- list()
excluded <- character(0)

t0_all <- Sys.time()

for (sp in species) {
  message(sprintf("\n── [%s] %s", sp, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()

  zf_f <- file.path(cache_dir, sprintf("zerofilled_%s.rds", safe_name(sp)))
  if (!file.exists(zf_f)) {
    message("  Skipping: zerofilled cache not found at ", zf_f)
    excluded <- c(excluded, sp); next
  }

  obs <- readRDS(zf_f)[, c("checklist_id", "observation_count"), drop = FALSE]
  zf  <- merge(master, obs, by = "checklist_id", all.x = TRUE)
  zf$observation_count[is.na(zf$observation_count)] <- "0"
  zf$species_observed <- zf$observation_count != "0"

  df_base <- ebirdabund:::clean_ebird(zf, max_count = max_count)
  set.seed(seed_sub)
  df_base <- ebirdabund:::subsample_hex(df_base, spacing_km = hex_km)

  n_pos <- sum(df_base$observation_count > 0L)
  if (n_pos < 50L) {
    message(sprintf("  Skipping: insufficient detections after subsampling (%d).",
                    n_pos))
    excluded <- c(excluded, sp); next
  }

  # Derived effort: speed in km/h, computed post-clean (distance already
  # imputed to 0 for stationary; duration always >= 10 min so denominator > 0).
  df_var <- df_base
  df_var$speed_kph <- df_var$effort_distance_km /
                      (df_var$duration_minutes / 60)

  # Paired CV: identical fold IDs across both formulas
  set.seed(seed_fold)
  fold_ids <- sample(rep_len(seq_len(k_folds), nrow(df_base)))

  hab_cols     <- detect_hab_cols(df_base)
  formula_base <- ebirdabund:::build_gam_formula(df_base, hab_cols)
  formula_var  <- stats::update(formula_base,
                                . ~ . + s(log1p(speed_kph), k = 4))

  cv_base <- tryCatch(
    evaluate_model_cv(df_base, formula = formula_base, fold_ids = fold_ids),
    error = function(e) { message("  baseline CV failed: ",
                                  conditionMessage(e)); NULL })
  cv_var  <- tryCatch(
    evaluate_model_cv(df_var,  formula = formula_var,  fold_ids = fold_ids),
    error = function(e) { message("  variant CV failed: ",
                                  conditionMessage(e)); NULL })

  if (is.null(cv_base) || is.null(cv_var)) {
    excluded <- c(excluded, sp); next
  }

  metrics <- c("spearman_r", "pearson_r", "mae", "rmse", "holdout_dev_expl")
  row <- data.frame(species = sp, n_pos = n_pos, n_rows = nrow(df_base),
                    stringsAsFactors = FALSE)
  for (m in metrics) {
    row[[paste0(m, "_base")]]  <- unname(cv_base$summary[m])
    row[[paste0(m, "_var")]]   <- unname(cv_var$summary[m])
    row[[paste0(m, "_delta")]] <- unname(cv_var$summary[m] -
                                          cv_base$summary[m])
  }
  results[[sp]] <- row

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf(
    "  done in %5.1fs |  ΔSpearman = %+0.4f |  ΔRMSE = %+0.4f |  ΔDevExpl = %+0.4f",
    dt,
    row$spearman_r_delta, row$rmse_delta, row$holdout_dev_expl_delta))
}

per_species <- do.call(rbind, results)
rownames(per_species) <- NULL

# Per-metric paired Wilcoxon (deltas tested in the better direction)
metric_dir <- list(spearman_r = TRUE, pearson_r = TRUE, mae = FALSE,
                   rmse = FALSE, holdout_dev_expl = TRUE)
summary_rows <- lapply(names(metric_dir), function(m) {
  deltas <- per_species[[paste0(m, "_delta")]]
  higher_better <- metric_dir[[m]]
  alt <- if (higher_better) "greater" else "less"
  p_val <- tryCatch(
    stats::wilcox.test(deltas, alternative = alt, exact = FALSE)$p.value,
    error = function(e) NA_real_)
  wins <- if (higher_better) sum(deltas > 0, na.rm = TRUE)
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
csv_f <- file.path(out_dir, "speed_effort_results.csv")
txt_f <- file.path(out_dir, "speed_effort_summary.txt")
utils::write.csv(per_species, csv_f, row.names = FALSE)

con <- file(txt_f, "w")
writeLines(c(
  "Effort-term test: baseline + s(log1p(speed_kph), k = 4)",
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
message(sprintf("\nTotal wall time: %.1f min",
                as.numeric(difftime(Sys.time(), t0_all, units = "mins"))))
