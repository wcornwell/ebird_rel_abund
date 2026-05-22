# test_gamma_sweep_smoke.R
#
# Smoke test for test_gamma_sweep.R: 2 species (one high-N, one mid-N)
# x 2 gammas (1.4 baseline + 6.4 BIC-like extreme) x 5 folds = 20 fits.
# Should complete in a few minutes and verify the pipeline end-to-end
# before launching the full sweep.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
})

GAMMAS         <- c(1.4, 6.4)
SPECIES        <- c("Superb Fairywren", "Rose Robin")
CACHE_DIR      <- "ebirdabund_cache_nsw_buffer"
OUTPUT_DIR     <- "gamma_sweep/smoke"
K_FOLDS        <- 5L
SEED           <- 42L
HEX_SPACING_KM <- 5
MAX_COUNT      <- 200L

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

master_f <- file.path(CACHE_DIR, "sampling_master.rds")
message("Loading sampling master ...")
master <- readRDS(master_f)

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
    stop(sprintf("Insufficient detections (%d).", n_pos))

  set.seed(SEED + 1L)
  fold_ids <- sample(rep_len(seq_len(K_FOLDS), nrow(df)))

  hab_cols <- detect_hab_cols(df)
  formula  <- build_gam_formula(df, hab_cols)

  list(df = df, fold_ids = fold_ids, formula = formula, n_pos = n_pos)
}

rows    <- list()
t_start <- Sys.time()

for (sp in SPECIES) {
  message(sprintf("\n== [%s] ============================================", sp))
  t_sp <- Sys.time()
  prepared <- build_species_data(sp)
  message(sprintf("  n_pos = %d, prep took %.1fs",
                  prepared$n_pos,
                  as.numeric(difftime(Sys.time(), t_sp, units = "secs"))))

  for (g in GAMMAS) {
    message(sprintf("  gamma = %.2f ...", g))
    t_g <- Sys.time()
    cv <- evaluate_model_cv(prepared$df,
                            formula  = prepared$formula,
                            k        = K_FOLDS,
                            fold_ids = prepared$fold_ids,
                            gamma    = g)
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
    message(sprintf("    done in %.1fs  Spearman=%.3f DevExpl=%.3f",
                    as.numeric(difftime(Sys.time(), t_g, units = "secs")),
                    unname(s["spearman_r"]), unname(s["holdout_dev_expl"])))
  }
}

per_species <- do.call(rbind, rows)
out_f <- file.path(OUTPUT_DIR, "smoke_metrics.csv")
utils::write.csv(per_species, out_f, row.names = FALSE)
message("\nTotal elapsed: ",
        format(difftime(Sys.time(), t_start, units = "mins")))
message("Wrote ", out_f)
print(per_species)
