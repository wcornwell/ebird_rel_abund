# test_count_threshold.R
# ─────────────────────────────────────────────────────────────────────────────
# Two questions:
#  (1) Does excluding X-count checklists from the SHARED sampling pool (the
#      pre-zerofill fix) change which checklists each species sees?
#  (2) Where should the high-count cap sit? Sweep it from 50 → Inf and look
#      at CV skill for 3 species with very different count distributions:
#        - Budgerigar              (very long tail, max ≈ 70,000)
#        - White-browed Woodswallow (long tail, max ≈ 3,500)
#        - Maned Duck              (tighter, max ≈ 1,000, q99.9 ≈ 170)
#
# Outputs:
#   analysis/test_count_threshold_results.csv  — per-species × cap CV metrics
#   analysis/count_dist_<species>.png           — log-log count histogram
#   analysis/cv_sweep_<species>.png             — CV metrics vs cap
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf); library(terra); library(geodata); library(ggplot2); library(dplyr)
})

OUT_DIR <- "analysis"
dir.create(OUT_DIR, showWarnings = FALSE)

EBD <- c(
  "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt",
  "ebirdabund/raw_data/ebd_AU-ACT_unv_smp_relMar-2026/ebd_AU-ACT_unv_smp_relMar-2026.txt",
  "ebirdabund/raw_data/ebd_AU-VIC_unv_smp_relMar-2026/ebd_AU-VIC_unv_smp_relMar-2026.txt",
  "ebirdabund/raw_data/ebd_AU-QLD_unv_smp_relMar-2026/ebd_AU-QLD_unv_smp_relMar-2026.txt",
  "ebirdabund/raw_data/ebd_AU-SA_unv_smp_relMar-2026/ebd_AU-SA_unv_smp_relMar-2026.txt"
)
SAMP <- sub("\\.txt$", "_sampling.txt", EBD)
COV_CACHE      <- "ebirdabund_cache"
ZEROFILL_CACHE <- "ebirdabund_cache_nsw_buffer"

SPECIES <- c("Budgerigar", "White-browed Woodswallow", "Maned Duck")
CAPS    <- c(50L, 100L, 200L, 500L, 1000L, .Machine$integer.max)
CAP_LABEL <- function(c) if (c == .Machine$integer.max) "Inf" else as.character(c)

# ── Study polygon: NSW + 100 km buffer (matches run_batch_nsw.R) ─────────────
message("Loading NSW + 100 km buffer polygon...")
aus     <- geodata::gadm("AUS", 1, path = COV_CACHE)
nsw     <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(nsw, 3577), 100000), 4326
)

message("Loading covariates...")
cov <- prepare_covariates(polygon, cache_dir = COV_CACHE)

# ─────────────────────────────────────────────────────────────────────────────
# Part 1: X-count exclusion impact (descriptive only)
# ─────────────────────────────────────────────────────────────────────────────
# The shared-pool fix drops any checklist with X for ANY species. Existing
# caches still include those checklists as zero-detection rows for non-focal
# species. Compare per-species n before/after the new filter using the raw
# cache + global x_count_ids.rds if available.
x_ids_f <- file.path(ZEROFILL_CACHE, "x_count_ids.rds")
if (file.exists(x_ids_f)) {
  x_ids <- readRDS(x_ids_f)
  message(sprintf("\nX-count checklists (any species): %d (global set).",
                  length(x_ids)))
  cat("\n── Part 1: per-species checklist counts ─────────────────────────\n")
  cat(sprintf("  %-30s %12s %12s %10s\n",
              "Species", "old_n_pool", "new_n_pool", "delta"))
  for (sp in SPECIES) {
    cache_f <- file.path(ZEROFILL_CACHE,
                         sprintf("zerofilled_%s.rds",
                                 gsub("[^a-z0-9]+", "_", tolower(sp))))
    if (!file.exists(cache_f)) next
    zf <- readRDS(cache_f)
    # OLD behaviour: drop only focal-species X rows (current clean_ebird).
    old_n <- sum(zf$observation_count != "X")
    # NEW behaviour: drop any checklist in the global X-count set.
    new_n <- sum(!zf$checklist_id %in% x_ids)
    cat(sprintf("  %-30s %12d %12d %10d\n",
                sp, old_n, new_n, new_n - old_n))
  }
} else {
  message("\nNo x_count_ids.rds in ", ZEROFILL_CACHE,
          " — run the pre-cache step in run_batch_nsw.R first, ",
          "or call find_x_count_checklists(EBD) here.")
}

# ─────────────────────────────────────────────────────────────────────────────
# Part 2: count-cap sweep
# ─────────────────────────────────────────────────────────────────────────────
results <- list()

for (sp in SPECIES) {
  message("\n══════════════════════════════════════════════════")
  message("Species: ", sp)
  message("══════════════════════════════════════════════════")
  safe <- gsub("[^a-z0-9]+", "_", tolower(sp))

  # ── Raw count distribution diagnostic (from the cached zerofill) ─────────
  cache_f <- file.path(ZEROFILL_CACHE, sprintf("zerofilled_%s.rds", safe))
  if (file.exists(cache_f)) {
    zf <- readRDS(cache_f)
    raw_pos <- suppressWarnings(as.integer(zf$observation_count))
    raw_pos <- raw_pos[!is.na(raw_pos) & raw_pos > 0]
    q <- quantile(raw_pos, c(0.5, 0.9, 0.95, 0.99, 0.995, 0.999, 1))
    message(sprintf(
      "  Positive counts: n=%d  median=%g  q90=%g  q95=%g  q99=%g  q99.9=%g  max=%g",
      length(raw_pos), q[1], q[2], q[3], q[4], q[6], q[7]
    ))
    hist_df <- data.frame(count = raw_pos)
    p_hist <- ggplot(hist_df, aes(x = count)) +
      geom_histogram(bins = 60, fill = "#2c7bb6") +
      scale_x_log10() + scale_y_log10() +
      geom_vline(xintercept = c(50, 100, 200, 500, 1000),
                 colour = "red", linetype = 2, linewidth = 0.4) +
      labs(title = paste("Positive count distribution:", sp),
           subtitle = "Red dashed lines: 50, 100, 200, 500, 1000 cap thresholds",
           x = "Observation count (log10)", y = "Frequency (log10)") +
      theme_minimal(base_size = 11)
    ggsave(file.path(OUT_DIR, sprintf("count_dist_%s.png", safe)),
           p_hist, width = 8, height = 5, dpi = 130)
  }

  for (cap in CAPS) {
    cap_lab <- CAP_LABEL(cap)
    message(sprintf("\n── cap = %s ────────────────────────────────", cap_lab))

    fit <- tryCatch(
      fit_species_model(
        polygon      = polygon,
        ebird_zip    = EBD,
        sampling_txt = SAMP,
        species      = sp,
        cov_stack    = cov,
        cache_dir    = ZEROFILL_CACHE,
        max_count    = cap
      ),
      ebirdabund_excluded = function(e) {
        message("  [SKIP] ", conditionMessage(e)); NULL
      },
      error = function(e) {
        message("  [FAILED] ", conditionMessage(e)); NULL
      }
    )
    if (is.null(fit)) next

    df       <- fit$data
    dev_expl <- suppressWarnings(summary(fit$model)$dev.expl)
    n_pos    <- sum(df$observation_count > 0L)
    max_obs  <- max(df$observation_count, na.rm = TRUE)

    cv <- evaluate_model_cv(df, k = 5L)

    results[[paste(sp, cap_lab, sep = "|")]] <- data.frame(
      species          = sp,
      cap              = cap_lab,
      cap_num          = if (cap == .Machine$integer.max) Inf else cap,
      n                = nrow(df),
      n_pos            = n_pos,
      max_obs          = max_obs,
      dev_expl_full    = dev_expl,
      cv_spearman      = unname(cv$summary["spearman_r"]),
      cv_pearson       = unname(cv$summary["pearson_r"]),
      cv_rmse          = unname(cv$summary["rmse"]),
      cv_holdout_devx  = unname(cv$summary["holdout_dev_expl"]),
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "  n=%d  n_pos=%d  max=%d  dev.expl=%.3f  cv_spearman=%.3f  cv_holdout_dev=%.3f",
      nrow(df), n_pos, max_obs, dev_expl,
      cv$summary["spearman_r"], cv$summary["holdout_dev_expl"]
    ))
  }
}

# ── Summary table ───────────────────────────────────────────────────────────
summary_df <- do.call(rbind, results)
rownames(summary_df) <- NULL
csv_out <- file.path(OUT_DIR, "test_count_threshold_results.csv")
write.csv(summary_df, csv_out, row.names = FALSE)
message("\nSaved: ", csv_out)

cat("\n══════════════════════════════════════════════════\nSUMMARY\n",
    "══════════════════════════════════════════════════\n", sep = "")
for (sp in SPECIES) {
  sub <- summary_df[summary_df$species == sp, ]
  if (nrow(sub) == 0L) next
  cat(sprintf("\n%s\n", sp))
  cat(sprintf("  %5s  %7s  %7s  %7s  %8s  %8s  %8s  %8s\n",
              "cap", "n", "n_pos", "max", "dev.expl",
              "cv_spear", "cv_pear", "cv_hdev"))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %5s  %7d  %7d  %7d  %8.3f  %8.3f  %8.3f  %8.3f\n",
                sub$cap[i], sub$n[i], sub$n_pos[i], sub$max_obs[i],
                sub$dev_expl_full[i], sub$cv_spearman[i],
                sub$cv_pearson[i], sub$cv_holdout_devx[i]))
  }
}

# ── Per-species sweep plots ─────────────────────────────────────────────────
for (sp in SPECIES) {
  sub <- summary_df[summary_df$species == sp, ]
  if (nrow(sub) == 0L) next
  safe <- gsub("[^a-z0-9]+", "_", tolower(sp))

  plot_df <- tidyr::pivot_longer(
    sub[, c("cap_num", "dev_expl_full", "cv_spearman",
            "cv_pearson", "cv_holdout_devx")],
    cols = -cap_num, names_to = "metric", values_to = "value"
  )

  p <- ggplot(plot_df, aes(x = cap_num, y = value, colour = metric)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    scale_x_log10(breaks = c(50, 100, 200, 500, 1000, 1e6),
                  labels = c("50", "100", "200", "500", "1000", "Inf")) +
    labs(title = paste("Count-cap sweep:", sp),
         x = "Upper count cap (log)", y = "Metric", colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUT_DIR, sprintf("cv_sweep_%s.png", safe)),
         p, width = 8, height = 5, dpi = 130)
}

message("\nDone.")
