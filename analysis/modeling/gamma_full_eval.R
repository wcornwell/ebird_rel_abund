# gamma_full_eval.R
#
# Full evaluation of gamma = {1.4, 6.4} on the covariate test species set.
# Per species: paired k-fold CV (identical fold IDs across gammas) plus a
# full-data fit at each gamma, with overlaid partial-effect smooth plots.
#
# Outputs (under gamma_sweep/full/):
#   per_species_metrics.csv   - appended after each species
#   smooths_{species}.png     - one per species, both gammas overlaid
#   run.log                   - tee'd stdout/stderr
#
# Resume: a species is skipped if both its CSV row and PNG already exist.
# Runtime: ~3-4 hours serial for 19 species (12 fits per species).

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(ggplot2)
  library(patchwork)
})

GAMMAS         <- c(1.4, 6.4)
SPECIES_CSV    <- "analysis/covariate_test_species.csv"
CACHE_DIR      <- "ebirdabund_cache_nsw_buffer"
OUTPUT_DIR     <- "gamma_sweep/full"
METRICS_F      <- file.path(OUTPUT_DIR, "per_species_metrics.csv")
K_FOLDS        <- 5L
SEED           <- 42L
HEX_SPACING_KM <- 5
MAX_COUNT      <- 200L

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Resume helpers ──────────────────────────────────────────────────────────

done_species <- function() {
  if (!file.exists(METRICS_F)) return(character(0))
  d <- utils::read.csv(METRICS_F, stringsAsFactors = FALSE)
  have_rows <- table(d$species)
  expected  <- length(GAMMAS)
  csv_done  <- names(have_rows)[have_rows >= expected]
  png_done  <- vapply(csv_done, function(sp)
    file.exists(file.path(OUTPUT_DIR, sprintf("smooths_%s.png", safe_name(sp)))),
    logical(1L))
  csv_done[png_done]
}

append_rows <- function(rows) {
  df <- do.call(rbind, rows)
  if (file.exists(METRICS_F)) {
    utils::write.table(df, METRICS_F, sep = ",", append = TRUE,
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else {
    utils::write.csv(df, METRICS_F, row.names = FALSE)
  }
}

# ── Shared inputs ───────────────────────────────────────────────────────────

master_f <- file.path(CACHE_DIR, "sampling_master.rds")
if (!file.exists(master_f))
  stop("sampling_master.rds not found in ", CACHE_DIR)

message("Loading sampling master ...")
master <- readRDS(master_f)

test_set <- utils::read.csv(SPECIES_CSV, stringsAsFactors = FALSE)
all_species <- test_set$common_name

already <- done_species()
species  <- setdiff(all_species, already)

message(sprintf("Total species: %d   Already done: %d   To run: %d",
                length(all_species), length(already), length(species)))
if (length(already)) message("  Skipping: ", paste(already, collapse = ", "))

# ── Per-species helpers ─────────────────────────────────────────────────────

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

fit_full <- function(df, formula, g) {
  mgcv::bam(stats::as.formula(formula),
            data     = df,
            family   = mgcv::nb(),
            method   = "fREML",
            discrete = TRUE,
            gamma    = g,
            select   = TRUE)
}

smooth_overlay <- function(mods, vname) {
  ref   <- mods[[1]]
  train <- ref$model
  if (!vname %in% names(train)) return(NULL)
  x_raw <- train[[vname]]
  xvals <- seq(min(x_raw, na.rm = TRUE), max(x_raw, na.rm = TRUE),
               length.out = 200L)
  nd <- train[rep(1L, 200L), ]
  nd[[vname]] <- xvals
  term_names <- vapply(ref$smooth, function(s) s$term, character(1L))
  for (v in setdiff(term_names, vname)) {
    if (v %in% names(nd) && is.numeric(nd[[v]]))
      nd[[v]] <- mean(train[[v]], na.rm = TRUE)
  }

  rows <- lapply(names(mods), function(lab) {
    p <- mgcv::predict.gam(mods[[lab]], newdata = nd,
                           type = "terms", se.fit = TRUE)
    col_idx <- grep(vname, colnames(p$fit), fixed = TRUE)[1L]
    if (is.na(col_idx)) return(NULL)
    data.frame(
      x     = xvals,
      y     = p$fit[, col_idx],
      lo    = p$fit[, col_idx] - 2 * p$se.fit[, col_idx],
      hi    = p$fit[, col_idx] + 2 * p$se.fit[, col_idx],
      gamma = lab,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(df) || nrow(df) == 0L) return(NULL)

  ggplot(df, aes(x = x, colour = gamma, fill = gamma)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(aes(y = y), linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("1.4" = "steelblue", "6.4" = "firebrick")) +
    scale_fill_manual(values   = c("1.4" = "steelblue", "6.4" = "firebrick")) +
    labs(x = vname, y = NULL) +
    theme_minimal(base_size = 9)
}

save_smooth_panel <- function(mods, sp, n_pos) {
  ref_mod    <- mods[[1]]
  term_names <- vapply(ref_mod$smooth, function(s) s$term, character(1L))
  panels <- Filter(Negate(is.null),
                   lapply(term_names, function(v) smooth_overlay(mods, v)))
  if (length(panels) == 0L) return(invisible(NULL))

  # Strip individual legends, keep one shared on the bottom-right panel.
  panels <- lapply(seq_along(panels), function(i) {
    if (i == length(panels)) {
      panels[[i]] + theme(legend.position = c(0.78, 0.18),
                          legend.background =
                            element_rect(fill = alpha("white", 0.6),
                                         colour = NA),
                          legend.key.size = unit(0.3, "cm"),
                          legend.title = element_text(size = 7),
                          legend.text  = element_text(size = 7))
    } else {
      panels[[i]] + theme(legend.position = "none")
    }
  })

  edf <- vapply(mods, function(m) sprintf("%.1f", sum(m$edf)), character(1L))
  subtitle <- sprintf(
    "edf sum: gamma=1.4 -> %s, gamma=6.4 -> %s | n_pos = %d | n_rows = %d",
    edf[["1.4"]], edf[["6.4"]], n_pos, nrow(ref_mod$model)
  )

  combined <- wrap_plots(panels, ncol = 4L) +
    plot_annotation(
      title    = paste(sp, "- partial effects at gamma 1.4 vs 6.4"),
      subtitle = subtitle
    )

  out_f <- file.path(OUTPUT_DIR, sprintf("smooths_%s.png", safe_name(sp)))
  ggsave(out_f, combined,
         width = 12, height = 2.5 * ceiling(length(panels) / 4L),
         dpi = 110)
  message("    Wrote ", out_f)
}

# ── Sweep ───────────────────────────────────────────────────────────────────

excluded <- character(0)
t_start  <- Sys.time()

for (sp in species) {
  message(sprintf("\n== [%s] ============================================", sp))
  t_sp <- Sys.time()

  prepared <- tryCatch(build_species_data(sp),
                       error = function(e) {
                         message("  Skipping: ", conditionMessage(e))
                         excluded <<- c(excluded, sp)
                         NULL
                       })
  if (is.null(prepared)) next

  message(sprintf("  n_pos = %d, n_rows = %d, prep %.1fs",
                  prepared$n_pos, nrow(prepared$df),
                  as.numeric(difftime(Sys.time(), t_sp, units = "secs"))))

  rows  <- list()
  mods  <- list()
  fail  <- FALSE

  for (g in GAMMAS) {
    message(sprintf("  gamma = %.2f ...", g))
    t_g <- Sys.time()

    # CV
    cv <- tryCatch(
      evaluate_model_cv(prepared$df,
                        formula  = prepared$formula,
                        k        = K_FOLDS,
                        fold_ids = prepared$fold_ids,
                        gamma    = g),
      error = function(e) { message("    CV failed: ", conditionMessage(e)); NULL }
    )
    if (is.null(cv)) { fail <- TRUE; next }

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

    # Full-data fit for smooth plot
    mods[[as.character(g)]] <- tryCatch(
      fit_full(prepared$df, prepared$formula, g),
      error = function(e) {
        message("    Full fit failed: ", conditionMessage(e)); NULL
      }
    )

    message(sprintf("    done in %.1fs  Spearman=%.3f DevExpl=%.3f",
                    as.numeric(difftime(Sys.time(), t_g, units = "secs")),
                    unname(s["spearman_r"]), unname(s["holdout_dev_expl"])))
  }

  if (fail) {
    message("  At least one gamma failed - not writing partial results for ", sp)
    next
  }

  append_rows(rows)
  if (length(mods) == length(GAMMAS))
    save_smooth_panel(mods, sp, prepared$n_pos)

  # Free memory before next species (models are large).
  rm(mods, prepared); gc(full = TRUE)

  message(sprintf("  [%s] total %.1f min   elapsed %.1f min", sp,
                  as.numeric(difftime(Sys.time(), t_sp, units = "mins")),
                  as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

message(sprintf("\nDone. Total elapsed: %.1f min",
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
if (length(excluded))
  message("Excluded: ", paste(excluded, collapse = ", "))
