# plot_gamma_smooths.R
#
# Fit full-data GAMs at gamma = {1.4, 6.4} for the smoke species and overlay
# partial-effect curves per smooth term. Lets us see whether the BIC-like
# penalty is actually flattening any smooths or just leaving everything alone.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(ggplot2)
  library(patchwork)
})

GAMMAS         <- c(1.4, 6.4)
SPECIES        <- c("Superb Fairywren", "Rose Robin")
CACHE_DIR      <- "ebirdabund_cache_nsw_buffer"
OUTPUT_DIR     <- "gamma_sweep/smoke"
SEED           <- 42L
HEX_SPACING_KM <- 5
MAX_COUNT      <- 200L

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading sampling master ...")
master <- readRDS(file.path(CACHE_DIR, "sampling_master.rds"))

build_species_data <- function(sp) {
  zf_f <- file.path(CACHE_DIR, sprintf("zerofilled_%s.rds", safe_name(sp)))
  obs  <- readRDS(zf_f)[, c("checklist_id", "observation_count"), drop = FALSE]
  zf   <- merge(master, obs, by = "checklist_id", all.x = TRUE)
  zf$observation_count[is.na(zf$observation_count)] <- "0"
  zf$species_observed <- zf$observation_count != "0"
  df <- clean_ebird(zf, max_count = MAX_COUNT)
  set.seed(SEED)
  subsample_hex(df, spacing_km = HEX_SPACING_KM)
}

fit_at_gamma <- function(df, g) {
  hab_cols <- detect_hab_cols(df)
  formula  <- build_gam_formula(df, hab_cols)
  mgcv::bam(stats::as.formula(formula),
            data     = df,
            family   = mgcv::nb(),
            method   = "fREML",
            discrete = TRUE,
            gamma    = g,
            select   = TRUE)
}

# Build an overlay plot per smooth term: predict at varying x with other
# predictors held at training-data means, one line per gamma.
smooth_overlay <- function(mods, vname) {
  ref <- mods[[1]]
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
    theme_minimal(base_size = 9) +
    theme(legend.position = "none")
}

for (sp in SPECIES) {
  message(sprintf("\n== [%s] ============================================", sp))
  df <- build_species_data(sp)
  message(sprintf("  n = %d rows, fitting at gammas %s ...",
                  nrow(df), paste(GAMMAS, collapse = ", ")))

  mods <- list()
  for (g in GAMMAS) {
    t0 <- Sys.time()
    mods[[as.character(g)]] <- fit_at_gamma(df, g)
    message(sprintf("    gamma = %.1f done in %.1fs (edf sum = %.1f)",
                    g,
                    as.numeric(difftime(Sys.time(), t0, units = "secs")),
                    sum(mods[[as.character(g)]]$edf)))
  }

  ref_mod    <- mods[[1]]
  term_names <- vapply(ref_mod$smooth, function(s) s$term, character(1L))

  panels <- lapply(term_names, function(v) smooth_overlay(mods, v))
  panels <- Filter(Negate(is.null), panels)

  legend_plot <- ggplot(data.frame(g = names(mods), v = 1)) +
    geom_line(aes(x = g, y = v, colour = g, group = 1), linewidth = 1) +
    scale_colour_manual("gamma",
                        values = c("1.4" = "steelblue", "6.4" = "firebrick")) +
    theme_void(base_size = 9) +
    theme(legend.position = "right")
  shared_legend <- cowplot::get_legend(legend_plot +
                                         theme(legend.position = "right"))

  edf_summary <- vapply(mods,
                        function(m) sprintf("%.1f", sum(m$edf)), character(1L))
  subtitle <- sprintf("edf sum: gamma=1.4 -> %s, gamma=6.4 -> %s | n=%d",
                      edf_summary[["1.4"]], edf_summary[["6.4"]], nrow(df))

  combined <- wrap_plots(panels, ncol = 4L) +
    plot_annotation(
      title    = paste(sp, "- partial effects at gamma 1.4 vs 6.4"),
      subtitle = subtitle
    )

  out_f <- file.path(OUTPUT_DIR,
                     sprintf("smooths_%s.png", safe_name(sp)))
  ggsave(out_f, combined,
         width = 12, height = 2.5 * ceiling(length(panels) / 4L),
         dpi = 110)
  message("  Wrote ", out_f)
}

message("\nDone.")
