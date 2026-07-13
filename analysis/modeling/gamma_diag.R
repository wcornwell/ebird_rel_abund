# gamma_diag.R
#
# Consolidated gamma / convergence diagnostics for the production NB GAM fit.
# One script, four modes (select with the first CLI arg; default "rare").
# Every mode reuses the same data-prep and fit-capture harness so the fit path
# is guaranteed identical to production — in particular the protocol reference
# level always comes from modal_protocol(), never a hardcoded name.
#
#   rare  — do the 10 rarest species that fitted at gamma=1.4 still converge at
#           the BIC-like gammas? Tests 1.4, log(n_pos)/2, log(nrow)/2 per
#           species.                          -> gamma_sweep/rare_diag.csv
#   g2    — do the species that broke at gamma=log(nrow)/2 converge at a fixed
#           gamma=2?                          -> gamma_sweep/g2_diag.csv
#   maxit — single-species maxit sweep (default Plum-headed Finch) at
#           gamma=log(nrow)/2: is it converging but running out of iterations?
#                                             -> gamma_sweep/maxit_diag.csv
#   error — verbose single-species fit at gamma=log(nrow)/2: print the formula,
#           per-covariate NA/unique/range table, and the actual error/warnings.
#           (console only, no CSV)
#
# Usage:  Rscript analysis/modeling/gamma_diag.R [mode] [species ...]
#         The optional trailing species names override the mode's default set,
#         e.g.  Rscript analysis/modeling/gamma_diag.R error "Superb Fruit-Dove"

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
})

CACHE_DIR      <- "ebirdabund_cache_nsw_buffer"
SEED           <- 42L
HEX_SPACING_KM <- 5
MAX_COUNT      <- 200L
TIME_KNOTS     <- list(time_observations_started = c(0, 24),
                       day_of_year               = c(0, 365))

# --- Mode configuration --------------------------------------------------
# gammas: function(n_pos, n_row) -> named numeric vector of gamma values.
# maxit : integer vector of control$maxit values (NULL = bam default only).
# output: CSV path, or NULL for console-only (verbose) diagnostics.
MODES <- list(
  rare = list(
    species = c("Hall's Babbler", "Inland Dotterel", "Common Ostrich",
                "Flock Bronzewing", "Australian Painted-Snipe",
                "Black-bellied Plover", "Pectoral Sandpiper",
                "Superb Fruit-Dove", "Black-breasted Kite",
                "Black-faced Cormorant"),
    gammas  = function(n_pos, n_row) c(g_1.4      = 1.4,
                                       g_log_npos = log(n_pos) / 2,
                                       g_log_nrow = log(n_row) / 2),
    maxit   = NULL,
    output  = "gamma_sweep/rare_diag.csv",
    verbose = FALSE
  ),
  g2 = list(
    species = c("Plum-headed Finch", "Superb Fruit-Dove"),
    gammas  = function(n_pos, n_row) c(g_2 = 2.0),
    maxit   = NULL,
    output  = "gamma_sweep/g2_diag.csv",
    verbose = FALSE
  ),
  maxit = list(
    species = "Plum-headed Finch",
    gammas  = function(n_pos, n_row) c(g_log_nrow = log(n_row) / 2),
    maxit   = c(200L, 500L, 1000L, 2000L),
    output  = "gamma_sweep/maxit_diag.csv",
    verbose = FALSE
  ),
  error = list(
    species = "Plum-headed Finch",
    gammas  = function(n_pos, n_row) c(g_log_nrow = log(n_row) / 2),
    maxit   = NULL,
    output  = NULL,
    verbose = TRUE
  )
)

# --- Shared harness ------------------------------------------------------

# Rebuild a species' subsampled model frame from the shared sampling master
# and its zerofilled cache — mirrors the production zero-fill exactly.
build_species_data <- function(sp, master) {
  zf_f <- file.path(CACHE_DIR, sprintf("zerofilled_%s.rds", safe_name(sp)))
  if (!file.exists(zf_f)) stop("zerofilled cache not found: ", zf_f)
  obs <- readRDS(zf_f)[, c("checklist_id", "observation_count"), drop = FALSE]
  zf  <- merge(master, obs, by = "checklist_id", all.x = TRUE)
  zf$observation_count[is.na(zf$observation_count)] <- "0"
  zf$species_observed <- zf$observation_count != "0"
  df <- clean_ebird(zf, max_count = MAX_COUNT)
  set.seed(SEED)
  subsample_hex(df, spacing_km = HEX_SPACING_KM)
}

# Habitat-column selection + formula + protocol relevel, matching production.
prep_fit <- function(df) {
  hab_cols <- grep(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height|nightlights|palsar_hv)",
    names(df), value = TRUE
  )
  hab_cols <- setdiff(hab_cols, "lc_shrubs")
  hab_cols <- hab_cols[vapply(hab_cols, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]
  df$protocol_type <- stats::relevel(df$protocol_type,
                                     ref = modal_protocol(df$protocol_type))
  list(df = df, formula = build_gam_formula(df, hab_cols), hab_cols = hab_cols)
}

# One discrete=TRUE select=TRUE fit. Warnings are logged and muffled so the
# fit runs to completion; a "did not converge" warning flags non-convergence;
# errors are caught and their message recorded.
fit_capture <- function(df, formula, gamma_val, maxit_val = NULL) {
  ctrl      <- if (is.null(maxit_val)) list() else list(maxit = maxit_val)
  warns     <- character(0)
  converged <- TRUE
  err_msg   <- NA_character_
  t0 <- Sys.time()
  mod <- withCallingHandlers(
    tryCatch(
      mgcv::bam(formula, data = df, family = mgcv::nb(),
                method = "fREML", discrete = TRUE, knots = TIME_KNOTS,
                gamma = gamma_val, select = TRUE, control = ctrl),
      error = function(e) {
        err_msg   <<- conditionMessage(e)
        converged <<- FALSE
        NULL
      }
    ),
    warning = function(w) {
      msg   <- conditionMessage(w)
      warns <<- c(warns, msg)
      if (grepl("did not converge", msg, fixed = TRUE)) converged <<- FALSE
      invokeRestart("muffleWarning")
    }
  )
  list(
    converged = converged && !is.null(mod),
    edf_sum   = if (is.null(mod)) NA_real_ else sum(mod$edf),
    elapsed_s = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    warnings  = paste(unique(warns), collapse = " | "),
    error_msg = err_msg
  )
}

# Verbose per-covariate sanity table (error mode).
print_covariate_table <- function(df, hab_cols) {
  message("\nNA counts per covariate column:")
  cols <- c(hab_cols, "day_of_year", "time_observations_started",
            "duration_minutes", "effort_distance_km", "number_observers",
            "observer_expertise", "observation_count")
  for (col in cols) {
    if (!col %in% names(df)) next
    message(sprintf("  %-30s NAs=%d  unique=%d  range=[%g, %g]",
                    col, sum(is.na(df[[col]])),
                    length(unique(stats::na.omit(df[[col]]))),
                    suppressWarnings(min(df[[col]], na.rm = TRUE)),
                    suppressWarnings(max(df[[col]], na.rm = TRUE))))
  }
}

# --- Dispatch ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1L]] else "rare"
if (!mode %in% names(MODES))
  stop("unknown mode '", mode, "' — choose one of: ",
       paste(names(MODES), collapse = ", "))
cfg     <- MODES[[mode]]
species <- if (length(args) >= 2L) args[-1L] else cfg$species

message(sprintf("Mode: %s   species: %d", mode, length(species)))
message("Loading sampling master ...")
master <- readRDS(file.path(CACHE_DIR, "sampling_master.rds"))

if (!is.null(cfg$output))
  dir.create(dirname(cfg$output), showWarnings = FALSE, recursive = TRUE)

results <- list()
t_start <- Sys.time()

for (sp in species) {
  message(sprintf("\n== %s ==", sp))
  df0 <- tryCatch(build_species_data(sp, master), error = function(e) {
    message("  Skipping (data prep failed): ", conditionMessage(e)); NULL
  })
  if (is.null(df0)) next

  n_pos <- sum(df0$observation_count > 0L)
  n_row <- nrow(df0)
  message(sprintf("  n_pos = %d   n_rows = %d", n_pos, n_row))

  p       <- prep_fit(df0)
  gammas  <- cfg$gammas(n_pos, n_row)
  maxits  <- if (is.null(cfg$maxit)) NA_integer_ else cfg$maxit

  if (isTRUE(cfg$verbose)) {
    message("\nFormula:")
    message(deparse(p$formula, width.cutoff = 80))
    print_covariate_table(p$df, p$hab_cols)
    message("\nFitting (errors/warnings will be printed) ...")
  }

  for (nm in names(gammas)) {
    g <- gammas[[nm]]
    for (m in maxits) {
      maxit_val <- if (is.na(m)) NULL else m
      label <- if (is.na(m)) sprintf("%s (gamma=%.2f)", nm, g)
               else          sprintf("%s (gamma=%.2f, maxit=%d)", nm, g, m)
      message(sprintf("  %s ...", label))
      res <- fit_capture(p$df, p$formula, g, maxit_val = maxit_val)
      message(sprintf("    converged=%s  edf=%.1f  time=%.1fs",
                      res$converged,
                      if (is.na(res$edf_sum)) -1 else res$edf_sum,
                      res$elapsed_s))
      if (isTRUE(cfg$verbose)) {
        if (nzchar(res$warnings)) message("    warnings: ", res$warnings)
        if (!is.na(res$error_msg)) message("    error: ", res$error_msg)
      }
      results[[length(results) + 1L]] <- data.frame(
        species     = sp,
        n_pos       = n_pos,
        n_rows      = n_row,
        gamma_label = nm,
        gamma       = g,
        maxit       = if (is.na(m)) NA_integer_ else m,
        converged   = res$converged,
        edf_sum     = res$edf_sum,
        elapsed_s   = res$elapsed_s,
        warnings    = res$warnings,
        error_msg   = res$error_msg,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (!is.null(cfg$output) && length(results) > 0L) {
  out <- do.call(rbind, results)
  utils::write.csv(out, cfg$output, row.names = FALSE)
  message(sprintf("\nWrote %s", cfg$output))

  cat("\n── Convergence summary ───────────────────────────────\n")
  conv_tab <- aggregate(converged ~ gamma_label, data = out,
                        FUN = function(x) sprintf("%d/%d", sum(x), length(x)))
  print(conv_tab, row.names = FALSE)
}

message(sprintf("\nTotal elapsed: %.1f min",
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
