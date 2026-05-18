#' Test whether a covariate swap or addition improves model fit
#'
#' Paired k-fold cross-validation across a fixed set of test species, comparing
#' a baseline GAM (effort + current habitat covariates) to a variant in which
#' one covariate is **swapped** for an alternative source or one or more new
#' covariates are **added**. Reuses the existing `sampling_master.rds` and
#' `zerofilled_*.rds` caches — never invalidates or rewrites them. The variant
#' covariate values are extracted once and cached separately so repeat tests
#' are instant.
#'
#' Designed for fast iteration on candidate covariates (e.g. swapping the
#' tree-height layer). See `analysis/test_tree_height_alternatives.R` for an
#' example call.
#'
#' @param variant_name Short identifier used in cache filenames and output
#'   paths (e.g. `"tree_height_eth_lang"`).
#' @param variant_source Either:
#'   * a character path to a raster file (`.tif`/`.vrt`/etc.) — extracted with
#'     `terra::extract(method = "bilinear")` at each checklist location;
#'   * a function `function(bbox, template)` returning a `terra::SpatRaster`,
#'     where `bbox` is the polygon bbox in EPSG:4326 and `template` is the
#'     existing covariate stack (useful for resolution/CRS reference). Used
#'     for derived layers (ratios, transformations) or on-the-fly downloads.
#' @param change A length-1 list with either `swap = "<col>"` or
#'   `add = "<col>"`:
#'   * `swap`: drop the baseline column with that name and replace it with the
#'     extracted variant. The variant inherits whatever transform the baseline
#'     formula applies (e.g. `log1p(tree_height)`), so the comparison is
#'     apples-to-apples. The extracted column is renamed to match the swapped
#'     name automatically.
#'   * `add`: append the extracted column to the master and add an extra
#'     `s(<col>, k = 4)` term to the variant formula. No transform is applied;
#'     pre-transform inside `variant_source` if needed.
#' @param species Character vector of common names to test. `NULL` (default)
#'   loads the fixed test set from `species_csv`.
#' @param species_csv Path to a CSV with a `common_name` column. Default
#'   `analysis/covariate_test_species.csv`.
#' @param cache_dir Region cache directory holding `sampling_master.rds` and
#'   `zerofilled_*.rds`. Default `"ebirdabund_cache_nsw_buffer"`.
#' @param variant_cache_dir Where the per-variant extracted values are cached.
#'   Defaults to `file.path(cache_dir, "variant_extracts")`.
#' @param output_dir Root for results. One subdirectory is created per
#'   `variant_name`. Default `"covariate_tests"`.
#' @param polygon,cov_stack Required only when `variant_source` is a function;
#'   passed through verbatim.
#' @param k Number of CV folds. Default `5`.
#' @param seed Random seed for fold assignment (paired across baseline and
#'   variant). Default `42`.
#' @param hex_spacing_km Hex-cell size for spatiotemporal subsampling.
#'   Default `5`.
#' @param max_count Mega-flock cap, same as `fit_species_model()`. Default
#'   `200`.
#' @param recompute_extract Force re-extraction of the variant cache even if a
#'   cached file exists.
#'
#' @return Invisibly, a list with:
#' \describe{
#'   \item{per_species}{Data frame: one row per species, columns for each
#'     metric in both baseline and variant plus deltas.}
#'   \item{summary}{Data frame: per-metric paired Wilcoxon p-value, mean delta,
#'     and win/tie/loss counts (variant vs baseline).}
#'   \item{output_dir}{Resolved path to the directory containing CSV / TXT /
#'     PDF artefacts.}
#' }
#'
#' Side-effect files under `{output_dir}/{variant_name}/`:
#' \itemize{
#'   \item `per_species_metrics.csv`
#'   \item `summary.txt`
#'   \item `smooth_comparison.pdf` (baseline vs variant partial-effect plots
#'     for the affected covariate, one page per species)
#'   \item `metrics_full.rds` (raw cv objects for re-analysis)
#' }
#'
#' @examples
#' \dontrun{
#' # Swap tree_height (currently Meta CHMv2) for ETH-Lang 2020
#' test_covariate_variant(
#'   variant_name   = "tree_height_eth_lang",
#'   variant_source = "data/eth_lang_2020_nsw.tif",
#'   change         = list(swap = "tree_height")
#' )
#' }
#'
#' @export
test_covariate_variant <- function(variant_name,
                                   variant_source,
                                   change,
                                   species           = NULL,
                                   species_csv       = "analysis/covariate_test_species.csv",
                                   cache_dir         = "ebirdabund_cache_nsw_buffer",
                                   variant_cache_dir = NULL,
                                   output_dir        = "covariate_tests",
                                   polygon           = NULL,
                                   cov_stack         = NULL,
                                   k                 = 5L,
                                   seed              = 42L,
                                   hex_spacing_km    = 5,
                                   max_count         = 200L,
                                   recompute_extract = FALSE) {

  validate_change(change)
  if (is.null(variant_cache_dir))
    variant_cache_dir <- file.path(cache_dir, "variant_extracts")

  master_f <- file.path(cache_dir, "sampling_master.rds")
  if (!file.exists(master_f))
    stop("sampling_master.rds not found in cache_dir: ", cache_dir,
         "\nRun the main batch pipeline first to build it.")

  if (is.null(species)) {
    if (!file.exists(species_csv))
      stop("species_csv not found: ", species_csv)
    test_set <- utils::read.csv(species_csv, stringsAsFactors = FALSE)
    if (!"common_name" %in% names(test_set))
      stop("species_csv must have a `common_name` column.")
    species <- test_set$common_name
  }

  out_dir <- file.path(output_dir, variant_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n╔══════════════════════════════════════════════════════════════╗")
  message(sprintf("║  Covariate variant test: %-36s║", variant_name))
  message(sprintf("║  Change: %-52s║",
                  paste0(names(change), " = ", unlist(change))))
  message(sprintf("║  Species: %-51d║", length(species)))
  message("╚══════════════════════════════════════════════════════════════╝\n")

  master <- readRDS(master_f)

  # ── Variant extraction (cached) ──────────────────────────────────────────
  variant_extract <- extract_variant_covariate(
    variant_source    = variant_source,
    variant_name      = variant_name,
    master            = master,
    variant_cache_dir = variant_cache_dir,
    polygon           = polygon,
    cov_stack         = cov_stack,
    recompute         = recompute_extract
  )

  # ── Per-species paired CV ────────────────────────────────────────────────
  per_species_rows <- list()
  full_fits        <- list()
  excluded         <- character(0)

  for (sp in species) {
    message(sprintf("\n── [%s] ─────────────────────────────────────────",
                    sp))
    res <- tryCatch(
      paired_cv_one_species(
        species          = sp,
        master           = master,
        variant_extract  = variant_extract,
        change           = change,
        cache_dir        = cache_dir,
        k                = k,
        seed             = seed,
        hex_spacing_km   = hex_spacing_km,
        max_count        = max_count
      ),
      error = function(e) {
        message("  Skipping ", sp, ": ", conditionMessage(e))
        excluded <<- c(excluded, sp)
        NULL
      }
    )
    if (is.null(res)) next

    per_species_rows[[sp]] <- assemble_metric_row(sp, res$cv_base, res$cv_var,
                                                  res$n_pos)
    full_fits[[sp]]        <- list(base = res$fit_base, var = res$fit_var)
  }

  if (length(per_species_rows) == 0L)
    stop("No species produced a comparison. See messages above.")

  per_species <- do.call(rbind, per_species_rows)
  rownames(per_species) <- NULL

  # ── Summary stats ────────────────────────────────────────────────────────
  summary_df <- summarise_deltas(per_species)

  # ── Write artefacts ──────────────────────────────────────────────────────
  per_species_csv <- file.path(out_dir, "per_species_metrics.csv")
  utils::write.csv(per_species, per_species_csv, row.names = FALSE)
  message("\nWrote ", per_species_csv)

  summary_txt <- file.path(out_dir, "summary.txt")
  write_summary_txt(summary_df, per_species, variant_name, change,
                    excluded, summary_txt)
  message("Wrote ", summary_txt)

  metrics_rds <- file.path(out_dir, "metrics_full.rds")
  saveRDS(list(per_species = per_species, summary = summary_df,
               full_fits = full_fits, change = change),
          metrics_rds)
  message("Wrote ", metrics_rds)

  pdf_f <- file.path(out_dir, "smooth_comparison.pdf")
  tryCatch(
    render_smooth_comparison_pdf(full_fits, change, pdf_f),
    error = function(e) message("  Smooth PDF skipped: ", conditionMessage(e))
  )

  invisible(list(per_species = per_species, summary = summary_df,
                 output_dir = out_dir))
}


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

# Validate the `change` argument shape.
validate_change <- function(change) {
  if (!is.list(change) || length(change) != 1L)
    stop("`change` must be a length-1 list: list(swap = '<col>') or ",
         "list(add = '<col>').")
  kind <- names(change)
  if (!kind %in% c("swap", "add"))
    stop("`change` must be named 'swap' or 'add', got: ", kind)
  val <- change[[1L]]
  if (!is.character(val) || length(val) != 1L || nchar(val) == 0L)
    stop("`change[[1]]` must be a single non-empty column name.")
  invisible(TRUE)
}


# Extract variant covariate at checklist locations, with caching.
# Cache file: {variant_cache_dir}/{variant_name}.rds
#   data.frame(checklist_id, <variant_value_col>)
# Cache validity: every checklist_id in `master` must appear in the cache.
# If not, the cache is recomputed.
extract_variant_covariate <- function(variant_source,
                                       variant_name,
                                       master,
                                       variant_cache_dir,
                                       polygon  = NULL,
                                       cov_stack = NULL,
                                       recompute = FALSE) {

  dir.create(variant_cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_f <- file.path(variant_cache_dir, sprintf("%s.rds", variant_name))

  if (file.exists(cache_f) && !recompute) {
    cached <- readRDS(cache_f)
    if (all(master$checklist_id %in% cached$checklist_id)) {
      message("  Variant extract loaded from cache: ", cache_f)
      return(cached[cached$checklist_id %in% master$checklist_id, , drop = FALSE])
    }
    message("  Variant cache stale (missing checklist IDs); recomputing.")
  }

  message("  Extracting variant covariate from ",
          if (is.function(variant_source)) "user function" else variant_source)

  # Materialise the variant raster
  if (is.character(variant_source)) {
    if (!file.exists(variant_source))
      stop("variant_source path does not exist: ", variant_source)
    rast <- terra::rast(variant_source)
  } else if (is.function(variant_source)) {
    if (is.null(polygon))
      stop("`polygon` must be supplied when `variant_source` is a function.")
    bbox <- as.numeric(sf::st_bbox(sf::st_transform(polygon, 4326)))
    rast <- variant_source(bbox = bbox, template = cov_stack)
    if (!inherits(rast, "SpatRaster"))
      stop("variant_source function must return a terra::SpatRaster.")
  } else {
    stop("`variant_source` must be a character path or a function.")
  }

  # Build points in raster CRS so terra::extract doesn't reproject the raster
  pts_wgs <- terra::vect(
    data.frame(x = master$longitude, y = master$latitude),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  pts <- if (terra::crs(rast) != terra::crs(pts_wgs))
    terra::project(pts_wgs, terra::crs(rast)) else pts_wgs

  vals <- terra::extract(rast, pts, method = "bilinear")
  vals <- vals[, -1, drop = FALSE]  # drop ID column

  if (ncol(vals) != 1L)
    stop("v1 supports single-column variants only; got ", ncol(vals),
         " columns from variant_source.")

  # Always store under a fixed internal column name so apply_variant() can
  # find it even when the variant's raster layer is named the same as the
  # swap target (otherwise merge() collides and produces foo.x / foo.y).
  # The original raster layer name is preserved as an attribute for
  # diagnostics.
  out <- data.frame(
    checklist_id = master$checklist_id,
    variant_value = vals[[1L]],
    stringsAsFactors = FALSE
  )
  attr(out, "variant_layer_name") <- names(vals)[1L]

  saveRDS(out, cache_f)
  message(sprintf("  Variant extract cached (%d rows) → %s", nrow(out), cache_f))
  out
}


# Apply the swap/add change to a training data frame.
# Returns the modified data frame, preserving rows.
apply_variant <- function(df, variant_extract, change) {
  kind <- names(change)
  target <- change[[1L]]
  INTERNAL_COL <- "variant_value"
  if (!INTERNAL_COL %in% names(variant_extract))
    stop("variant_extract must contain a '", INTERNAL_COL, "' column.")

  # Join the variant column onto df under its internal name (no collision
  # possible because '__variant_value__' is reserved).
  df_var <- merge(df, variant_extract[, c("checklist_id", INTERNAL_COL),
                                       drop = FALSE],
                  by = "checklist_id", all.x = TRUE)
  df_var <- df_var[match(df$checklist_id, df_var$checklist_id), , drop = FALSE]

  if (kind == "swap") {
    if (!target %in% names(df))
      stop("swap target '", target, "' not present in training data columns.")
    df_var[[target]] <- df_var[[INTERNAL_COL]]
  } else {  # add
    if (target %in% names(df_var) && target != INTERNAL_COL)
      stop("add target '", target,
           "' already exists in training data; use swap instead.")
    df_var[[target]] <- df_var[[INTERNAL_COL]]
  }
  df_var[[INTERNAL_COL]] <- NULL

  # NA handling: drop any rows where the variant column ended up NA
  n_before <- nrow(df_var)
  df_var <- df_var[!is.na(df_var[[target]]), , drop = FALSE]
  if (nrow(df_var) < n_before)
    message(sprintf(
      "  Dropped %d rows with NA in variant column '%s' (%d remaining).",
      n_before - nrow(df_var), target, nrow(df_var)))

  df_var
}


# Run paired CV for one species: baseline vs variant, same folds.
# Returns list(cv_base, cv_var, n_pos, fit_base, fit_var).
paired_cv_one_species <- function(species,
                                  master,
                                  variant_extract,
                                  change,
                                  cache_dir,
                                  k,
                                  seed,
                                  hex_spacing_km,
                                  max_count) {

  zf_f <- file.path(cache_dir, sprintf("zerofilled_%s.rds", safe_name(species)))
  if (!file.exists(zf_f))
    stop("zerofilled cache not found: ", zf_f)

  obs <- readRDS(zf_f)[, c("checklist_id", "observation_count"), drop = FALSE]
  zf  <- merge(master, obs, by = "checklist_id", all.x = TRUE)
  zf$observation_count[is.na(zf$observation_count)] <- "0"
  zf$species_observed <- zf$observation_count != "0"

  df_base <- clean_ebird(zf, max_count = max_count)

  # Subsample with a deterministic seed so baseline and variant operate on
  # the exact same rows (variant is derived by column substitution below).
  set.seed(seed)
  df_base <- subsample_hex(df_base, spacing_km = hex_spacing_km)

  n_pos <- sum(df_base$observation_count > 0L)
  if (n_pos < 50L)
    stop(sprintf("Insufficient detections after subsampling (%d, need >= 50).",
                 n_pos))

  # Apply the variant by column substitution / addition
  df_var <- apply_variant(df_base, variant_extract, change)

  # Pair: keep only rows present in both (apply_variant may drop NA rows)
  shared_ids <- intersect(df_base$checklist_id, df_var$checklist_id)
  df_base <- df_base[df_base$checklist_id %in% shared_ids, , drop = FALSE]
  df_var  <- df_var [df_var$checklist_id  %in% shared_ids, , drop = FALSE]
  df_base <- df_base[match(shared_ids, df_base$checklist_id), , drop = FALSE]
  df_var  <- df_var [match(shared_ids, df_var$checklist_id),  , drop = FALSE]

  # Same fold IDs for both fits → paired CV
  set.seed(seed + 1L)
  fold_ids <- sample(rep_len(seq_len(k), nrow(df_base)))

  hab_cols_base <- detect_hab_cols(df_base)
  formula_base  <- build_gam_formula(df_base, hab_cols_base)

  hab_cols_var <- if (names(change) == "add")
    c(hab_cols_base, change$add) else hab_cols_base
  formula_var  <- build_gam_formula(df_var, hab_cols_var)

  message("  Running baseline CV ...")
  cv_base <- evaluate_model_cv(df_base, formula = formula_base,
                               k = k, fold_ids = fold_ids)
  message("  Running variant  CV ...")
  cv_var  <- evaluate_model_cv(df_var,  formula = formula_var,
                               k = k, fold_ids = fold_ids)

  # Full-data fits for the smooth-comparison PDF
  message("  Fitting full-data models for diagnostics ...")
  fit_base <- safe_full_fit(df_base, formula_base)
  fit_var  <- safe_full_fit(df_var,  formula_var)

  list(cv_base = cv_base, cv_var = cv_var, n_pos = n_pos,
       fit_base = fit_base, fit_var = fit_var)
}


# Detect habitat covariate columns the same way fit_gam() does.
detect_hab_cols <- function(df) {
  hab <- grep(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height)",
    names(df), value = TRUE
  )
  hab <- setdiff(hab, "lc_shrubs")
  hab[vapply(hab, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]
}


# A lighter-weight version of fit_gam() with a single attempt, used only for
# generating smooth-comparison plots. We don't want this to abort the whole
# species if it doesn't converge — return NULL on failure.
safe_full_fit <- function(df, formula) {
  protocols <- levels(df$protocol_type)
  ref_proto <- if ("Traveling Count" %in% protocols) "Traveling Count" else
    protocols[1L]
  df$protocol_type <- stats::relevel(df$protocol_type, ref = ref_proto)
  time_knots <- list(time_observations_started = c(0, 24),
                     day_of_year               = c(0, 365))
  tryCatch(
    mgcv::bam(formula, data = df, family = mgcv::nb(),
              method = "fREML", discrete = TRUE,
              knots = time_knots, gamma = 1.4, select = TRUE),
    error = function(e) NULL
  )
}


# Build a one-row metric frame for a species (baseline + variant + deltas).
assemble_metric_row <- function(species, cv_base, cv_var, n_pos) {
  metrics <- c("spearman_r", "pearson_r", "mae", "rmse", "holdout_dev_expl")
  b <- cv_base$summary[metrics]
  v <- cv_var$summary[metrics]
  row <- data.frame(
    species = species,
    n_pos   = n_pos,
    stringsAsFactors = FALSE
  )
  for (m in metrics) {
    row[[paste0(m, "_base")]]  <- unname(b[m])
    row[[paste0(m, "_var")]]   <- unname(v[m])
    row[[paste0(m, "_delta")]] <- unname(v[m] - b[m])
  }
  row
}


# Per-metric paired Wilcoxon + win/tie/loss counts.
summarise_deltas <- function(per_species) {
  # higher_is_better defines sign convention for "win"
  metrics <- list(
    spearman_r        = TRUE,
    pearson_r         = TRUE,
    mae               = FALSE,
    rmse              = FALSE,
    holdout_dev_expl  = TRUE
  )

  rows <- lapply(names(metrics), function(m) {
    base <- per_species[[paste0(m, "_base")]]
    var  <- per_species[[paste0(m, "_var")]]
    deltas <- var - base
    higher_better <- metrics[[m]]

    # Wilcoxon: alternative = "greater" when higher_is_better, else "less"
    # i.e. test whether variant deltas are systematically in the "good"
    # direction
    alt <- if (higher_better) "greater" else "less"
    p_val <- tryCatch(
      stats::wilcox.test(deltas, alternative = alt, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    if (higher_better) {
      wins   <- sum(deltas > 0, na.rm = TRUE)
      losses <- sum(deltas < 0, na.rm = TRUE)
    } else {
      wins   <- sum(deltas < 0, na.rm = TRUE)
      losses <- sum(deltas > 0, na.rm = TRUE)
    }
    ties <- sum(deltas == 0, na.rm = TRUE)

    data.frame(
      metric        = m,
      higher_better = higher_better,
      mean_delta    = mean(deltas, na.rm = TRUE),
      median_delta  = stats::median(deltas, na.rm = TRUE),
      wilcox_p      = p_val,
      wins          = wins,
      ties          = ties,
      losses        = losses,
      n             = sum(!is.na(deltas)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# Human-readable summary file.
write_summary_txt <- function(summary_df, per_species, variant_name, change,
                              excluded, path) {
  con <- file(path, "w")
  on.exit(close(con))

  writeLines(c(
    sprintf("Covariate variant test: %s", variant_name),
    sprintf("Change: %s = %s", names(change), unlist(change)),
    sprintf("Species evaluated: %d", nrow(per_species)),
    if (length(excluded))
      sprintf("Species excluded: %d (%s)", length(excluded),
              paste(excluded, collapse = ", "))
    else "Species excluded: 0",
    "",
    "Per-metric summary (variant vs baseline):",
    sprintf("  %-18s  %10s  %10s  %10s  %5s  %5s  %5s",
            "metric", "mean_delta", "median_dlt", "wilcox_p",
            "win", "tie", "loss")
  ), con)

  for (i in seq_len(nrow(summary_df))) {
    r <- summary_df[i, ]
    writeLines(sprintf(
      "  %-18s  %10.4f  %10.4f  %10.4f  %5d  %5d  %5d",
      r$metric, r$mean_delta, r$median_delta, r$wilcox_p,
      r$wins, r$ties, r$losses
    ), con)
  }

  writeLines(c("", "Convention: 'win' = variant moved in the better direction",
               "(higher Spearman/Pearson/DevExpl, lower MAE/RMSE).",
               "Wilcox p tests deltas against zero in the better direction."),
             con)
}


# Render baseline vs variant partial-effect plots for the affected covariate.
# One page per species; left = baseline smooth, right = variant smooth.
render_smooth_comparison_pdf <- function(full_fits, change, path) {
  kind   <- names(change)
  target <- change[[1L]]

  ok <- vapply(full_fits, function(f) !is.null(f$base) && !is.null(f$var),
               logical(1))
  if (!any(ok)) {
    message("  No converged full-data fits; skipping smooth PDF.")
    return(invisible(NULL))
  }
  full_fits <- full_fits[ok]

  grDevices::pdf(path, width = 10, height = 5)
  on.exit(grDevices::dev.off())

  for (sp in names(full_fits)) {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    smooth_term <- find_smooth_term(full_fits[[sp]]$base, target)
    plot_one_smooth(full_fits[[sp]]$base, smooth_term,
                    sprintf("%s — baseline", sp))
    smooth_term <- find_smooth_term(full_fits[[sp]]$var, target)
    plot_one_smooth(full_fits[[sp]]$var, smooth_term,
                    sprintf("%s — variant (%s = %s)", sp, kind, target))
  }
  message("Wrote ", path)
}


# Find the position of the smooth term whose s.label mentions `target`.
# Handles wrapped names like s(log1p(tree_height)).
find_smooth_term <- function(mod, target) {
  if (is.null(mod)) return(NA_integer_)
  labs <- vapply(mod$smooth, function(s) s$label, character(1))
  hit  <- grep(target, labs, fixed = TRUE)
  if (length(hit) == 0L) NA_integer_ else hit[1L]
}


plot_one_smooth <- function(mod, idx, title) {
  if (is.null(mod) || is.na(idx)) {
    graphics::plot.new()
    graphics::title(main = title, sub = "(no smooth — fit failed)")
    return(invisible(NULL))
  }
  tryCatch(
    mgcv::plot.gam(mod, select = idx, shade = TRUE, scheme = 1,
                   main = title, rug = TRUE),
    error = function(e) {
      graphics::plot.new()
      graphics::title(main = title, sub = paste("plot error:",
                                                 conditionMessage(e)))
    }
  )
}
