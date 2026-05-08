#' Estimate relative abundance for multiple species in parallel
#'
#' Runs \code{\link{estimate_abundance}} for each species in
#' \code{species_list} in parallel.  Before launching the full batch, it runs
#' the first species sequentially as a timing benchmark and prints an estimated
#' wall-time for the remaining species.
#'
#' It is strongly recommended to call \code{\link{prepare_covariates}} first
#' and pass the result via \code{cov_stack}, otherwise each worker will
#' re-download the raster stack.
#'
#' @param polygon An \code{sf} object (POLYGON or MULTIPOLYGON) defining the
#'   study area.
#' @param ebird_zip Path to the raw eBird EBD \code{.zip} or \code{.txt} file.
#' @param sampling_txt Path to the eBird sampling-events \code{.txt} file.
#' @param species_list Character vector of species common names exactly as they
#'   appear in eBird (duplicates are silently removed).
#' @param cov_stack A \code{terra::SpatRaster} returned by
#'   \code{\link{prepare_covariates}}.  If \code{NULL}, each worker calls
#'   \code{prepare_covariates} independently.
#' @param cache_dir Cache directory for zero-filled eBird data.  Default
#'   \code{"ebirdabund_cache"}.
#' @param grid_res_km Prediction grid resolution in km.  Default \code{1}.
#' @param hex_spacing_km Hex-cell diameter for spatiotemporal subsampling.
#'   Default \code{5}.
#' @param peak_time Optional decimal-hour override for prediction start time.
#' @param use_range If \code{TRUE} (default), cells outside the eBird species
#'   range are set to \code{NA}.  See \code{\link{predict_species_map}} for
#'   details.
#' @param range_resolution Resolution passed to \code{ebirdst::load_ranges};
#'   \code{"27km"} (default) or \code{"9km"}.
#' @param n_cores Number of parallel workers.  Defaults to
#'   \code{parallel::detectCores() - 1} (minimum 1).
#' @param output_dir If non-\code{NULL}, each species' abundance map
#'   (\code{.png}) and prediction raster (\code{.tif}) are saved here using
#'   the species name as the file stem.
#'
#' @return A named list with one element per species (in the same order as
#'   \code{species_list}).  Each element is either the list returned by
#'   \code{\link{estimate_abundance}} or a \code{simpleError} if that species
#'   failed, so a single failure never aborts the batch.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(ebirdabund)
#'
#' nsw <- sf::st_as_sf(geodata::gadm("AUS", level = 1,
#'                                    path = "ebirdabund_cache")[
#'          geodata::gadm("AUS", level = 1, path = "ebirdabund_cache")$NAME_1
#'            == "New South Wales", ])
#'
#' cov <- prepare_covariates(nsw)
#'
#' results <- estimate_abundance_batch(
#'   polygon      = nsw,
#'   ebird_zip    = "raw_data/ebd_AU-NSW.txt",
#'   sampling_txt = "raw_data/ebd_AU-NSW_sampling.txt",
#'   species_list = c("Superb Fairywren", "Superb Parrot", "Gang-gang Cockatoo"),
#'   cov_stack    = cov,
#'   output_dir   = "species_maps"
#' )
#'
#' # Check which species succeeded
#' ok <- !vapply(results, inherits, TRUE, "error")
#' message(sum(ok), " of ", length(results), " species succeeded.")
#' }
#'
#' @export
estimate_abundance_batch <- function(
    polygon,
    ebird_zip,
    sampling_txt,
    species_list,
    taxonomy         = NULL,
    cov_stack        = NULL,
    cache_dir        = "ebirdabund_cache",
    grid_res_km      = 1,
    hex_spacing_km   = 5,
    peak_time        = NULL,
    use_range        = TRUE,
    botw_path        = NULL,
    range_resolution = "27km",
    border           = NULL,
    n_cores          = max(1L, parallel::detectCores() - 1L),
    output_dir       = NULL,
    log_file         = NULL) {

  # ── Validation ───────────────────────────────────────────────────────────────
  if (!is.character(species_list) || length(species_list) == 0L) {
    stop("`species_list` must be a non-empty character vector.")
  }
  species_list <- unique(trimws(species_list))
  n_sp <- length(species_list)
  n_cores <- min(max(1L, as.integer(n_cores)), parallel::detectCores())

  # grid_res_km may be a vector; each value gets its own subdirectory
  grid_res_km <- sort(unique(as.numeric(grid_res_km)))

  if (!is.null(output_dir)) {
    for (res_km in grid_res_km) {
      dir.create(file.path(output_dir, paste0(res_km, "km")),
                 showWarnings = FALSE, recursive = TRUE)
    }
    dir.create(file.path(output_dir, "models"),  showWarnings = FALSE)
    dir.create(file.path(output_dir, "smooths"), showWarnings = FALSE)
  }

  # ── Helpers ──────────────────────────────────────────────────────────────────
  fmt_duration <- function(s) {
    s <- round(s)
    if (s < 60)   return(sprintf("%ds", s))
    if (s < 3600) return(sprintf("%dm %02ds", s %/% 60L, s %% 60L))
    return(sprintf("%dh %02dm %02ds", s %/% 3600L, (s %% 3600L) %/% 60L, s %% 60L))
  }

  # Fit model once, predict + save at every requested resolution.
  # Returns a lightweight summary list (no rasters over the socket).
  run_species <- function(sp, cov) {
    model_fit <- fit_species_model(
      polygon        = polygon,
      ebird_zip      = ebird_zip,
      sampling_txt   = sampling_txt,
      species        = sp,
      cov_stack      = cov,
      cache_dir      = cache_dir,
      hex_spacing_km = hex_spacing_km
    )

    if (!is.null(output_dir)) {
      saveRDS(model_fit$model,
              file.path(output_dir, "models", paste0(safe_name(sp), ".rds")))
      smooth_plot <- plot_gam_smooths(model_fit$model, sp)
      if (!is.null(smooth_plot)) {
        n_smooths <- length(model_fit$model$smooth)
        ggplot2::ggsave(
          file.path(output_dir, "smooths", paste0(safe_name(sp), ".png")),
          smooth_plot,
          width  = 14,
          height = ceiling(n_smooths / 4) * 3 + 1.5,
          dpi    = 150
        )
      }
    }

    # ── Per-species stats from training data ─────────────────────────────────
    max_obs_count <- max(model_fit$data$observation_count, na.rm = TRUE)

    model_sum        <- NA_real_
    max_modeled_abd  <- NA_real_
    spatial_cv       <- NA_real_
    range_source_out <- NA_character_
    excluded         <- FALSE
    cached_range_vect <- NULL
    cached_peak_doy   <- NULL   # computed by first resolution, reused by rest
    cached_peak_time  <- NULL

    for (res_km in grid_res_km) {
      pred <- predict_species_map(
        model_fit        = model_fit,
        polygon          = polygon,
        species          = sp,
        sci_name         = sci_lookup[[sp]],
        grid_res_km      = res_km,
        peak_doy         = cached_peak_doy,
        peak_time        = cached_peak_time,
        use_range        = use_range,
        botw_path        = botw_path,
        range_resolution = range_resolution,
        border           = border,
        range_vect       = cached_range_vect
      )
      if (is.null(cached_range_vect) && !is.null(pred$range_vect))
        cached_range_vect <- pred$range_vect
      if (is.null(cached_peak_doy)) {
        cached_peak_doy  <- pred$peak_doy
        cached_peak_time <- pred$peak_time
      }

      # After the first resolution, capture raster-based stats and check
      # whether predicted abundance is meaningful.
      if (res_km == grid_res_km[1L]) {
        abd_vals        <- terra::values(pred$predictions[["abd"]], na.rm = TRUE)
        model_sum       <- sum(abd_vals)
        max_modeled_abd <- if (length(abd_vals) > 0L) max(abd_vals) else NA_real_
        spatial_cv      <- if (length(abd_vals) > 1L && mean(abd_vals) > 0)
          round(sd(abd_vals) / mean(abd_vals), 4) else NA_real_
        range_source_out <- pred$range_source

        if (is.na(model_sum) || model_sum <= 1e-5) {
          excluded <- TRUE
          message(sprintf(
            "  [EXCLUDED] %s: negligible polygon abundance (sum = %.2e)",
            sp, if (is.na(model_sum)) 0 else model_sum
          ))
          break
        }
      }

      if (!is.null(output_dir)) {
        out_dir <- file.path(output_dir, paste0(res_km, "km"))
        stem    <- file.path(out_dir, safe_name(sp))
        ggplot2::ggsave(paste0(stem, ".png"), pred$plot,
                        width = 10, height = 9, dpi = 150)
        terra::writeRaster(pred$predictions, paste0(stem, ".tif"),
                           overwrite = TRUE)
      }
    }

    dev_expl <- withCallingHandlers(
      tryCatch(summary(model_fit$model)$dev.expl, error = function(e) NA_real_),
      warning = function(w) invokeRestart("muffleWarning")
    )
    if (is.null(dev_expl) || is.nan(dev_expl) || is.na(dev_expl)) {
      dev_expl <- tryCatch({
        mod   <- model_fit$model
        theta <- exp(mod$family$getTheta())
        y     <- mod$y
        fv    <- mod$fitted.values
        dr    <- 2 * (ifelse(y > 0, y * log(y / fv), 0) -
                      (y + theta) * log((y + theta) / (fv + theta)))
        de <- 1 - sum(dr, na.rm = TRUE) / mod$null.deviance
        if (is.nan(de) || is.na(de)) NA_real_ else de
      }, error = function(e) NA_real_)
    }

    list(n_checklists     = nrow(model_fit$data),
         n_positive       = sum(model_fit$data$observation_count > 0L),
         dev_expl         = dev_expl,
         model_sum        = model_sum,
         peak_doy         = cached_peak_doy,
         peak_time        = cached_peak_time,
         max_obs_count    = max_obs_count,
         max_modeled_abd  = max_modeled_abd,
         range_source     = range_source_out,
         spatial_cv       = spatial_cv,
         excluded         = excluded,
         exclusion_reason = if (excluded) "negligible_abundance" else NA_character_)
  }

  # Named vector: common_name -> scientific_name (NA if unknown)
  sci_lookup <- if (!is.null(taxonomy)) {
    setNames(taxonomy$scientific_name, taxonomy$common_name)
  } else {
    setNames(rep(NA_character_, length(species_list)), species_list)
  }

  # Format one result as a single-row data.frame for the log.
  make_log_row <- function(sp, res) {
    base <- data.frame(
      common_name      = sp,
      scientific_name  = unname(sci_lookup[sp]),
      run_date         = format(Sys.Date(), "%Y-%m-%d"),
      stringsAsFactors = FALSE
    )
    if (inherits(res, "error")) {
      cbind(base, data.frame(
        status           = "failed",
        exclusion_reason = NA_character_,
        n_checklists     = NA_integer_,
        n_positive       = NA_integer_,
        dev_expl         = NA_real_,
        model_sum        = NA_real_,
        peak_doy         = NA_integer_,
        peak_time        = NA_real_,
        max_obs_count    = NA_integer_,
        max_modeled_abd  = NA_real_,
        range_source     = NA_character_,
        spatial_cv       = NA_real_,
        error_message    = conditionMessage(res),
        stringsAsFactors = FALSE
      ))
    } else {
      cbind(base, data.frame(
        status           = if (isTRUE(res$excluded)) "excluded" else "ok",
        exclusion_reason = res$exclusion_reason,
        n_checklists     = res$n_checklists,
        n_positive       = res$n_positive,
        dev_expl         = res$dev_expl,
        model_sum        = res$model_sum,
        peak_doy         = res$peak_doy,
        peak_time        = res$peak_time,
        max_obs_count    = res$max_obs_count,
        max_modeled_abd  = res$max_modeled_abd,
        range_source     = res$range_source,
        spatial_cv       = res$spatial_cv,
        error_message    = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Append one row to the log CSV; write header only on first write.
  write_log_row <- function(sp, res) {
    if (is.null(log_file)) return(invisible(NULL))
    row <- make_log_row(sp, res)
    if (!file.exists(log_file)) {
      write.csv(row, log_file, row.names = FALSE)
    } else {
      write.table(row, log_file, sep = ",", col.names = FALSE,
                  row.names = FALSE, append = TRUE, quote = TRUE)
    }
  }

  run_one <- function(species, cov) {
    tryCatch(
      run_species(species, cov),
      ebirdabund_excluded = function(e) {
        message(sprintf("  [EXCLUDED] %s: %s", species, conditionMessage(e)))
        list(excluded         = TRUE,
             exclusion_reason = e$exclusion_reason,
             n_checklists     = e$n_checklists,
             n_positive       = e$n_positive,
             dev_expl         = NA_real_,
             model_sum        = NA_real_,
             peak_doy         = NA_integer_,
             peak_time        = NA_real_,
             max_obs_count    = NA_integer_,
             max_modeled_abd  = NA_real_,
             range_source     = NA_character_,
             spatial_cv       = NA_real_)
      },
      error = function(e) {
        message(sprintf("  [FAILED] %s: %s", species, conditionMessage(e)))
        e
      }
    )
  }

  results <- vector("list", n_sp)
  names(results) <- species_list

  # ── Timing benchmark: run species 1 serially ─────────────────────────────────
  message(sprintf(
    "\n── Timing benchmark [1/%d]: %s", n_sp, species_list[1]
  ))
  t_batch_start <- proc.time()[["elapsed"]]

  results[[1]] <- run_one(species_list[1], cov_stack)
  write_log_row(species_list[1], results[[1]])

  t_single <- proc.time()[["elapsed"]] - t_batch_start

  # ── Time estimate for remaining species ──────────────────────────────────────
  n_remaining <- n_sp - 1L

  if (n_remaining > 0L) {
    workers   <- min(n_cores, n_remaining)
    n_batches <- ceiling(n_remaining / workers)
    est_remain <- n_batches * t_single
    est_total  <- t_single + est_remain

    message(sprintf(paste(
      "\n┌─ Time estimate ──────────────────────────────────────────────┐",
      "│  Per-species sample time : %-34s│",
      "│  Remaining species       : %-34s│",
      "│  Parallel workers        : %-34s│",
      "│  Estimated batches       : %-34s│",
      "│  Est. time remaining     : %-34s│",
      "│  Est. total wall time    : %-34s│",
      "└──────────────────────────────────────────────────────────────┘",
      sep = "\n"
    ),
    fmt_duration(t_single),
    n_remaining,
    workers,
    n_batches,
    fmt_duration(est_remain),
    fmt_duration(est_total)
    ))
  }

  if (n_sp == 1L) return(results)

  # ── Parallel batch ────────────────────────────────────────────────────────────
  workers <- min(n_cores, n_remaining)
  message(sprintf(
    "\n── Parallel batch: %d species on %d worker(s) ──────────────────",
    n_remaining, workers
  ))

  # terra SpatRasters must be wrapped to survive PSOCK serialisation
  wrapped_cov <- if (!is.null(cov_stack)) terra::wrap(cov_stack) else NULL

  # Workers load the package; support both installed and source-dir dev workflows
  pkg_src <- normalizePath("ebirdabund", mustWork = FALSE)

  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, "pkg_src", envir = environment())
  parallel::clusterEvalQ(cl, {
    if (requireNamespace("ebirdabund", quietly = TRUE)) {
      library(ebirdabund)
    } else if (dir.exists(pkg_src)) {
      devtools::load_all(pkg_src, quiet = TRUE)
    } else {
      stop("Cannot load ebirdabund on worker: install it or run from the package source directory.")
    }
  })

  parallel::clusterExport(
    cl,
    c("polygon", "ebird_zip", "sampling_txt", "wrapped_cov",
      "cache_dir", "grid_res_km", "hex_spacing_km", "peak_time",
      "use_range", "botw_path", "range_resolution", "border", "output_dir",
      "sci_lookup", "run_species"),
    envir = environment()
  )

  # Workers fit once and predict at all resolutions, saving to disk immediately.
  # Only a lightweight summary is returned over the socket (no raster transfer).
  # Process in chunks of `workers` so results are logged after each chunk —
  # a crash loses at most one chunk's worth of records.
  remaining_spp <- species_list[-1]
  chunks <- split(remaining_spp,
                  ceiling(seq_along(remaining_spp) / workers))

  for (chunk in chunks) {
    chunk_res <- tryCatch(
      parallel::parLapplyLB(
        cl,
        chunk,
        function(sp) {
          cov <- if (!is.null(wrapped_cov)) terra::unwrap(wrapped_cov) else NULL
          tryCatch(run_species(sp, cov), error = function(e) e)
        }
      ),
      error = function(e) {
        # Worker process died (e.g. SIGPIPE, OOM, GDAL segfault). Rebuild the
        # cluster so subsequent chunks can continue, and treat every species in
        # this chunk as failed.
        message(sprintf(
          "  [CLUSTER ERROR] chunk failed (%s). Rebuilding cluster and marking %d species as failed.",
          conditionMessage(e), length(chunk)
        ))
        tryCatch(parallel::stopCluster(cl), error = function(e2) invisible(NULL))
        cl <<- parallel::makeCluster(workers)
        parallel::clusterExport(cl, "pkg_src", envir = environment())
        parallel::clusterEvalQ(cl, {
          if (requireNamespace("ebirdabund", quietly = TRUE)) {
            library(ebirdabund)
          } else if (dir.exists(pkg_src)) {
            devtools::load_all(pkg_src, quiet = TRUE)
          }
        })
        parallel::clusterExport(
          cl,
          c("polygon", "ebird_zip", "sampling_txt", "wrapped_cov",
            "cache_dir", "grid_res_km", "hex_spacing_km", "peak_time",
            "use_range", "botw_path", "range_resolution", "border", "output_dir",
            "sci_lookup", "run_species"),
          envir = parent.env(environment())
        )
        lapply(chunk, function(sp) simpleError(conditionMessage(e)))
      }
    )
    for (i in seq_along(chunk)) {
      sp  <- chunk[i]
      res <- chunk_res[[i]]
      results[[sp]] <- res
      write_log_row(sp, res)
    }
  }

  # ── Summary ───────────────────────────────────────────────────────────────────
  t_total <- proc.time()[["elapsed"]] - t_batch_start
  n_ok  <- sum(!vapply(results, inherits, TRUE, "error"))
  n_err <- n_sp - n_ok

  message(sprintf(
    "\n── Batch complete: %d/%d succeeded, %d failed, wall time %s ────",
    n_ok, n_sp, n_err, fmt_duration(t_total)
  ))
  if (n_err > 0L) {
    failed <- species_list[vapply(results, inherits, TRUE, "error")]
    message("  Failed species: ", paste(failed, collapse = ", "))
  }

  results
}
