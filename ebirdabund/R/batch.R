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
    n_cores          = max(1L, parallel::detectCores() - 1L),
    output_dir       = NULL) {

  # ── Validation ───────────────────────────────────────────────────────────────
  if (!is.character(species_list) || length(species_list) == 0L) {
    stop("`species_list` must be a non-empty character vector.")
  }
  species_list <- unique(trimws(species_list))
  n_sp <- length(species_list)
  n_cores <- min(max(1L, as.integer(n_cores)), parallel::detectCores())

  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # ── Helpers ──────────────────────────────────────────────────────────────────
  fmt_duration <- function(s) {
    s <- round(s)
    if (s < 60)   return(sprintf("%ds", s))
    if (s < 3600) return(sprintf("%dm %02ds", s %/% 60L, s %% 60L))
    return(sprintf("%dh %02dm %02ds", s %/% 3600L, (s %% 3600L) %/% 60L, s %% 60L))
  }

  save_species <- function(result, species) {
    if (is.null(output_dir) || inherits(result, "error")) return(invisible(NULL))
    stem <- file.path(output_dir, safe_name(species))
    ggplot2::ggsave(paste0(stem, ".png"), result$plot,
                    width = 10, height = 9, dpi = 150)
    terra::writeRaster(result$predictions, paste0(stem, ".tif"), overwrite = TRUE)
  }

  # Named vector: common_name -> scientific_name (NA if unknown)
  sci_lookup <- if (!is.null(taxonomy)) {
    setNames(taxonomy$scientific_name, taxonomy$common_name)
  } else {
    setNames(rep(NA_character_, length(species_list)), species_list)
  }

  run_one <- function(species, cov) {
    tryCatch(
      estimate_abundance(
        polygon          = polygon,
        ebird_zip        = ebird_zip,
        sampling_txt     = sampling_txt,
        species          = species,
        sci_name         = sci_lookup[[species]],
        cov_stack        = cov,
        cache_dir        = cache_dir,
        grid_res_km      = grid_res_km,
        hex_spacing_km   = hex_spacing_km,
        peak_time        = peak_time,
        use_range        = use_range,
        botw_path        = botw_path,
        range_resolution = range_resolution
      ),
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
  save_species(results[[1]], species_list[1])

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
      "use_range", "botw_path", "range_resolution", "output_dir",
      "sci_lookup"),
    envir = environment()
  )

  # Workers save outputs to disk immediately on completion so progress is
  # visible in output_dir and survives a crash. Only a lightweight summary
  # is returned over the socket (no raster serialisation).
  remaining <- parallel::parLapplyLB(
    cl,
    species_list[-1],
    function(sp) {
      cov <- if (!is.null(wrapped_cov)) terra::unwrap(wrapped_cov) else NULL
      tryCatch(
        {
          res <- estimate_abundance(
            polygon          = polygon,
            ebird_zip        = ebird_zip,
            sampling_txt     = sampling_txt,
            species          = sp,
            sci_name         = sci_lookup[[sp]],
            cov_stack        = cov,
            cache_dir        = cache_dir,
            grid_res_km      = grid_res_km,
            hex_spacing_km   = hex_spacing_km,
            peak_time        = peak_time,
            use_range        = use_range,
            botw_path        = botw_path,
            range_resolution = range_resolution
          )
          if (!is.null(output_dir)) {
            stem <- file.path(output_dir, safe_name(sp))
            ggplot2::ggsave(paste0(stem, ".png"), res$plot,
                            width = 10, height = 9, dpi = 150)
            terra::writeRaster(res$predictions, paste0(stem, ".tif"),
                               overwrite = TRUE)
          }
          list(n_checklists = nrow(res$data),
               dev_expl     = summary(res$model)$dev.expl)
        },
        error = function(e) e
      )
    }
  )

  for (i in seq_along(remaining)) {
    sp  <- species_list[i + 1L]
    res <- remaining[[i]]
    results[[sp]] <- res
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
