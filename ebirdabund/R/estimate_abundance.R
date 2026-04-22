#' Estimate relative abundance within a polygon from raw eBird data
#'
#' Runs the full pipeline: eBird filtering -> zero-filling -> habitat covariate
#' extraction -> spatiotemporal subsampling -> GAM fitting -> abundance
#' prediction.
#'
#' It is strongly recommended to call \code{\link{prepare_covariates}} first
#' and pass the result via \code{cov_stack}.  This separates the one-time
#' raster download (which can exceed 1 GB) from the per-species modelling.
#'
#' @param polygon An \code{sf} object (POLYGON or MULTIPOLYGON) defining the
#'   study area. Any CRS is accepted; it is reprojected internally.
#' @param ebird_zip Path to the raw eBird EBD \code{.zip} \emph{or}
#'   \code{.txt} file.
#' @param sampling_txt Path to the eBird sampling-events \code{.txt} file.
#' @param species Common name exactly as it appears in eBird
#'   (e.g. \code{"Superb Fairywren"}).
#' @param cov_stack A \code{terra::SpatRaster} returned by
#'   \code{\link{prepare_covariates}}. If \code{NULL} (default),
#'   \code{prepare_covariates} is called automatically using \code{cache_dir}.
#' @param cache_dir Directory for the zero-filled eBird cache. Defaults to
#'   \code{"ebirdabund_cache"} in the working directory.
#' @param grid_res_km Prediction grid resolution in km. Default \code{1}.
#' @param hex_spacing_km Hex-cell diameter for spatiotemporal subsampling in
#'   km. Default \code{5}.
#' @param peak_time Optional decimal-hour override for observation start time
#'   used in predictions. Estimated from the model when \code{NULL}.
#'
#' @return A named list with:
#' \describe{
#'   \item{model}{The fitted \code{mgcv::gam} object.}
#'   \item{predictions}{A \code{terra::SpatRaster} with layers \code{abd}
#'     and \code{abd_se}.}
#'   \item{data}{The subsampled data.frame used for model fitting.}
#'   \item{plot}{A \code{ggplot2} map of log-transformed relative abundance.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(ebirdabund)
#'
#' nsw <- sf::st_as_sf(geodata::gadm("AUS", level = 1,
#'                                    path = "ebirdabund_cache")[
#'          geodata::gadm("AUS", level=1, path="ebirdabund_cache")$NAME_1
#'            == "New South Wales", ])
#'
#' # Step 1 — once per region
#' cov <- prepare_covariates(nsw)
#'
#' # Step 2 — once per species
#' result <- estimate_abundance(
#'   polygon      = nsw,
#'   ebird_zip    = "raw_data/ebd_AU-NSW.txt",
#'   sampling_txt = "raw_data/ebd_AU-NSW_sampling.txt",
#'   species      = "Superb Fairywren",
#'   cov_stack    = cov
#' )
#' result$plot
#' }
#'
#' @export
estimate_abundance <- function(polygon,
                               ebird_zip,
                               sampling_txt,
                               species,
                               cov_stack      = NULL,
                               cache_dir      = "ebirdabund_cache",
                               grid_res_km    = 1,
                               hex_spacing_km = 5,
                               peak_time      = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────
  if (!inherits(polygon, "sf") && !inherits(polygon, "sfc")) {
    stop("`polygon` must be an sf or sfc object.")
  }
  if (!file.exists(ebird_zip) &&
        !file.exists(tools::file_path_sans_ext(ebird_zip))) {
    stop("eBird data file not found: ", ebird_zip)
  }
  if (!file.exists(sampling_txt)) {
    stop("Sampling events file not found: ", sampling_txt)
  }
  if (!is.character(species) || nchar(trimws(species)) == 0) {
    stop("`species` must be a non-empty character string.")
  }

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Step 1: Load and filter eBird data ───────────────────────────────────
  message("\n── Step 1/6: Loading eBird data ─────────────────────────────")
  ebird_df <- load_ebird(
    polygon, ebird_zip, sampling_txt, species, cache_dir
  )

  # ── Step 2: Habitat covariate stack ──────────────────────────────────────
  message("\n── Step 2/6: Habitat covariates ─────────────────────────────")
  if (is.null(cov_stack)) {
    message(
      "  (tip: call prepare_covariates(polygon) first to skip this step)"
    )
    cov_stack <- prepare_covariates(polygon, cache_dir)
  } else {
    message("  Using supplied covariate stack.")
  }

  # ── Step 3: Extract covariates at checklist locations ────────────────────
  message("\n── Step 3/6: Extracting covariates ──────────────────────────")
  ebird_df <- extract_covariates(ebird_df, cov_stack)

  cov_cols <- grep("^(lc_|elevation)", names(ebird_df), value = TRUE)
  n_before_drop <- nrow(ebird_df)
  ebird_df <- tidyr::drop_na(ebird_df, dplyr::all_of(cov_cols))
  n_dropped <- n_before_drop - nrow(ebird_df)
  if (n_dropped > 0) {
    message(sprintf(
      "  Dropped %d checklists with missing covariate values (%d remaining).",
      n_dropped, nrow(ebird_df)
    ))
  }

  # ── Step 4: Spatiotemporal subsampling ───────────────────────────────────
  message("\n── Step 4/6: Spatiotemporal subsampling ─────────────────────")
  ebird_ss <- subsample_hex(ebird_df, spacing_km = hex_spacing_km)

  if (nrow(ebird_ss) < 50) {
    warning(
      "Only ", nrow(ebird_ss), " checklists after subsampling — ",
      "model estimates may be unreliable."
    )
  }

  # ── Step 5: Fit GAM ───────────────────────────────────────────────────────
  message("\n── Step 5/6: Fitting GAM ────────────────────────────────────")
  gam_model <- fit_gam(ebird_ss)
  message(sprintf("Deviance explained: %.1f%%",
                  summary(gam_model)$dev.expl * 100))

  # ── Step 6: Predict ───────────────────────────────────────────────────────
  message("\n── Step 6/6: Predicting abundance ───────────────────────────")
  pred_surface <- make_prediction_surface(polygon, grid_res_km)
  pred_surface <- extract_covariates(pred_surface, cov_stack)
  pred_surface <- tidyr::drop_na(pred_surface, dplyr::all_of(cov_cols))

  r_pred <- predict_abundance(
    gam_model, pred_surface, polygon, grid_res_km, peak_time
  )

  message("\nDone.")

  list(
    model       = gam_model,
    predictions = r_pred,
    data        = ebird_ss,
    plot        = plot_abundance(r_pred, polygon, species)
  )
}
