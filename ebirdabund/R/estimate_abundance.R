#' Fit the abundance model for a species
#'
#' Loads and filters eBird data, extracts habitat covariates, applies
#' spatiotemporal subsampling, and fits a negative-binomial GAM.  The result
#' can be passed directly to [predict_species_map()].
#'
#' It is strongly recommended to call [prepare_covariates()] first and pass the
#' result via `cov_stack`.  This separates the one-time raster download (which
#' can exceed 1 GB) from the per-species modelling.
#'
#' @param polygon An `sf` object (POLYGON or MULTIPOLYGON) defining the study
#'   area. Any CRS is accepted; it is reprojected internally.
#' @param ebird_zip Path to the raw eBird EBD `.zip` *or* `.txt` file.
#' @param sampling_txt Path to the eBird sampling-events `.txt` file.
#' @param species Common name exactly as it appears in eBird
#'   (e.g. `"Superb Fairywren"`).
#' @param cov_stack A `terra::SpatRaster` returned by [prepare_covariates()].
#'   If `NULL` (default), [prepare_covariates()] is called automatically using
#'   `cache_dir`.
#' @param cache_dir Directory for the zero-filled eBird cache. Defaults to
#'   `"ebirdabund_cache"` in the working directory.
#' @param hex_spacing_km Hex-cell diameter for spatiotemporal subsampling in
#'   km. Default `5`.
#'
#' @return A list with:
#' \describe{
#'   \item{model}{The fitted `mgcv::gam` object.}
#'   \item{data}{The subsampled data.frame used for model fitting.}
#'   \item{cov_stack}{The `terra::SpatRaster` covariate stack (needed by
#'     [predict_species_map()]).}
#'   \item{cov_cols}{Character vector of habitat covariate column names.}
#' }
#'
#' @seealso [predict_species_map()], [estimate_abundance()]
#' @export
fit_species_model <- function(polygon,
                              ebird_zip,
                              sampling_txt,
                              species,
                              cov_stack      = NULL,
                              cache_dir      = "ebirdabund_cache",
                              hex_spacing_km = 5) {

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
  message("\n── Step 1/4: Loading eBird data ─────────────────────────────")
  ebird_df <- load_ebird(polygon, ebird_zip, sampling_txt, species, cache_dir)

  # ── Step 2: Habitat covariate stack ──────────────────────────────────────
  message("\n── Step 2/4: Habitat covariates ─────────────────────────────")
  if (is.null(cov_stack)) {
    message(
      "  (tip: call prepare_covariates(polygon) first to skip this step)"
    )
    cov_stack <- prepare_covariates(polygon, cache_dir)
  } else {
    message("  Using supplied covariate stack.")
  }

  # ── Step 3: Extract covariates at checklist locations ────────────────────
  message("\n── Step 3/4: Extracting covariates ──────────────────────────")
  ebird_df <- extract_covariates(ebird_df, cov_stack)

  cov_cols      <- grep("^(lc_|elevation|precip_|temp_)", names(ebird_df), value = TRUE)
  n_before_drop <- nrow(ebird_df)
  ebird_df      <- tidyr::drop_na(ebird_df, dplyr::all_of(cov_cols))
  n_dropped     <- n_before_drop - nrow(ebird_df)
  if (n_dropped > 0) {
    message(sprintf(
      "  Dropped %d checklists with missing covariate values (%d remaining).",
      n_dropped, nrow(ebird_df)
    ))
  }

  # ── Step 4: Spatiotemporal subsampling + GAM ──────────────────────────────
  message("\n── Step 4/4: Spatiotemporal subsampling ─────────────────────")
  ebird_ss <- subsample_hex(ebird_df, spacing_km = hex_spacing_km)

  if (nrow(ebird_ss) < 50) {
    warning(
      "Only ", nrow(ebird_ss), " checklists after subsampling — ",
      "model estimates may be unreliable."
    )
  }

  message("\n── Fitting GAM ──────────────────────────────────────────────")
  gam_model <- fit_gam(ebird_ss)
  message(sprintf("Deviance explained: %.1f%%",
                  summary(gam_model)$dev.expl * 100))

  list(
    model     = gam_model,
    data      = ebird_ss,
    cov_stack = cov_stack,
    cov_cols  = cov_cols
  )
}


#' Predict abundance from a fitted model
#'
#' Generates a gridded abundance surface from the output of
#' [fit_species_model()], optionally masking cells outside the eBird species
#' range.
#'
#' Range polygons are obtained via `ebirdst::load_ranges()` (requires the
#' **ebirdst** package and a downloaded species dataset).  All season polygons
#' are unioned into a single presence mask before cells are set to `NA`.  Set
#' `use_range = FALSE` to skip this step, or if **ebirdst** is not available.
#'
#' @param model_fit The list returned by [fit_species_model()].
#' @param polygon An `sf` object defining the study area (used for the initial
#'   polygon mask and the plot extent).
#' @param species Common name used for range lookup and the plot title.
#' @param grid_res_km Prediction grid resolution in km. Default `1`.
#' @param peak_time Optional decimal-hour override for observation start time
#'   used in predictions. Defaults to the median observation start time in the
#'   training data (a "typical checklist" reference).
#' @param use_range If `TRUE` (default), cells outside the eBird species range
#'   (from `ebirdst::load_ranges`) are set to `NA`.  Falls back to polygon-only
#'   masking with a warning if the range cannot be loaded.
#' @param range_resolution Resolution passed to `ebirdst::load_ranges`;
#'   `"27km"` (default) or `"9km"`.
#'
#' @return A list with:
#' \describe{
#'   \item{predictions}{A `terra::SpatRaster` with layers `abd` and `abd_se`.}
#'   \item{plot}{A `ggplot2` map of log-transformed relative abundance.}
#' }
#'
#' @seealso [fit_species_model()], [estimate_abundance()]
#' @export
predict_species_map <- function(model_fit,
                                polygon,
                                species,
                                sci_name         = NULL,
                                grid_res_km      = 1,
                                peak_time        = NULL,
                                use_range        = TRUE,
                                botw_path        = NULL,
                                range_resolution = "9km") {

  gam_model <- model_fit$model
  cov_stack <- model_fit$cov_stack
  cov_cols  <- model_fit$cov_cols

  # ── Step 1: Build and predict over grid ──────────────────────────────────
  message("\n── Step 1/2: Predicting abundance ───────────────────────────")
  pred_surface <- make_prediction_surface(polygon, grid_res_km)
  pred_surface <- extract_covariates(pred_surface, cov_stack)
  pred_surface <- tidyr::drop_na(pred_surface, dplyr::all_of(cov_cols))

  r_pred <- predict_abundance(
    gam_model, pred_surface, polygon, grid_res_km, peak_time
  )

  # ── Step 2: Mask to species range ────────────────────────────────────────
  if (use_range) {
    message("\n── Step 2/2: Masking to species range ───────────────────────")

    range_sf <- NULL

    # Prefer BOTW if a path and scientific name are supplied
    if (!is.null(botw_path) && !is.null(sci_name) && !is.na(sci_name)) {
      range_sf <- load_range_botw(sci_name, botw_path)
      if (!is.null(range_sf)) {
        message(sprintf("  BOTW: %d polygon(s) loaded for '%s'.",
                        nrow(range_sf), sci_name))
      } else {
        warning("BOTW: no range found for '", sci_name,
                "'; falling back to ebirdst.")
      }
    }

    # Fall back to ebirdst if BOTW not available or not found
    if (is.null(range_sf)) {
      range_sf <- load_range_ebirdst(species, range_resolution)
      if (!is.null(range_sf) && nrow(range_sf) > 0) {
        message(sprintf("  ebirdst: %d season polygon(s) loaded.", nrow(range_sf)))
      }
    }

    if (!is.null(range_sf) && nrow(range_sf) > 0) {
      # BirdLife polygons often have near-duplicate vertices that the S2
      # spherical engine rejects. Disable S2 so GEOS handles the repair.
      old_s2 <- sf::sf_use_s2()
      sf::sf_use_s2(FALSE)
      range_sf <- sf::st_make_valid(sf::st_union(
        sf::st_transform(range_sf, 4326)
      ))
      sf::sf_use_s2(old_s2)
      range_vect <- terra::vect(range_sf)
      r_pred <- terra::mask(r_pred, range_vect)
    } else {
      warning("Could not load range for '", species,
              "'; abundance map is unmasked.")
    }
  }

  message("\nDone.")

  list(
    predictions = r_pred,
    plot        = plot_abundance(r_pred, polygon, species)
  )
}


#' Estimate relative abundance within a polygon
#'
#' Convenience wrapper that runs the full pipeline in one call by combining
#' [fit_species_model()] and [predict_species_map()].  For more control — or to
#' reuse a fitted model across multiple prediction settings — call those two
#' functions directly.
#'
#' It is strongly recommended to call [prepare_covariates()] first and pass the
#' result via `cov_stack`.  This separates the one-time raster download (which
#' can exceed 1 GB) from the per-species modelling.
#'
#' @param polygon An `sf` object (POLYGON or MULTIPOLYGON) defining the
#'   study area. Any CRS is accepted; it is reprojected internally.
#' @param ebird_zip Path to the raw eBird EBD `.zip` *or* `.txt` file.
#' @param sampling_txt Path to the eBird sampling-events `.txt` file.
#' @param species Common name exactly as it appears in eBird
#'   (e.g. `"Superb Fairywren"`).
#' @param cov_stack A `terra::SpatRaster` returned by [prepare_covariates()].
#'   If `NULL` (default), [prepare_covariates()] is called automatically using
#'   `cache_dir`.
#' @param cache_dir Directory for the zero-filled eBird cache. Defaults to
#'   `"ebirdabund_cache"` in the working directory.
#' @param grid_res_km Prediction grid resolution in km. Default `1`.
#' @param hex_spacing_km Hex-cell diameter for spatiotemporal subsampling in
#'   km. Default `5`.
#' @param peak_time Optional decimal-hour override for observation start time
#'   used in predictions. Estimated from the model when `NULL`.
#' @param use_range If `TRUE` (default), cells outside the eBird species range
#'   are set to `NA`. See [predict_species_map()] for details.
#' @param range_resolution Resolution for `ebirdst::load_ranges`; `"27km"`
#'   (default) or `"9km"`.
#'
#' @return A named list with:
#' \describe{
#'   \item{model}{The fitted `mgcv::gam` object.}
#'   \item{predictions}{A `terra::SpatRaster` with layers `abd` and `abd_se`.}
#'   \item{data}{The subsampled data.frame used for model fitting.}
#'   \item{plot}{A `ggplot2` map of log-transformed relative abundance.}
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
#' # Step 2 — modelling (once per species)
#' model_fit <- fit_species_model(
#'   polygon      = nsw,
#'   ebird_zip    = "raw_data/ebd_AU-NSW.txt",
#'   sampling_txt = "raw_data/ebd_AU-NSW_sampling.txt",
#'   species      = "Superb Fairywren",
#'   cov_stack    = cov
#' )
#'
#' # Step 3 — prediction (masked to the species range)
#' result <- predict_species_map(
#'   model_fit,
#'   polygon = nsw,
#'   species = "Superb Fairywren"
#' )
#' result$plot
#' }
#'
#' @export
estimate_abundance <- function(polygon,
                               ebird_zip,
                               sampling_txt,
                               species,
                               sci_name         = NULL,
                               cov_stack        = NULL,
                               cache_dir        = "ebirdabund_cache",
                               grid_res_km      = 1,
                               hex_spacing_km   = 5,
                               peak_time        = NULL,
                               use_range        = TRUE,
                               botw_path        = NULL,
                               range_resolution = "27km") {

  model_fit <- fit_species_model(
    polygon        = polygon,
    ebird_zip      = ebird_zip,
    sampling_txt   = sampling_txt,
    species        = species,
    cov_stack      = cov_stack,
    cache_dir      = cache_dir,
    hex_spacing_km = hex_spacing_km
  )

  pred_out <- predict_species_map(
    model_fit        = model_fit,
    polygon          = polygon,
    species          = species,
    sci_name         = sci_name,
    grid_res_km      = grid_res_km,
    peak_time        = peak_time,
    use_range        = use_range,
    botw_path        = botw_path,
    range_resolution = range_resolution
  )

  list(
    model       = model_fit$model,
    predictions = pred_out$predictions,
    data        = model_fit$data,
    plot        = pred_out$plot
  )
}
