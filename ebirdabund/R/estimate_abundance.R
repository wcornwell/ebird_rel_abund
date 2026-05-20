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
#' @param max_count Upper cap on observation count; checklists above this are
#'   excluded (mega-flock filter). Default `200`.
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
                              hex_spacing_km = 5,
                              max_count      = 200L) {

  # ── Input validation ─────────────────────────────────────────────────────
  if (!inherits(polygon, "sf") && !inherits(polygon, "sfc")) {
    stop("`polygon` must be an sf or sfc object.")
  }
  missing_ebd <- ebird_zip[!file.exists(ebird_zip) &
                             !file.exists(tools::file_path_sans_ext(ebird_zip))]
  if (length(missing_ebd) > 0) {
    stop("eBird data file(s) not found: ", paste(missing_ebd, collapse = ", "))
  }
  missing_samp <- sampling_txt[!file.exists(sampling_txt)]
  if (length(missing_samp) > 0) {
    stop("Sampling events file(s) not found: ", paste(missing_samp, collapse = ", "))
  }
  if (!is.character(species) || nchar(trimws(species)) == 0) {
    stop("`species` must be a non-empty character string.")
  }

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Step 1: Load and filter eBird data ───────────────────────────────────
  message("\n── Step 1/4: Loading eBird data ─────────────────────────────")
  ebird_df <- load_ebird(polygon, ebird_zip, sampling_txt, species, cache_dir,
                         max_count = max_count)

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
  # Skipped if load_ebird() used the sampling master fast path (covariates
  # are already attached and drop_na was already applied).
  cov_cols <- grep(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height|nightlights|palsar_hv)", names(ebird_df), value = TRUE
  )
  if (length(cov_cols) == 0L) {
    message("\n── Step 3/4: Extracting covariates ──────────────────────────")
    ebird_df <- extract_covariates(ebird_df, cov_stack)
    cov_cols <- grep(
      "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height|nightlights|palsar_hv)", names(ebird_df), value = TRUE
    )
    n_before_drop <- nrow(ebird_df)
    ebird_df      <- tidyr::drop_na(ebird_df, dplyr::all_of(cov_cols))
    n_dropped     <- n_before_drop - nrow(ebird_df)
    if (n_dropped > 0) {
      message(sprintf(
        "  Dropped %d checklists with missing covariate values (%d remaining).",
        n_dropped, nrow(ebird_df)
      ))
    }
  } else {
    message("\n── Step 3/4: Covariates from sampling master (skipping extraction).")
  }

  # ── Step 4: Spatiotemporal subsampling + GAM ──────────────────────────────
  message("\n── Step 4/4: Spatiotemporal subsampling ─────────────────────")
  ebird_ss <- subsample_hex(ebird_df, spacing_km = hex_spacing_km)

  n_pos_ss <- sum(ebird_ss$observation_count > 0L)
  if (n_pos_ss < 50L) {
    cond <- structure(
      class = c("ebirdabund_excluded", "error", "condition"),
      list(
        message          = sprintf(
          "Insufficient detections: only %d positive checklists after subsampling (need >= 50).",
          n_pos_ss
        ),
        exclusion_reason = "insufficient_detections",
        n_checklists     = nrow(ebird_ss),
        n_positive       = n_pos_ss
      )
    )
    stop(cond)
  }

  message("\n── Fitting GAM ──────────────────────────────────────────────")
  gam_model <- fit_gam(ebird_ss)
  message(sprintf("Deviance explained: %.1f%%",
                  suppressWarnings(summary(gam_model)$dev.expl) * 100))

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
                                peak_doy         = NULL,
                                peak_time        = NULL,
                                use_range        = TRUE,
                                botw_path        = NULL,
                                range_resolution = "9km",
                                border           = NULL,
                                range_vect       = NULL) {

  gam_model <- model_fit$model
  cov_stack <- model_fit$cov_stack
  cov_cols  <- model_fit$cov_cols

  # ── Step 1: Build and predict over grid ──────────────────────────────────
  message("\n── Step 1/2: Predicting abundance ───────────────────────────")
  pred_surface <- make_prediction_surface(polygon, grid_res_km)
  pred_surface <- extract_covariates(pred_surface, cov_stack)
  pred_surface <- tidyr::drop_na(pred_surface, dplyr::all_of(cov_cols))

  abd_result <- predict_abundance(
    gam_model, pred_surface, polygon, grid_res_km, peak_doy, peak_time
  )
  r_pred    <- abd_result$predictions
  peak_doy  <- abd_result$peak_doy
  peak_time <- abd_result$peak_time

  # ── Step 2: Mask to species range ────────────────────────────────────────
  range_source <- "unmasked"

  if (use_range) {
    message("\n── Step 2/2: Masking to species range ───────────────────────")

    # If a pre-built range vector was supplied (e.g. reused across resolutions),
    # use it directly to avoid re-reading BOTW and re-running st_union.
    if (!is.null(range_vect)) {
      r_pred       <- terra::mask(r_pred, range_vect)
      src          <- attr(range_vect, "range_source")
      range_source <- if (!is.null(src)) src else "cached"
      message(sprintf("  %s: cached range applied.", range_source))
    } else {

      union_range <- function(range) {
        old_s2 <- sf::sf_use_s2()
        sf::sf_use_s2(FALSE)
        u <- sf::st_make_valid(sf::st_union(sf::st_transform(range, 4326)))
        sf::sf_use_s2(old_s2)
        terra::vect(u)
      }

      try_mask <- function(r, range_sf, source_label) {
        if (is.null(range_sf) || nrow(range_sf) == 0L) return(NULL)
        rv  <- union_range(range_sf)
        r_m <- terra::mask(r, rv)
        s   <- terra::global(r_m[["abd"]], "sum", na.rm = TRUE)[[1L]]
        if (is.na(s) || s <= 1e-5) {
          warning(source_label, " range for '", species,
                  "' does not overlap study region — trying next source.")
          return(NULL)
        }
        message(sprintf("  %s: %d polygon(s) applied.", source_label, nrow(range_sf)))
        list(raster = r_m, vect = rv)
      }

      mask_result  <- NULL

      # 1. Try ebirdst first
      range_ebirdst <- load_range_ebirdst(species, range_resolution)
      mask_result   <- try_mask(r_pred, range_ebirdst, "ebirdst")
      if (!is.null(mask_result)) range_source <- "ebirdst"

      # 2. Fall back to BOTW
      if (is.null(mask_result) && !is.null(botw_path) &&
          !is.null(sci_name) && !is.na(sci_name)) {
        range_botw  <- load_range_botw(sci_name, botw_path)
        mask_result <- try_mask(r_pred, range_botw, "BOTW")
        if (!is.null(mask_result)) range_source <- "BOTW"
      }

      if (!is.null(mask_result)) {
        r_pred     <- mask_result$raster
        range_vect <- mask_result$vect
        attr(range_vect, "range_source") <- range_source
      } else {
        warning("No usable range found for '", species, "' — map is unmasked.")
        range_vect <- NULL
      }
    }
  }

  message("\nDone.")

  list(
    predictions  = r_pred,
    plot         = plot_abundance(r_pred, polygon, species, border = border),
    range_source = range_source,
    range_vect   = range_vect,
    peak_doy     = peak_doy,
    peak_time    = peak_time
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
#' @param max_count Upper cap on observation count (mega-flock filter).
#'   Default `200`.
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
                               max_count        = 200L,
                               peak_doy         = NULL,
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
    hex_spacing_km = hex_spacing_km,
    max_count      = max_count
  )

  pred_out <- predict_species_map(
    model_fit        = model_fit,
    polygon          = polygon,
    species          = species,
    sci_name         = sci_name,
    grid_res_km      = grid_res_km,
    peak_doy         = peak_doy,
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
