# Helper for the PALSAR HV experiment. Pull ALOS-2 PALSAR-2 yearly mosaic
# gamma0 (HV polarization) from Microsoft Planetary Computer's COG mirror,
# convert DN -> dB, aggregate to the master template.
#
# Access: STAC search -> per-asset SAS signing -> VSICURL on signed URL.
# No subscription key required; signed URLs expire after ~1 hour but that's
# fine since we cache the merged result to disk.
#
# If the experiment proves useful, this should be moved to
# ebirdabund/R/covariates.R and wired into download_covariates().

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

#' Load PALSAR HV gamma0 for a bounding box, aggregated to a template
#'
#' @param bb numeric(4): xmin, ymin, xmax, ymax in WGS84.
#' @param ext terra::ext object for the study region (used for final crop).
#' @param template SpatRaster to align resolution / extent to.
#' @param cache_dir directory for caching the merged dB raster.
#' @param year integer; PALSAR mosaic year (2015-2021 available on MPC).
#' @return SpatRaster of HV gamma0 (dB), aligned to template.
load_palsar_hv <- function(bb, ext, template, cache_dir, year = 2020L,
                           overview_level = 3L) {
  # PALSAR COGs are 4500x4500 at 25 m native, with overviews at:
  #   L0 native (4500x4500, 25 m), L1 (2250, 50 m), L2 (1125, 100 m),
  #   L3 (563, 200 m), L4 (282, 400 m).
  # We aggregate to a 1 km template, so reading at L3 (~200 m) is plenty fine
  # and downloads ~64x less data than native. Note: the COG-builtin average
  # resampling averages DN values, not power — slightly different from
  # averaging in linear power space, but well-behaved for relative comparisons.
  cache_file <- file.path(cache_dir, sprintf(
    "palsar_hv_%d_%.2f_%.2f_%.2f_%.2f.tif",
    year, bb[1], bb[2], bb[3], bb[4]
  ))

  if (!file.exists(cache_file)) {
    Sys.setenv(
      VSI_CACHE                    = "TRUE",
      VSI_CACHE_SIZE               = "536870912",
      GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
      GDAL_HTTP_MULTIPLEX          = "YES",
      GDAL_HTTP_VERSION            = "2"
    )

    stac_url <- sprintf(
      paste0("https://planetarycomputer.microsoft.com/api/stac/v1/search?",
             "collections=alos-palsar-mosaic&bbox=%f,%f,%f,%f",
             "&datetime=%d-01-01/%d-12-31&limit=500"),
      bb[1], bb[2], bb[3], bb[4], year, year
    )
    items <- fromJSON(stac_url, simplifyVector = FALSE)$features
    message(sprintf("  PALSAR STAC: %d items for year %d", length(items), year))
    if (length(items) == 0L) stop("No PALSAR items in bbox/year")

    tiles <- list()
    t0 <- Sys.time()
    for (i in seq_along(items)) {
      it   <- items[[i]]
      href <- it$assets$HV$href
      signed <- fromJSON(sprintf(
        "https://planetarycomputer.microsoft.com/api/sas/v1/sign?href=%s",
        URLencode(href, reserved = TRUE)
      ))$href
      r <- tryCatch(
        terra::rast(paste0("/vsicurl/", signed),
                    opts = sprintf("OVERVIEW_LEVEL=%d", overview_level)),
        error = function(e) { message("    open failed: ", e$message); NULL }
      )
      if (is.null(r)) next
      cr <- tryCatch(
        terra::crop(r, ext, snap = "out"),
        error = function(e) NULL
      )
      if (!is.null(cr) && prod(dim(cr)[1:2]) > 0L) {
        tiles[[length(tiles) + 1L]] <- cr
      }
      if (i %% 25L == 0L) {
        elapsed <- as.numeric(Sys.time() - t0, units = "secs")
        message(sprintf("    %d/%d items in %.0fs", i, length(items), elapsed))
      }
    }
    if (length(tiles) == 0L) stop("No PALSAR tiles loaded")

    merged <- if (length(tiles) == 1L) tiles[[1L]]
              else terra::mosaic(terra::sprc(tiles), fun = "mean")

    # JAXA calibration: gamma0 (dB) = 20 * log10(DN) - 83.0; DN == 0 is no-data.
    merged <- terra::ifel(merged == 0, NA, 20 * log10(merged) - 83.0)

    terra::writeRaster(
      merged, cache_file, overwrite = TRUE,
      gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER")
    )
    message(sprintf("  PALSAR HV cached: %s (%d tiles)",
                    basename(cache_file), length(tiles)))
  }

  raw      <- terra::rast(cache_file)
  raw_crop <- terra::crop(raw, ext)

  # Aggregate to template resolution. SAR backscatter must be averaged in
  # linear power, not dB — dB is a log scale, so arithmetic-mean of dB
  # values underweights bright pixels.
  ratio <- round(terra::res(template)[1] / terra::res(raw_crop)[1])
  if (ratio > 1L) {
    raw_lin <- 10^(raw_crop / 10)
    agg_lin <- terra::aggregate(raw_lin, fact = ratio, fun = "mean",
                                na.rm = TRUE)
    raw_crop <- 10 * log10(agg_lin)
  }
  terra::resample(raw_crop, template, method = "bilinear")
}
