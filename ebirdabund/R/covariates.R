#' Download and cache habitat covariates for a study region
#'
#' Downloads ESA WorldCover land-cover layers and SRTM elevation from the
#' internet (via \code{geodata}), crops them to the polygon extent, and
#' saves the combined raster stack to \code{cache_dir} as a GeoTIFF.
#' On subsequent calls the cached file is loaded directly — no
#' re-downloading needed.
#'
#' Run this \emph{once per study region} before calling
#' \code{\link{estimate_abundance}}.
#'
#' @param polygon An \code{sf} object defining the study area.
#' @param cache_dir Directory for cached rasters. Created if absent.
#'   Defaults to \code{"ebirdabund_cache"}.
#' @param buffer_deg Degrees of buffer added around the polygon extent before
#'   downloading. Default \code{0.5}.
#'
#' @return A \code{terra::SpatRaster} with layers \code{lc_trees},
#'   \code{lc_grassland}, \code{lc_shrubs}, \code{lc_cropland},
#'   \code{lc_built}, \code{lc_water}, \code{elevation}, \code{precip_annual},
#'   \code{temp_annual}.
#'
#' @export
prepare_covariates <- function(polygon,
                               cache_dir  = "ebirdabund_cache",
                               buffer_deg = 0.5) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  bb  <- as.numeric(sf::st_bbox(sf::st_transform(polygon, 4326)))
  ext <- terra::ext(
    bb[1] - buffer_deg, bb[3] + buffer_deg,
    bb[2] - buffer_deg, bb[4] + buffer_deg
  )

  # Cache key: rounded bbox so tiny reprojection differences don't miss cache
  key        <- paste(round(bb, 2), collapse = "_")
  stack_path <- file.path(cache_dir, paste0("cov_stack_", key, ".tif"))

  if (file.exists(stack_path)) {
    message("Loading cached covariate stack from ", stack_path)
    return(terra::rast(stack_path))
  }

  message("Building covariate stack (downloads cached after first run)...")
  lc_vars <- c("trees", "grassland", "shrubs", "cropland", "built", "water")

  lc_layers <- lapply(lc_vars, function(v) {
    message("  landcover: ", v)
    terra::crop(geodata::landcover(v, path = cache_dir), ext)
  })

  message("  elevation")
  elev_raw <- geodata::elevation_30s(country = "AUS", path = cache_dir)
  elev     <- terra::resample(terra::crop(elev_raw, ext), lc_layers[[1]])

  message("  climate: precipitation and temperature")
  # WorldClim tiles are 30-degree wide; collect all tiles that intersect the
  # extent (a single-centre lookup silently drops any tile boundary crossing).
  tile_lon_mins <- seq(floor(bb[1] / 30) * 30, floor(bb[3] / 30) * 30, by = 30)
  tile_lat_mins <- seq(floor(bb[2] / 30) * 30, floor(bb[4] / 30) * 30, by = 30)
  wc_tiles <- list()
  for (tlon in tile_lon_mins) {
    for (tlat in tile_lat_mins) {
      wc_tiles[[length(wc_tiles) + 1L]] <- geodata::worldclim_tile(
        lon  = tlon + 15, lat  = tlat + 15,
        var  = "bio", res = 2.5, path = cache_dir
      )
    }
  }
  wc      <- if (length(wc_tiles) == 1L) wc_tiles[[1L]] else
               do.call(terra::mosaic, wc_tiles)
  wc_crop <- terra::crop(wc, ext)
  # Extract annual precipitation (BIO12) and annual temperature (BIO1)
  # Use indices rather than names — tile resolution suffix varies (30s vs 2.5m)
  precip_annual <- terra::resample(wc_crop[[12]], lc_layers[[1]])
  temp_annual   <- terra::resample(wc_crop[[1]],  lc_layers[[1]])

  stack        <- terra::rast(c(lc_layers, list(elev, precip_annual, temp_annual)))
  names(stack) <- c(paste0("lc_", lc_vars), "elevation", "precip_annual", "temp_annual")

  terra::writeRaster(stack, stack_path, overwrite = TRUE)
  message("Covariate stack saved to ", stack_path)

  stack
}

# Extract covariate values at lon/lat locations in df.
# df must have columns 'longitude' and 'latitude' (WGS84 decimal degrees).
# Returns df with covariate columns appended.
extract_covariates <- function(df, cov_stack) {
  pts <- terra::vect(
    data.frame(x = df$longitude, y = df$latitude),
    geom = c("x", "y"),
    crs  = "EPSG:4326"
  )

  vals <- terra::extract(cov_stack, pts, method = "bilinear")
  vals <- vals[, -1, drop = FALSE]  # drop ID column

  dplyr::bind_cols(df, vals)
}

# Build a regular prediction grid inside polygon.
# Returns a data.frame with longitude and latitude columns only
# (covariates added separately).
make_prediction_surface <- function(polygon, grid_res_km) {
  polygon_wgs84 <- sf::st_transform(polygon, 4326)
  res_deg <- grid_res_km / 111.32

  # Cover the full bounding box — no polygon filter here.
  # Masking is applied later in predict_abundance(). Predicting at the full
  # bbox avoids NA gaps in the output raster where the polygon is concave.
  grid_pts <- sf::st_make_grid(
    polygon_wgs84,
    cellsize = res_deg,
    what     = "centers"
  ) |> sf::st_as_sf()

  coords <- sf::st_coordinates(grid_pts)
  data.frame(longitude = coords[, 1], latitude = coords[, 2])
}
