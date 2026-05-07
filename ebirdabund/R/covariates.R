# Load JRC Global Surface Water occurrence tiles (Pekel et al. 2021) via VSICURL.
# Tiles are 10°×10°, named by their upper-left (top-left) corner.
# bb: c(xmin, ymin, xmax, ymax) in WGS84; ext: terra::ext; template: reference raster.
load_jrc_water <- function(bb, ext, template) {
  lon_lefts <- seq(floor(bb[1] / 10) * 10, floor(bb[3] / 10) * 10, by = 10)
  # lat_tops: northernmost latitude of each tile (ceiling handles S hemisphere correctly)
  lat_tops  <- seq(ceiling(bb[4] / 10) * 10, ceiling(bb[2] / 10) * 10, by = -10)
  base_url  <- "https://storage.googleapis.com/global-surface-water/downloads2021/occurrence/"
  tiles <- list()
  for (lon0 in lon_lefts) {
    for (lat0 in lat_tops) {
      lon_str <- sprintf("%d%s", abs(lon0), if (lon0 >= 0) "E" else "W")
      lat_str <- sprintf("%d%s", abs(lat0), if (lat0 >= 0) "N" else "S")
      fname   <- sprintf("occurrence_%s_%sv1_4_2021.tif", lon_str, lat_str)
      r <- tryCatch(terra::crop(terra::rast(paste0("/vsicurl/", base_url, fname)), ext),
                    error = function(e) NULL)
      if (!is.null(r)) tiles[[length(tiles) + 1L]] <- r
    }
  }
  if (length(tiles) == 0L) stop("No JRC GSW tiles found for extent")
  merged <- if (length(tiles) == 1L) tiles[[1L]] else do.call(terra::mosaic, tiles)
  terra::resample(merged, template)
}

# Load SoilGrids clay content (0-5 cm mean, 250 m) via VSICURL.
# SoilGrids uses Homolosine projection — we crop in native CRS before reprojecting.
# Returns percent clay (0–100).
load_clay <- function(bb, ext, template) {
  clay_vsi   <- geodata::soil_world_vsi(var = "clay", depth = 5, stat = "mean")
  native_crs <- terra::crs(clay_vsi)
  bb_pts <- terra::project(
    terra::vect(data.frame(x = c(bb[1], bb[3]), y = c(bb[2], bb[4])),
                geom = c("x", "y"), crs = "EPSG:4326"),
    native_crs
  )
  e  <- as.vector(terra::ext(bb_pts))  # xmin, xmax, ymin, ymax in native CRS
  de <- max(e[2] - e[1], e[4] - e[3]) * 0.05
  ext_native <- terra::ext(e[1] - de, e[2] + de, e[3] - de, e[4] + de)
  clay_crop  <- terra::crop(clay_vsi, ext_native)
  clay_wgs84 <- terra::project(clay_crop, "EPSG:4326")
  # SoilGrids clay is in g/kg (0–1000); divide by 10 → percent
  terra::resample(terra::crop(clay_wgs84, ext) / 10, template)
}

# Load OzTreeMap canopy height (Pucino et al. 2025, CSIRO) via the EASI WCS endpoint.
# Publicly accessible — no credentials required.
# The server limits each request to 32 tiles, so we split into chunk_deg×chunk_deg
# chunks and mosaic. Each chunk is snapped to the template pixel grid before
# mosaicking so terra::mosaic() sees consistent resolutions.
# variant: "median" (highest accuracy) or "best_pick" (vegetation-class-specific).
load_oztreemap <- function(bb, ext, template, cache_dir, variant = "median",
                           chunk_deg = 4) {
  cache_file <- file.path(cache_dir, sprintf(
    "oztreemap_%s_%.2f_%.2f_%.2f_%.2f.tif", variant, bb[1], bb[2], bb[3], bb[4]
  ))

  if (!file.exists(cache_file)) {
    # Clip to OzTreeMap's actual WCS coverage to avoid 500 errors on empty tiles
    oz_bb  <- c(112.92, -43.76, 154.10, -9.86)
    eff_bb <- c(max(bb[1], oz_bb[1]), max(bb[2], oz_bb[2]),
                min(bb[3], oz_bb[3]), min(bb[4], oz_bb[4]))
    if (eff_bb[1] >= eff_bb[3] || eff_bb[2] >= eff_bb[4]) {
      warning("Study extent does not overlap OzTreeMap coverage — tree_height will be NA")
      r <- template; terra::values(r) <- NA_real_; return(r)
    }

    # Snap to nearest 1/120° (standard 30-arc-second grid) so the WCS request
    # always gets an exact value; minor rounding in terra ext/res arithmetic
    # produces e.g. 0.008334 which the server rejects with HTTP 400.
    res      <- round(terra::res(template)[1] / (1 / 120)) * (1 / 120)
    wcs_base <- "https://ows.csiro.easi-eo.solutions/"

    lon_breaks <- seq(floor(eff_bb[1] / chunk_deg) * chunk_deg,
                      ceiling(eff_bb[3] / chunk_deg) * chunk_deg,
                      by = chunk_deg)
    lat_breaks <- seq(floor(eff_bb[2] / chunk_deg) * chunk_deg,
                      ceiling(eff_bb[4] / chunk_deg) * chunk_deg,
                      by = chunk_deg)

    n_chunks <- (length(lon_breaks) - 1L) * (length(lat_breaks) - 1L)
    message(sprintf("  Downloading OzTreeMap %s via WCS (%d chunks)...",
                    variant, n_chunks))

    # WCS server rejects requests wider than ~480 pixels per side, so we cannot
    # add a meaningful overlap buffer without exceeding the limit for interior
    # chunks.  Seam artefacts at 1 km / k=4 GAM scale are negligible.
    overlap <- 0

    tiles <- list()
    for (i in seq_len(length(lon_breaks) - 1L)) {
      for (j in seq_len(length(lat_breaks) - 1L)) {
        cb <- c(
          max(lon_breaks[i],     eff_bb[1]), max(lat_breaks[j],     eff_bb[2]),
          min(lon_breaks[i + 1], eff_bb[3]), min(lat_breaks[j + 1], eff_bb[4])
        )
        if (cb[1] >= cb[3] || cb[2] >= cb[4]) next

        # Buffered bbox for the WCS request; clipped to dataset bounds
        dl <- c(
          max(cb[1] - overlap, oz_bb[1]), max(cb[2] - overlap, oz_bb[2]),
          min(cb[3] + overlap, oz_bb[3]), min(cb[4] + overlap, oz_bb[4])
        )

        wcs_url <- paste0(
          wcs_base,
          "?service=WCS&version=1.0.0&request=GetCoverage",
          "&coverage=oztreemap_chm_", variant,
          "&format=GeoTIFF",
          sprintf("&bbox=%.4f,%.4f,%.4f,%.4f", dl[1], dl[2], dl[3], dl[4]),
          "&crs=EPSG:4326",
          sprintf("&resx=%.6f&resy=%.6f", res, res)
        )

        tmp <- tempfile(fileext = ".tif")
        ok  <- tryCatch({
          utils::download.file(wcs_url, tmp, mode = "wb", quiet = TRUE)
          file.exists(tmp) && file.size(tmp) > 1000L
        }, error = function(e) FALSE)

        if (ok) {
          r <- tryCatch(terra::rast(tmp), error = function(e) NULL)
          if (!is.null(r)) {
            # WCS returns 2-band GeoTIFF (best_pick + median); keep requested one
            if (terra::nlyr(r) > 1L) r <- r[[variant]]
            # Resample onto template-aligned grid, then crop to non-overlapping
            # core so adjacent tiles don't double-count the buffer strip.
            sub_tmpl <- terra::crop(template, terra::ext(r), snap = "out")
            r_res <- terra::resample(r, sub_tmpl, method = "bilinear")
            core_ext <- terra::ext(cb[1], cb[3], cb[2], cb[4])
            tiles[[length(tiles) + 1L]] <- terra::crop(r_res, core_ext)
          }
        }
        if (file.exists(tmp)) file.remove(tmp)
      }
    }

    if (length(tiles) == 0L) {
      warning("OzTreeMap WCS: no chunks succeeded — tree_height will be NA")
      r <- template
      terra::values(r) <- NA_real_
      return(r)
    }

    merged <- if (length(tiles) == 1L) {
      tiles[[1L]]
    } else {
      do.call(terra::mosaic, tiles)
    }
    # WCS GeoTIFF may include a second mask band — keep only band 1
    if (terra::nlyr(merged) > 1L) merged <- merged[[1L]]
    # Smooth Sentinel-2 tile-boundary artefacts before caching.
    # 3x3 focal mean at 1 km scale; negligible effect on k=4 GAM smooth.
    merged <- terra::clamp(
      terra::focal(merged, w = 3, fun = "mean",
                   na.policy = "omit", na.rm = TRUE),
      lower = 0, values = TRUE
    )
    terra::writeRaster(merged, cache_file, overwrite = TRUE)
    message(sprintf("  OzTreeMap cached: %s", basename(cache_file)))
  }

  raw <- terra::rast(cache_file)
  # Final resample to template can nudge near-zero cells slightly negative.
  resampled <- terra::resample(terra::crop(raw, ext), template)
  terra::clamp(resampled, lower = 0, values = TRUE)
}

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
#'   \code{temp_annual}, \code{pop_density} (log10 persons per km²,
#'   WorldPop 2020 via \code{geodata::population()}), and
#'   \code{water_occ} (JRC Global Surface Water occurrence 0–100,
#'   Pekel et al. 2016 updated 2021, streamed via VSICURL).
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

  # Cache key: rounded bbox + version tag (bump version when layer set changes)
  key        <- paste(round(bb, 2), collapse = "_")
  stack_path <- file.path(cache_dir, paste0("cov_stack_v4_", key, ".tif"))

  if (file.exists(stack_path)) {
    message("Loading cached covariate stack from ", stack_path)
    return(terra::rast(stack_path))
  }

  message("Building covariate stack (downloads cached after first run)...")
  # lc_water dropped in favour of water_occ (JRC GSW, 30 m, see below)
  lc_vars <- c("trees", "grassland", "shrubs", "cropland", "built")

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
  wc <- if (length(wc_tiles) == 1L) {
    wc_tiles[[1L]]
  } else {
    do.call(terra::mosaic, wc_tiles)
  }
  wc_crop <- terra::crop(wc, ext)
  # Extract annual precipitation (BIO12) and annual temperature (BIO1)
  # Use indices rather than names — tile resolution suffix varies (30s vs 2.5m)
  precip_annual <- terra::resample(wc_crop[[12]], lc_layers[[1]])
  temp_annual   <- terra::resample(wc_crop[[1]],  lc_layers[[1]])

  message("  population density (WorldPop 2020)")
  pop_raw     <- geodata::population(year = 2020, res = "5", path = cache_dir)
  pop_cropped <- terra::crop(pop_raw, ext)
  # log10-transform: persons/km² spans orders of magnitude; +1 avoids log(0)
  pop_density <- terra::resample(log10(pop_cropped + 1), lc_layers[[1]])

  message("  water occurrence (JRC Global Surface Water, Pekel et al. 2021)")
  water_occ <- load_jrc_water(bb, ext, lc_layers[[1]])

  message("  clay content (SoilGrids 250 m, 0-5 cm)")
  clay <- load_clay(bb, ext, lc_layers[[1]])

  message("  canopy height (OzTreeMap 30 m, Pucino et al. 2025)")
  tree_height <- load_oztreemap(bb, ext, lc_layers[[1]], cache_dir)

  stack <- terra::rast(c(
    lc_layers,
    list(elev, precip_annual, temp_annual, pop_density,
         water_occ, clay, tree_height)
  ))
  names(stack) <- c(
    paste0("lc_", lc_vars),
    "elevation", "precip_annual", "temp_annual", "pop_density", "water_occ",
    "clay", "tree_height"
  )

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

  # OzTreeMap returns NA where there is no tree canopy (open grassland,
  # cropland, urban, water). Treat as height = 0, not missing data.
  if ("tree_height" %in% names(vals)) {
    vals$tree_height[is.na(vals$tree_height)] <- 0
  }
  # SoilGrids clay returns NA at ocean/water. Treat as 0 (no soil).
  if ("clay" %in% names(vals)) {
    vals$clay[is.na(vals$clay)] <- 0
  }

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
