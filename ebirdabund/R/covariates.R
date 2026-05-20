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

# WGS84 lon/lat → slippy-map tile (x, y) at zoom z (Bing/Google/OSM scheme).
.lonlat_to_tile <- function(lon, lat, z) {
  n <- 2^z
  x <- floor((lon + 180) / 360 * n)
  lat_rad <- lat * pi / 180
  y <- floor((1 - log(tan(lat_rad) + 1 / cos(lat_rad)) / pi) / 2 * n)
  list(x = as.integer(x), y = as.integer(y))
}

# Tile (x, y) at zoom z → Bing quadkey base-4 string (length z).
.xy_to_quadkey <- function(x, y, z) {
  q <- character(z)
  for (i in seq.int(z, 1L)) {
    digit <- 0L
    mask  <- bitwShiftL(1L, i - 1L)
    if (bitwAnd(x, mask) != 0L) digit <- digit + 1L
    if (bitwAnd(y, mask) != 0L) digit <- digit + 2L
    q[z - i + 1L] <- as.character(digit)
  }
  paste0(q, collapse = "")
}

# Load Meta/WRI Canopy Height Maps v2 (DINOv3, 2024 imagery, ~1 m native) from
# the public AWS Open Data bucket (CC-BY 4.0, no auth). Tiles are 10-digit Bing
# quadkeys at zoom 10 in EPSG:3857, stored as COGs with built-in overviews
# (16384²/8192²/4096²/2048²/1024²/512² = 2.4/4.8/9.5/19/38/76 m).
#
# We force OVERVIEW_LEVEL=4 (1024×1024 ≈ 38 m effective) — the overview level
# closest to the ~30 m output we actually want. Reading at this level lets us
# use Meta's pre-built overview aggregation rather than averaging ourselves,
# and is ~150× faster than letting GDAL auto-pick. Each tile is projected onto
# a sub-template aligned to a master ~38 m WGS84 grid, then mosaicked.
load_meta_chmv2 <- function(bb, ext, template, cache_dir,
                            zoom = 10L, overview_level = 4L) {
  cache_file <- file.path(cache_dir, sprintf(
    "meta_chmv2_%.2f_%.2f_%.2f_%.2f.tif", bb[1], bb[2], bb[3], bb[4]
  ))

  if (!file.exists(cache_file)) {
    Sys.setenv(
      VSI_CACHE                          = "TRUE",
      VSI_CACHE_SIZE                     = "536870912",
      GDAL_DISABLE_READDIR_ON_OPEN       = "EMPTY_DIR",
      CPL_VSIL_CURL_ALLOWED_EXTENSIONS   = ".tif",
      GDAL_HTTP_MULTIPLEX                = "YES",
      GDAL_HTTP_VERSION                  = "2",
      GDAL_HTTP_MERGE_CONSECUTIVE_RANGES = "YES",
      GDAL_NUM_THREADS                   = "ALL_CPUS"
    )

    base_url <- paste0(
      "/vsicurl/https://dataforgood-fb-data.s3.amazonaws.com/",
      "forests/v2/global/dinov3_global_chm_v2_ml3/chm/"
    )

    # Quadkey set covering bbox (clamped to Web Mercator's valid lat range).
    lat_lo <- max(bb[2], -85.0511)
    lat_hi <- min(bb[4],  85.0511)
    nw <- .lonlat_to_tile(bb[1], lat_hi, zoom)
    se <- .lonlat_to_tile(bb[3], lat_lo, zoom)
    keys <- character(0)
    for (y in seq.int(nw$y, se$y)) {
      for (x in seq.int(nw$x, se$x)) {
        keys <- c(keys, .xy_to_quadkey(x, y, zoom))
      }
    }

    # Master ~38 m WGS84 grid covering ext (matches OVERVIEW_LEVEL=4 native
    # resolution at mid-latitudes). Sub-templates snap to this grid so adjacent
    # tiles mosaic cleanly without drift.
    chm_tmpl <- terra::rast(ext, resolution = 38 / 111320, crs = "EPSG:4326")

    message(sprintf(
      "  Streaming Meta CHMv2: %d candidate quadkeys (zoom %d, OL=%d)...",
      length(keys), zoom, overview_level
    ))
    t_start  <- Sys.time()
    tiles    <- list()
    n_loaded <- 0L
    n_missing <- 0L

    for (k in keys) {
      url <- paste0(base_url, k, ".tif")
      r <- tryCatch(
        suppressWarnings(suppressMessages(
          terra::rast(url, opts = sprintf("OVERVIEW_LEVEL=%d", overview_level))
        )),
        error = function(e) NULL
      )
      if (is.null(r)) {
        n_missing <- n_missing + 1L
        next
      }

      ext_v   <- terra::vect(terra::ext(r), crs = terra::crs(r))
      ext_wgs <- terra::ext(terra::project(ext_v, "EPSG:4326"))
      sub_tmpl <- tryCatch(
        terra::crop(chm_tmpl, ext_wgs, snap = "out"),
        error = function(e) NULL
      )
      if (is.null(sub_tmpl) || prod(dim(sub_tmpl)[1:2]) == 0L) {
        n_missing <- n_missing + 1L
        next
      }

      r_proj <- tryCatch(
        suppressWarnings(terra::project(r, sub_tmpl, method = "bilinear",
                                        threads = TRUE)),
        error = function(e) NULL
      )
      if (!is.null(r_proj)) {
        tiles[[length(tiles) + 1L]] <- r_proj
        n_loaded <- n_loaded + 1L
      } else {
        n_missing <- n_missing + 1L
      }

      if ((n_loaded + n_missing) %% 100L == 0L) {
        elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
        message(sprintf(
          "    %d loaded / %d missing of %d  —  %.0f sec elapsed",
          n_loaded, n_missing, length(keys), elapsed
        ))
      }
    }

    if (length(tiles) == 0L) {
      warning("Meta CHMv2: no tiles loaded — tree_height will be NA")
      r <- template
      terra::values(r) <- NA_real_
      return(r)
    }

    merged <- if (length(tiles) == 1L) {
      tiles[[1L]]
    } else {
      terra::mosaic(terra::sprc(tiles), fun = "mean")
    }
    # CHMv2 produces a small number of unphysical outliers — ~0.001% of pixels
    # exceed plausible canopy height. Most are model misclassifications of tall
    # narrow non-vegetation: wind turbines on ridge crests (NSW turbines are
    # 100–200 m to rotor tip), transmission towers, silos. Verified that they
    # are NOT cliff-edge artefacts (outlier slope distribution matches
    # background). NSW old-growth eucalypts top out around 70–75 m, so anything
    # above 60 m is set to NA rather than clamped — clamping would still
    # contaminate the cell mean with a non-vegetation pixel; NA + area-mean
    # downsample drops the bad pixel and averages over real canopy nearby.
    merged <- terra::clamp(merged, lower = 0, values = TRUE)
    merged <- terra::ifel(merged > 60, NA, merged)
    terra::writeRaster(
      merged, cache_file, overwrite = TRUE,
      gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER")
    )
    message(sprintf("  Meta CHMv2 cached: %s (%d/%d tiles loaded)",
                    basename(cache_file), n_loaded, length(keys)))
  }

  raw      <- terra::rast(cache_file)
  raw_crop <- terra::crop(raw, ext)
  # 90th-percentile area aggregation rather than mean. Captures emergent
  # canopy structure (tall trees in each cell — hollow-nesting parrot habitat,
  # raptor perches, etc.) without being dominated by isolated artefact
  # pixels: residual wind turbines / towers in the 40–60 m range sit in the
  # top ~0.1% of a 26×26 cell window and so are excluded by the p90 cut.
  # mean would smooth them in; max would let one bad pixel dominate.
  ratio <- round(terra::res(template)[1] / terra::res(raw_crop)[1])
  agg <- terra::aggregate(
    raw_crop, fact = ratio,
    fun = function(x, ...) stats::quantile(x, 0.9, na.rm = TRUE)
  )
  terra::resample(agg, template, method = "bilinear")
}

# Load Falchi/Cinzano World Atlas of Artificial Night Sky Brightness 2015
# from a local GeoTIFF. Global WGS84 at ~30 arc-seconds (~1 km), Float32
# radiance proxy (mcd/m²) with a sentinel NoData value over oceans.
# Crops in native CRS then resamples to template.
load_nightlights <- function(path, ext, template) {
  if (!file.exists(path)) {
    stop("Nightlights raster not found at: ", path)
  }
  r <- terra::rast(path)
  terra::resample(terra::crop(r, ext), template)
}

# Load ALOS-2 PALSAR-2 yearly mosaic gamma0 (HV polarization, L-band SAR)
# from Microsoft Planetary Computer's COG mirror.
#
# Access: STAC search → per-asset SAS signing → VSICURL on signed URL.
# No subscription key required; signed URLs expire after ~1 hour but the
# merged dB raster is cached to disk so subsequent calls skip the network.
#
# PALSAR COGs are 4500x4500 at 25 m native, with built-in overviews at:
#   L0 native (4500x4500, 25 m), L1 (2250, 50 m), L2 (1125, 100 m),
#   L3 (563, 200 m), L4 (282, 400 m).
# We aggregate to the ~1 km template, so reading at L3 (~200 m) is plenty
# fine and downloads ~64x less data than native.
#
# JAXA calibration: gamma0 (dB) = 20 * log10(DN) - 83.0; DN == 0 is no-data.
# SAR backscatter must be averaged in linear power, not dB, when aggregating
# (dB is a log scale; arithmetic-mean of dB underweights bright pixels).
load_palsar_hv <- function(bb, ext, template, cache_dir, year = 2020L,
                           overview_level = 3L) {
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
    items <- jsonlite::fromJSON(stac_url, simplifyVector = FALSE)$features
    message(sprintf("  PALSAR STAC: %d items for year %d",
                    length(items), year))
    if (length(items) == 0L) stop("No PALSAR items in bbox/year")

    tiles <- list()
    t0 <- Sys.time()
    for (i in seq_along(items)) {
      it     <- items[[i]]
      href   <- it$assets$HV$href
      signed <- jsonlite::fromJSON(sprintf(
        "https://planetarycomputer.microsoft.com/api/sas/v1/sign?href=%s",
        utils::URLencode(href, reserved = TRUE)
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
        message(sprintf("    %d/%d items in %.0fs",
                        i, length(items), elapsed))
      }
    }
    if (length(tiles) == 0L) stop("No PALSAR tiles loaded")

    merged <- if (length(tiles) == 1L) tiles[[1L]]
              else terra::mosaic(terra::sprc(tiles), fun = "mean")

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

  # Aggregate to template resolution in linear power, then back to dB.
  ratio <- round(terra::res(template)[1] / terra::res(raw_crop)[1])
  if (ratio > 1L) {
    raw_lin <- 10^(raw_crop / 10)
    agg_lin <- terra::aggregate(raw_lin, fact = ratio, fun = "mean",
                                na.rm = TRUE)
    raw_crop <- 10 * log10(agg_lin)
  }
  terra::resample(raw_crop, template, method = "bilinear")
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
#'   Pekel et al. 2016 updated 2021, streamed via VSICURL),
#'   \code{clay} (SoilGrids 0–5 cm clay percent, 250 m), and
#'   \code{tree_height} (Meta/WRI Canopy Height Maps v2 DINOv3, 2024 imagery,
#'   ~38 m effective via OVERVIEW_LEVEL=4 of the source COG pyramid),
#'   \code{nightlights} (Falchi/Cinzano World Atlas of Artificial Night Sky
#'   Brightness 2015, ~1 km, read from a local GeoTIFF), and
#'   \code{palsar_hv} (ALOS-2 PALSAR-2 yearly mosaic gamma0 HV polarization
#'   in dB, ~200 m via OVERVIEW_LEVEL=3 of the MPC COG mirror, aggregated
#'   to template in linear power then converted back to dB).
#'
#' @param nightlights_path Path to the World Atlas 2015 nightlights GeoTIFF.
#'   Defaults to \code{file.path(cache_dir, "nightlights/World_Atlas_2015.tif")}.
#' @param palsar_year Year of the ALOS-2 PALSAR-2 yearly mosaic (Microsoft
#'   Planetary Computer offers 2015–2021). Default \code{2020L}.
#'
#' @export
prepare_covariates <- function(polygon,
                               cache_dir        = "ebirdabund_cache",
                               buffer_deg       = 0.5,
                               nightlights_path = NULL,
                               palsar_year      = 2020L) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  bb  <- as.numeric(sf::st_bbox(sf::st_transform(polygon, 4326)))
  ext <- terra::ext(
    bb[1] - buffer_deg, bb[3] + buffer_deg,
    bb[2] - buffer_deg, bb[4] + buffer_deg
  )

  # Cache key: rounded bbox + version tag (bump version when layer set changes)
  # v5: tree_height switched from OzTreeMap (CSIRO WCS) to Meta/WRI CHMv2.
  # v6: added nightlights (Falchi/Cinzano World Atlas 2015).
  # v7: added palsar_hv (ALOS-2 PALSAR-2 yearly mosaic, MPC COG mirror).
  key        <- paste(round(bb, 2), collapse = "_")
  stack_path <- file.path(cache_dir, paste0("cov_stack_v7_", key, ".tif"))

  if (is.null(nightlights_path)) {
    nightlights_path <- file.path(cache_dir, "nightlights", "World_Atlas_2015.tif")
  }

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

  message("  canopy height (Meta/WRI CHMv2 DINOv3, ~38 m via OL=4)")
  tree_height <- load_meta_chmv2(bb, ext, lc_layers[[1]], cache_dir)

  message("  nightlights (Falchi/Cinzano World Atlas 2015)")
  nightlights <- load_nightlights(nightlights_path, ext, lc_layers[[1]])

  message(sprintf("  PALSAR HV gamma0 (ALOS-2 yearly mosaic, %d)", palsar_year))
  palsar_hv <- load_palsar_hv(bb, ext, lc_layers[[1]], cache_dir,
                              year = palsar_year)

  stack <- terra::rast(c(
    lc_layers,
    list(elev, precip_annual, temp_annual, pop_density,
         water_occ, clay, tree_height, nightlights, palsar_hv)
  ))
  names(stack) <- c(
    paste0("lc_", lc_vars),
    "elevation", "precip_annual", "temp_annual", "pop_density", "water_occ",
    "clay", "tree_height", "nightlights", "palsar_hv"
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

  # Meta CHMv2 returns NA where there is no tree canopy or no valid scene
  # (open grassland, cropland, urban, water). Treat as height = 0.
  if ("tree_height" %in% names(vals)) {
    vals$tree_height[is.na(vals$tree_height)] <- 0
  }
  # SoilGrids clay returns NA at ocean/water. Fill with median of non-NA
  # values so water cells don't pull predictions toward 0.
  if ("clay" %in% names(vals)) {
    clay_med <- median(vals$clay, na.rm = TRUE)
    vals$clay[is.na(vals$clay)] <- clay_med
  }
  # World Atlas NoData (oceans, polar gap) → 0 brightness.
  if ("nightlights" %in% names(vals)) {
    vals$nightlights[is.na(vals$nightlights)] <- 0
  }
  # PALSAR HV returns NA at ocean/coastal pixels (DN==0 masked at calibration)
  # and at mosaic tile edges. Median-impute so checklists near the coast
  # don't pull predictions to an extreme value; water_occ already encodes
  # the wet/dry signal separately.
  if ("palsar_hv" %in% names(vals)) {
    palsar_med <- stats::median(vals$palsar_hv, na.rm = TRUE)
    vals$palsar_hv[is.na(vals$palsar_hv)] <- palsar_med
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
