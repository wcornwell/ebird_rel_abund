#!/usr/bin/env Rscript
# Pull a 1 km x 1 km window at native ~1 m CHMv2 resolution around a point
# and save it as a PNG. Reuses the quadkey helpers from ebirdabund/R/covariates.R.

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(tidyterra)
})

PT_LON <- 141.704623
PT_LAT <- -31.087603
HALF_M <- 500          # half-window in metres (=> 1 km x 1 km)
ZOOM   <- 10L
OUT_TIF <- "analysis/chmv2_native_point.tif"
OUT_PNG <- "analysis/chmv2_native_point.png"

# Quadkey helpers (copied from ebirdabund/R/covariates.R — they're not exported).
lonlat_to_tile <- function(lon, lat, z) {
  n <- 2^z
  x <- floor((lon + 180) / 360 * n)
  lat_rad <- lat * pi / 180
  y <- floor((1 - log(tan(lat_rad) + 1 / cos(lat_rad)) / pi) / 2 * n)
  list(x = as.integer(x), y = as.integer(y))
}
xy_to_quadkey <- function(x, y, z) {
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

# GDAL env (same as load_meta_chmv2).
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

# Build window bbox in WGS84 around the point.
dlat <- HALF_M / 111320
dlon <- HALF_M / (111320 * cos(PT_LAT * pi / 180))
bb <- c(PT_LON - dlon, PT_LAT - dlat, PT_LON + dlon, PT_LAT + dlat)

# Cover the window with z=10 quadkeys (likely just one, but handle edges).
nw <- lonlat_to_tile(bb[1], bb[4], ZOOM)
se <- lonlat_to_tile(bb[3], bb[2], ZOOM)
keys <- character(0)
for (y in seq.int(nw$y, se$y)) {
  for (x in seq.int(nw$x, se$x)) {
    keys <- c(keys, xy_to_quadkey(x, y, ZOOM))
  }
}
message(sprintf("Quadkeys covering window: %s", paste(keys, collapse = ", ")))

base_url <- paste0(
  "/vsicurl/https://dataforgood-fb-data.s3.amazonaws.com/",
  "forests/v2/global/dinov3_global_chm_v2_ml3/chm/"
)

# Open at native resolution (no OVERVIEW_LEVEL). Each tile is ~16384x16384 px
# in EPSG:3857 covering one z=10 slippy tile (~38 km at this latitude),
# so the native pixel size is ~2.3 m on the ground at lat -31.
# We crop in EPSG:3857 first (cheap, no resample), then project to WGS84.
pt_3857 <- project(vect(data.frame(x = PT_LON, y = PT_LAT),
                        geom = c("x", "y"), crs = "EPSG:4326"),
                   "EPSG:3857")
xy <- crds(pt_3857)
# In EPSG:3857 1 metre at the equator = 1 unit, but at lat phi the metric
# scale is cos(phi). So to span +/- HALF_M ground metres we need
# +/- HALF_M / cos(phi) units.
scale_3857 <- 1 / cos(PT_LAT * pi / 180)
half_units <- HALF_M * scale_3857
ext_3857 <- ext(xy[1] - half_units, xy[1] + half_units,
                xy[2] - half_units, xy[2] + half_units)

tiles <- list()
for (k in keys) {
  url <- paste0(base_url, k, ".tif")
  message(sprintf("Opening %s", url))
  r <- tryCatch(rast(url), error = function(e) { message("  ", e$message); NULL })
  if (is.null(r)) next
  message(sprintf("  Native res: %.3f x %.3f units (EPSG:3857)",
                  res(r)[1], res(r)[2]))
  cr <- tryCatch(crop(r, ext_3857), error = function(e) NULL)
  if (!is.null(cr) && prod(dim(cr)[1:2]) > 0L) tiles[[length(tiles) + 1L]] <- cr
}

if (length(tiles) == 0L) stop("No CHMv2 tiles loaded for this window")

merged <- if (length(tiles) == 1L) tiles[[1L]] else mosaic(sprc(tiles), fun = "mean")
message(sprintf("Merged raster: %d x %d px, res %.3f m (in 3857 units)",
                nrow(merged), ncol(merged), res(merged)[1]))

# Project to WGS84 so the saved tif and the plot are easy to interpret.
chm <- project(merged, "EPSG:4326", method = "bilinear", threads = TRUE)

writeRaster(chm, OUT_TIF, overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
message(sprintf("Wrote %s", OUT_TIF))

# Plot — keep raw native pixels; no >60 m mask so you can see anything weird.
df_summary <- as.data.frame(values(chm))
names(df_summary) <- "h"
df_summary <- df_summary[!is.na(df_summary$h), , drop = FALSE]
message(sprintf("Height: min %.2f, median %.2f, p90 %.2f, p99 %.2f, max %.2f m  (n=%d non-NA)",
                min(df_summary$h), median(df_summary$h),
                quantile(df_summary$h, 0.9), quantile(df_summary$h, 0.99),
                max(df_summary$h), nrow(df_summary)))

p <- ggplot() +
  geom_spatraster(data = chm) +
  scale_fill_viridis_c(name = "Canopy height (m)", na.value = "transparent",
                       option = "viridis") +
  geom_point(aes(x = PT_LON, y = PT_LAT), colour = "red", size = 2, shape = 4) +
  coord_sf(crs = 4326) +
  labs(
    title = sprintf("Meta CHMv2 (DINOv3) native resolution @ (%.4f, %.4f)",
                    PT_LAT, PT_LON),
    subtitle = sprintf("~1 km x 1 km window, %d x %d pixels",
                       nrow(chm), ncol(chm)),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave(OUT_PNG, p, width = 7, height = 7, dpi = 150)
message(sprintf("Wrote %s", OUT_PNG))
