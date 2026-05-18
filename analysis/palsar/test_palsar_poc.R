#!/usr/bin/env Rscript
# Proof-of-concept: pull ALOS-2 PALSAR-2 HV around our test point, convert DN
# to gamma0 (dB), and plot alongside the CHMv2 layer for the same window.
# Goal: see whether HV picks up the 0.5-2 m woody veg CHMv2 misses.

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(tidyterra)
  library(jsonlite)
  library(patchwork)
})

PT_LON  <- 141.704623
PT_LAT  <- -31.087603
HALF_M  <- 2500            # 5 km x 5 km window — bigger than the CHMv2 test
                           # because PALSAR is 25 m native (vs 1 m for CHMv2),
                           # so a 1 km window would only be 40x40 px.
YEAR    <- 2020L           # MPC coverage: 2015-2021
OUT_TIF <- "analysis/palsar_hv_poc.tif"
OUT_PNG <- "analysis/palsar_hv_poc.png"

# ---- 1. Build bbox in WGS84 around the point ------------------------------
dlat <- HALF_M / 111320
dlon <- HALF_M / (111320 * cos(PT_LAT * pi / 180))
bb <- c(PT_LON - dlon, PT_LAT - dlat, PT_LON + dlon, PT_LAT + dlat)
ext_wgs <- ext(bb[1], bb[3], bb[2], bb[4])

# ---- 2. STAC search for HV assets covering the bbox -----------------------
stac_url <- sprintf(
  paste0("https://planetarycomputer.microsoft.com/api/stac/v1/search?",
         "collections=alos-palsar-mosaic&bbox=%f,%f,%f,%f",
         "&datetime=%d-01-01/%d-12-31&limit=200"),
  bb[1], bb[2], bb[3], bb[4], YEAR, YEAR
)
items <- fromJSON(stac_url, simplifyVector = FALSE)$features
message(sprintf("STAC returned %d items for year %d", length(items), YEAR))
if (length(items) == 0L) stop("No PALSAR items in bbox/year")

# ---- 3. Sign each HV asset and load via VSICURL ---------------------------
sign_url <- function(href) {
  resp <- fromJSON(sprintf(
    "https://planetarycomputer.microsoft.com/api/sas/v1/sign?href=%s",
    URLencode(href, reserved = TRUE)
  ))
  resp$href
}

Sys.setenv(
  VSI_CACHE                    = "TRUE",
  VSI_CACHE_SIZE               = "268435456",
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  GDAL_HTTP_MULTIPLEX          = "YES",
  GDAL_HTTP_VERSION            = "2"
)

tiles <- list()
for (it in items) {
  href <- it$assets$HV$href
  message("  ", it$id)
  signed <- sign_url(href)
  r <- rast(paste0("/vsicurl/", signed))
  cr <- tryCatch(crop(r, ext_wgs, snap = "out"), error = function(e) NULL)
  if (!is.null(cr) && prod(dim(cr)[1:2]) > 0L) tiles[[length(tiles) + 1L]] <- cr
}
if (length(tiles) == 0L) stop("No PALSAR tiles loaded")

merged <- if (length(tiles) == 1L) tiles[[1L]] else mosaic(sprc(tiles), fun = "mean")
message(sprintf("Merged PALSAR HV: %d x %d px, res %.5f deg",
                nrow(merged), ncol(merged), res(merged)[1]))

# ---- 4. DN -> gamma0 dB ---------------------------------------------------
# JAXA calibration: gamma0 (dB) = 20 * log10(DN) - 83.0; DN == 0 is no-data.
merged <- ifel(merged == 0, NA, 20 * log10(merged) - 83.0)

# Crop to exact window
hv_db <- crop(merged, ext_wgs)
writeRaster(hv_db, OUT_TIF, overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
message(sprintf("Wrote %s", OUT_TIF))

# ---- 5. Stats -------------------------------------------------------------
v <- values(hv_db); v <- v[!is.na(v)]
message(sprintf(
  "HV gamma0: min %.1f, p10 %.1f, median %.1f, p90 %.1f, max %.1f dB (n=%d)",
  min(v), quantile(v, 0.1), median(v),
  quantile(v, 0.9), max(v), length(v)
))

# ---- 6. Side-by-side plot with CHMv2 (matched extent) ---------------------
# Pull native CHMv2 directly (the load_meta_chmv2 wrapper has a template-
# alignment edge case for this small window — easier to just inline the parts
# we need, same approach as analysis/test_chmv2_native_point.R).
lonlat_to_tile <- function(lon, lat, z) {
  n <- 2^z
  list(x = as.integer(floor((lon + 180) / 360 * n)),
       y = as.integer(floor((1 - log(tan(lat * pi / 180) +
                                     1 / cos(lat * pi / 180)) / pi) / 2 * n)))
}
xy_to_quadkey <- function(x, y, z) {
  q <- character(z)
  for (i in seq.int(z, 1L)) {
    digit <- 0L; mask <- bitwShiftL(1L, i - 1L)
    if (bitwAnd(x, mask) != 0L) digit <- digit + 1L
    if (bitwAnd(y, mask) != 0L) digit <- digit + 2L
    q[z - i + 1L] <- as.character(digit)
  }
  paste0(q, collapse = "")
}

chm <- tryCatch({
  nw <- lonlat_to_tile(bb[1], bb[4], 10L)
  se <- lonlat_to_tile(bb[3], bb[2], 10L)
  keys <- character()
  for (y in seq.int(nw$y, se$y))
    for (x in seq.int(nw$x, se$x))
      keys <- c(keys, xy_to_quadkey(x, y, 10L))
  base_url <- paste0("/vsicurl/https://dataforgood-fb-data.s3.amazonaws.com/",
                     "forests/v2/global/dinov3_global_chm_v2_ml3/chm/")
  # Crop in 3857 to avoid CRS reprojection surprises.
  pt_3857 <- project(vect(data.frame(x = PT_LON, y = PT_LAT),
                          geom = c("x", "y"), crs = "EPSG:4326"), "EPSG:3857")
  xy <- crds(pt_3857)
  scale_3857 <- 1 / cos(PT_LAT * pi / 180)
  half_u <- HALF_M * scale_3857
  ext_3857 <- ext(xy[1] - half_u, xy[1] + half_u,
                  xy[2] - half_u, xy[2] + half_u)
  chm_tiles <- list()
  for (k in keys) {
    r <- tryCatch(rast(paste0(base_url, k, ".tif")),
                  error = function(e) NULL)
    if (is.null(r)) next
    cr <- tryCatch(crop(r, ext_3857), error = function(e) NULL)
    if (!is.null(cr) && prod(dim(cr)[1:2]) > 0L)
      chm_tiles[[length(chm_tiles) + 1L]] <- cr
  }
  if (length(chm_tiles) == 0L) stop("no CHMv2 tiles")
  chm_native <- if (length(chm_tiles) == 1L) chm_tiles[[1L]]
                else mosaic(sprc(chm_tiles), fun = "mean")
  chm_native <- ifel(chm_native > 60, NA, clamp(chm_native, lower = 0))
  # Aggregate to PALSAR's 25 m: p90 first (matches the production pipeline),
  # then project to WGS84 and resample to hv_db's grid.
  ratio_3857 <- round(25 / res(chm_native)[1])
  chm_agg <- aggregate(chm_native, fact = ratio_3857,
                       fun = function(x, ...) stats::quantile(x, 0.9, na.rm = TRUE))
  chm_wgs <- project(chm_agg, hv_db, method = "bilinear", threads = TRUE)
  chm_wgs
}, error = function(e) { message("CHMv2 load failed: ", e$message); NULL })

p_hv <- ggplot() +
  geom_spatraster(data = hv_db) +
  scale_fill_viridis_c(name = "HV γ° (dB)", option = "magma",
                       na.value = "transparent") +
  geom_point(aes(x = PT_LON, y = PT_LAT), colour = "cyan",
             size = 2, shape = 4) +
  coord_sf(crs = 4326) +
  labs(title = sprintf("ALOS-2 PALSAR-2 HV %d (~25 m)", YEAR)) +
  theme_minimal(base_size = 10)

if (!is.null(chm)) {
  p_chm <- ggplot() +
    geom_spatraster(data = chm) +
    scale_fill_viridis_c(name = "Height (m)", na.value = "transparent",
                         limits = c(0, 20), oob = scales::squish) +
    geom_point(aes(x = PT_LON, y = PT_LAT), colour = "red",
               size = 2, shape = 4) +
    coord_sf(crs = 4326) +
    labs(title = "Meta CHMv2 (DINOv3) p90 @ 25 m") +
    theme_minimal(base_size = 10)
  p <- p_chm + p_hv +
    plot_annotation(
      title = sprintf("5 km window @ (%.4f, %.4f)", PT_LAT, PT_LON),
      subtitle = "CHMv2 misses sub-2 m woody veg; does PALSAR HV pick it up?"
    )
} else {
  p <- p_hv
}

ggsave(OUT_PNG, p, width = 12, height = 5.5, dpi = 150)
message(sprintf("Wrote %s", OUT_PNG))
