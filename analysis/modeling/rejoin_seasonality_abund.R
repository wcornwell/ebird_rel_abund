# rejoin_seasonality_abund.R
#
# Cheap re-join (no model fitting): refresh mean_abund / se_abund for EVERY row
# of rez_seasonality/seasonality_all.csv directly from the current post-protocol-
# fix abundance stack (species_maps/nsw_abundance_stack_3km.tif, rebuilt after
# commit 0f6e208), and add the new range_masked_in_rez flag. Supersedes the
# earlier per-species refit patch (patch_seasonality_mean_abund.R) now that the
# full stack is up to date.
#
# The seasonality metrics (is_seasonal, peak, window, dev_expl, ...) are NOT
# touched — they depend on the day-of-year curve shape, not the abundance stack.
# region_abd() here mirrors the updated version in rez_seasonality_batch.R:
#   * species absent from the stack (not modelled)      -> NA / NA / NA
#   * region wholly outside the range mask (all-NA crop)-> 0  / 0  / TRUE
#   * otherwise                                         -> mean / meanSE / FALSE

suppressPackageStartupMessages({
  library(sf); library(terra); library(geodata); library(yaml)
})
sf::sf_use_s2(FALSE)

cfg   <- yaml::read_yaml("config.yaml")
CACHE <- cfg$covariate_cache
STACK <- "species_maps/nsw_abundance_stack_3km.tif"
SE    <- "species_maps/nsw_abundance_se_stack_3km.tif"
REZ_PATH  <- "RenewableEnergyZones_Spatial/RenewableEnergyZones.shp"
SEASO_CSV <- "rez_seasonality/seasonality_all.csv"

STATEWIDE  <- "statewide"
TARGET_REZ <- c(central_west = "Central-West Orana",
                new_england  = "New England",
                south_west   = "South West")

safe_name <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# ── NSW-land polygon (mainland only, drops Lord Howe) + REZ region vectors ────
message("Building region polygons...")
aus       <- geodata::gadm(country = cfg$study_polygon$country, level = 1, path = CACHE)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
parts <- sf::st_cast(sf::st_geometry(region_sf), "POLYGON", warn = FALSE)
xmax  <- vapply(parts, function(g) sf::st_bbox(g)[["xmax"]], numeric(1))
nsw   <- sf::st_transform(sf::st_union(parts[xmax < 155]), 4326)

rez <- sf::st_read(REZ_PATH, quiet = TRUE) |> sf::st_transform(4326)
region_vect <- c(
  setNames(list(terra::vect(nsw)), STATEWIDE),
  setNames(lapply(names(TARGET_REZ),
                  function(k) terra::vect(rez[rez$REZ_Name == TARGET_REZ[[k]], ])),
           names(TARGET_REZ)))

message("Loading abundance + SE stacks...")
abd_stack <- terra::rast(STACK)
se_stack  <- terra::rast(SE)

region_abd <- function(stem, region) {
  if (!stem %in% names(abd_stack))
    return(c(mean_abund = NA_real_, se_abund = NA_real_, range_masked_in_rez = NA_real_))
  rv  <- region_vect[[region]]
  la  <- abd_stack[[stem]]; ls <- se_stack[[stem]]
  cr  <- terra::mask(terra::crop(la, rv), rv)
  ins <- terra::mask(terra::crop(terra::setValues(la, 1), rv), rv)
  nin <- as.numeric(terra::global(ins, "sum", na.rm = TRUE)[[1]])
  if (is.na(nin) || nin == 0)
    return(c(mean_abund = NA_real_, se_abund = NA_real_, range_masked_in_rez = NA_real_))
  n_modelled <- as.numeric(terra::global(!is.na(cr), "sum", na.rm = TRUE)[[1]])
  if (n_modelled == 0)
    return(c(mean_abund = 0, se_abund = 0, range_masked_in_rez = 1))
  se  <- terra::mask(terra::crop(ls, rv), rv)
  c(mean_abund = as.numeric(terra::global(cr, "sum", na.rm = TRUE)[[1]]) / nin,
    se_abund   = as.numeric(terra::global(se, "mean", na.rm = TRUE)[[1]]),
    range_masked_in_rez = 0)
}

# ── Recompute every row from the stack ────────────────────────────────────────
seaso <- read.csv(SEASO_CSV, stringsAsFactors = FALSE, check.names = FALSE)
file.copy(SEASO_CSV, paste0(SEASO_CSV, ".bak_pre_rejoin"), overwrite = TRUE)
message(sprintf("Re-joining abundance for %d rows...", nrow(seaso)))

# one region_abd call per (species, region); cache per species to avoid re-crop
stems  <- safe_name(seaso$common_name)
uniq   <- unique(data.frame(stem = stems, rez = seaso$rez, stringsAsFactors = FALSE))
vals   <- t(mapply(region_abd, uniq$stem, uniq$rez))
rownames(vals) <- paste(uniq$stem, uniq$rez, sep = "\r")
key    <- paste(stems, seaso$rez, sep = "\r")

seaso$mean_abund          <- as.numeric(vals[key, "mean_abund"])
seaso$se_abund            <- as.numeric(vals[key, "se_abund"])
rmask                     <- as.logical(vals[key, "range_masked_in_rez"])

# insert range_masked_in_rez immediately after se_abund, preserving all other cols
cols  <- names(seaso)
if (!"range_masked_in_rez" %in% cols) {
  seaso$range_masked_in_rez <- rmask
  pos   <- match("se_abund", cols)
  order <- append(cols, "range_masked_in_rez", after = pos)
  seaso <- seaso[, order]
} else {
  seaso$range_masked_in_rez <- rmask
}

write.csv(seaso, SEASO_CSV, row.names = FALSE)

# ── Summary ───────────────────────────────────────────────────────────────────
message(sprintf("\nWrote %s (backup: %s.bak_pre_rejoin)", SEASO_CSV, SEASO_CSV))
message(sprintf("  mean_abund NA (not modelled): %d", sum(is.na(seaso$mean_abund))))
message(sprintf("  range_masked_in_rez TRUE (masked -> 0): %d",
                sum(seaso$range_masked_in_rez %in% TRUE)))
message(sprintf("  se_abund max: %.4g (was 3022.9 pre-fix)", max(seaso$se_abund, na.rm = TRUE)))
