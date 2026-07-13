# patch_seasonality_mean_abund.R
#
# Interim patch, ahead of the full run_batch_nsw.R re-run: fits + predicts
# abundance (with the modal_protocol() fix applied) for the 22 species that
# were wrongly excluded/never-run under the old protocol_type reference-level
# bug, then back-fills their mean_abund/se_abund in the EXISTING
# rez_seasonality/seasonality_all.csv (whose seasonality metrics for these
# species are already correct — dev_expl, is_seasonal, window etc. only
# depend on the relative day-of-year curve shape, not the protocol reference
# level). Does NOT touch species_maps/ or batch_nsw_log.csv, and does NOT
# refresh the 394 already-"ok" species (their mean_abund still carries the
# old +-5-30% bias measured separately; the full batch re-run will fix that).

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf); library(terra); library(geodata); library(yaml)
})

cfg    <- yaml::read_yaml("config.yaml")
CACHE  <- cfg$covariate_cache
ZCACHE <- cfg$zerofill_cache
RAW    <- cfg$raw_data
ebd_files  <- file.path(RAW, cfg$ebd_files)
samp_files <- file.path(RAW, cfg$sampling_files)

REZ_PATH  <- "RenewableEnergyZones_Spatial/RenewableEnergyZones.shp"
SEASO_CSV <- "rez_seasonality/seasonality_all.csv"

STATEWIDE  <- "statewide"
TARGET_REZ <- c(central_west = "Central-West Orana",
                new_england  = "New England",
                south_west   = "South West")
REZ_LEVELS <- c(STATEWIDE, names(TARGET_REZ))

# ── Study + prediction polygons (mirrors run_batch_nsw.R exactly) ────────────
message("Building study polygon (matches production: polygon_pred, no Lord Howe Island)...")
aus       <- geodata::gadm(country = cfg$study_polygon$country, level = 1, path = CACHE)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
PP <- cfg$study_polygon$proj_processing
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_sf, PP), cfg$study_polygon$buffer_metres), 4326)

region_parts <- sf::st_cast(sf::st_geometry(region_sf), "POLYGON", warn = FALSE)
part_xmax    <- vapply(region_parts, function(g) sf::st_bbox(g)[["xmax"]], numeric(1))
region_main  <- sf::st_union(region_parts[part_xmax < 155])
polygon_pred <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_main, PP), cfg$study_polygon$buffer_metres), 4326)
nsw <- sf::st_transform(region_main, 4326)

message("Loading cached covariate stack (should hit cache, no downloads)...")
cov <- prepare_covariates(polygon, cache_dir = CACHE)

tax <- read.csv(cfg$taxonomy_file, stringsAsFactors = FALSE)
sci_lookup <- setNames(tax$scientific_name, tax$common_name)
botw_aliases <- if (file.exists("botw_name_aliases.csv"))
  read.csv("botw_name_aliases.csv", stringsAsFactors = FALSE) else NULL

# ── REZ region vectors (mirrors rez_seasonality_batch.R) ─────────────────────
rez      <- sf::st_read(REZ_PATH, quiet = TRUE) |> sf::st_transform(4326)
region_vect <- c(
  setNames(list(terra::vect(nsw)), STATEWIDE),
  setNames(lapply(names(TARGET_REZ),
                  function(k) terra::vect(rez[rez$REZ_Name == TARGET_REZ[[k]], ])),
           names(TARGET_REZ)))

region_abd <- function(r_pred, region) {
  rv <- region_vect[[region]]
  la <- r_pred[["abd"]]; ls <- r_pred[["abd_se"]]
  cr  <- terra::mask(terra::crop(la, rv), rv)
  ins <- terra::mask(terra::crop(terra::setValues(la, 1), rv), rv)
  nin <- as.numeric(terra::global(ins, "sum", na.rm = TRUE)[[1]])
  se  <- terra::mask(terra::crop(ls, rv), rv)
  c(mean_abund = as.numeric(terra::global(cr, "sum", na.rm = TRUE)[[1]]) / nin,
    se_abund   = as.numeric(terra::global(se, "mean", na.rm = TRUE)[[1]]))
}

SPECIES <- c(
  "Banded Stilt", "Black-tailed Godwit", "Rufous Fieldwren", "Redthroat",
  "Sanderling", "Cinnamon Quail-thrush", "Marbled Frogmouth", "Hooded Plover",
  "Red-browed Pardalote", "Black-bellied Plover", "Australian Bustard",
  "Superb Fruit-Dove", "Wandering Tattler", "Australian Masked-Owl",
  "Australian Painted-Snipe", "Regent Honeyeater", "Wood Sandpiper",
  "Radjah Shelduck", "Australasian Bittern", "Lewin's Rail",
  "Eyrean Grasswren", "Squatter Pigeon"
)

patch_rows <- list()
for (i in seq_along(SPECIES)) {
  sp <- SPECIES[i]
  message(sprintf("\n[%d/%d] %s", i, length(SPECIES), sp))
  res <- tryCatch({
    model_fit <- fit_species_model(
      polygon = polygon_pred, ebird_zip = ebd_files, sampling_txt = samp_files,
      species = sp, cov_stack = cov, cache_dir = ZCACHE, hex_spacing_km = 5)
    pred_out <- predict_species_map(
      model_fit = model_fit, polygon = polygon_pred, species = sp,
      sci_name = unname(sci_lookup[sp]), grid_res_km = 3,
      use_range = TRUE, botw_path = cfg$botw_path, botw_aliases = botw_aliases,
      range_resolution = "27km", border = nsw)
    r_pred <- pred_out$predictions
    rm(model_fit); gc(full = TRUE)
    do.call(rbind, lapply(REZ_LEVELS, function(region) {
      ab <- region_abd(r_pred, region)
      data.frame(common_name = sp, rez = region,
                mean_abund = ab[["mean_abund"]], se_abund = ab[["se_abund"]],
                stringsAsFactors = FALSE)
    }))
  }, error = function(e) {
    message(sprintf("  [FAILED] %s: %s", sp, conditionMessage(e)))
    NULL
  })
  if (!is.null(res)) patch_rows[[sp]] <- res
}

patch_df <- do.call(rbind, patch_rows)
message(sprintf("\nPatched %d species x %d regions = %d rows.",
                length(unique(patch_df$common_name)), length(REZ_LEVELS), nrow(patch_df)))

# ── Apply patch to seasonality_all.csv (with backup) ─────────────────────────
file.copy(SEASO_CSV, paste0(SEASO_CSV, ".bak_pre_protocol_fix_patch"), overwrite = TRUE)
seaso <- read.csv(SEASO_CSV, stringsAsFactors = FALSE)

for (i in seq_len(nrow(patch_df))) {
  idx <- which(seaso$common_name == patch_df$common_name[i] & seaso$rez == patch_df$rez[i])
  if (length(idx) != 1L) {
    message(sprintf("  [WARN] no unique match for %s / %s (n=%d) — skipped",
                    patch_df$common_name[i], patch_df$rez[i], length(idx)))
    next
  }
  seaso$mean_abund[idx] <- patch_df$mean_abund[i]
  seaso$se_abund[idx]   <- patch_df$se_abund[i]
  seaso$run_date[idx]   <- as.character(Sys.Date())
}

write.csv(seaso, SEASO_CSV, row.names = FALSE)
message(sprintf("\nPatched %s in place (backup at %s.bak_pre_protocol_fix_patch)",
                SEASO_CSV, SEASO_CSV))

message("\n== Patched values ==")
print(patch_df, row.names = FALSE)
