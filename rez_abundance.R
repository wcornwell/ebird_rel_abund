suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(ggplot2)
  library(dplyr)
})

STACK_PATH <- "species_maps/nsw_abundance_stack_3km.tif"
SE_PATH    <- "species_maps/nsw_abundance_se_stack_3km.tif"
REZ_PATH   <- "RenewableEnergyZones_Spatial/RenewableEnergyZones.shp"
CROSSWALK_PATH <- "taxonomy/species_wlab_crosswalk.csv"
TRAITS_PATH    <- "taxonomy/species_traits.csv"
CACHE_DIR  <- "ebirdabund_cache_nsw_buffer"
OUT_DIR    <- "rez_plots"
TOP_N      <- 50

dir.create(OUT_DIR, showWarnings = FALSE)

message("Loading data...")
stack    <- rast(STACK_PATH)
se_stack <- rast(SE_PATH)
rez      <- st_read(REZ_PATH, quiet = TRUE) |>
  st_transform(crs(stack))

# WLAB taxonomy crosswalk, keyed by layer stem (safe_name of the common name)
# so it joins cleanly to snake_case raster layer names.
cw <- read.csv(CROSSWALK_PATH, stringsAsFactors = FALSE)
cw$stem <- gsub("[^a-z0-9]+", "_", tolower(trimws(cw$common_name)))

# Species traits (AVONET migration + feeding guild, EltonTraits stratum),
# keyed the same way. Optional: skip trait columns if the file is absent.
traits <- if (file.exists(TRAITS_PATH)) {
  t <- read.csv(TRAITS_PATH, stringsAsFactors = FALSE)
  t$stem <- gsub("[^a-z0-9]+", "_", tolower(trimws(t$common_name)))
  t
} else NULL

# Pretty-print species names from snake_case layer names
pretty_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  # capitalise first letter only
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

# ── Pre-compute checklist frequency and mean count per species per REZ ─────────
# All cache files share the same checklist pool, so we derive REZ membership
# once from the first available file, then process each species in one pass.
message("Deriving REZ checklist membership from cache...")
first_cache <- list.files(CACHE_DIR, pattern = "zerofilled_.*\\.rds",
                          full.names = TRUE)[1]
chk_base <- readRDS(first_cache)[, c("checklist_id", "latitude", "longitude")]
chk_sf   <- st_as_sf(chk_base, coords = c("longitude", "latitude"),
                     crs = 4326, remove = FALSE)
rez_4326 <- st_transform(rez, 4326)
old_s2   <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)
chk_in_rez <- st_join(chk_sf, rez_4326[, "REZ_Name"], left = FALSE)
sf::sf_use_s2(old_s2)
rez_ids     <- split(chk_in_rez$checklist_id, chk_in_rez$REZ_Name)
n_per_rez   <- lengths(rez_ids)
message(sprintf("  Checklists per REZ: %s",
                paste(names(n_per_rez), n_per_rez, sep = " = ", collapse = ", ")))

sp_keys    <- names(stack)   # snake_case layer names = cache file stems
rez_names  <- rez$REZ_Name

freq_mat <- matrix(NA_real_, nrow = nlyr(stack), ncol = nrow(rez),
                   dimnames = list(sp_keys, rez_names))
cnt_mat  <- matrix(NA_real_, nrow = nlyr(stack), ncol = nrow(rez),
                   dimnames = list(sp_keys, rez_names))

message(sprintf("Computing checklist stats for %d species...", nlyr(stack)))
for (j in seq_len(nlyr(stack))) {
  sp      <- sp_keys[j]
  cache_f <- file.path(CACHE_DIR, paste0("zerofilled_", sp, ".rds"))
  if (!file.exists(cache_f)) next
  sp_data <- readRDS(cache_f)[, c("checklist_id", "observation_count")]
  sp_data$observation_count <- suppressWarnings(as.integer(sp_data$observation_count))

  for (rn in rez_names) {
    ids <- rez_ids[[rn]]
    obs <- sp_data$observation_count[sp_data$checklist_id %in% ids]
    n_det <- sum(obs > 0L, na.rm = TRUE)
    freq_mat[sp, rn] <- n_det / n_per_rez[[rn]]
    cnt_mat[sp, rn]  <- if (n_det > 0L) mean(obs[obs > 0L], na.rm = TRUE) else NA_real_
  }
}
message("  Done.")

message(sprintf("Computing mean abundance for %d species across %d REZs...",
                nlyr(stack), nrow(rez)))

for (i in seq_len(nrow(rez))) {
  rez_name <- rez$REZ_Name[i]
  message(sprintf("  %s", rez_name))

  rez_vect <- vect(rez[i, ])

  cropped <- crop(stack, rez_vect)
  masked  <- mask(cropped, rez_vect)

  # Cells inside the REZ polygon, regardless of per-species range NA.
  # Treat range-masked NAs as zero abundance: a species absent from part of
  # the REZ should pull the REZ mean down, not be excluded from the average.
  inside_r  <- mask(crop(setValues(stack[[1]], 1), rez_vect), rez_vect)
  n_inside  <- as.numeric(global(inside_r, fun = "sum", na.rm = TRUE)[[1]])

  sums  <- global(masked, fun = "sum", na.rm = TRUE)[[1]]
  means <- sums / n_inside
  names(means) <- names(stack)

  # Variance with NA-as-zero: E[x^2] - (E[x])^2 over all REZ cells
  sq_sums <- global(masked^2, fun = "sum", na.rm = TRUE)[[1]]
  sds     <- sqrt(pmax(sq_sums / n_inside - means^2, 0))
  names(sds) <- names(stack)

  # model_se is kept as mean over cells where the species was predicted —
  # it represents model uncertainty in the species' realised range, not in
  # cells that were range-masked to zero.
  se_cropped <- crop(se_stack, rez_vect)
  se_masked  <- mask(se_cropped, rez_vect)
  mean_se    <- global(se_masked, fun = "mean", na.rm = TRUE)[[1]]
  names(mean_se) <- names(stack)

  top_idx <- order(means, decreasing = TRUE)[
    seq_len(min(TOP_N, sum(!is.na(means))))
  ]
  df <- data.frame(
    species    = pretty_name(names(means)[top_idx]),
    abundance  = means[top_idx],
    model_se   = mean_se[top_idx],
    stringsAsFactors = FALSE
  )
  df$species <- factor(df$species, levels = rev(df$species))
  # Asymmetric 95% CI on the log scale, back-transformed
  z <- 1.96
  df$ci_lo <- df$abundance * exp(-z * df$model_se)
  df$ci_hi <- df$abundance * exp( z * df$model_se)

  p <- ggplot(df, aes(x = abundance, y = species)) +
    geom_col(fill = "#2c7bb6") +
    geom_errorbar(
      aes(xmin = ci_lo, xmax = ci_hi),
      orientation = "y", width = 0.4, colour = "grey30", linewidth = 0.4
    ) +
    labs(
      title    = sprintf(
        "Top %d species by mean relative abundance\n%s REZ",
        TOP_N, rez_name
      ),
      x        = "Mean relative abundance (95% CI, log scale)",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank()
    )

  safe_name <- gsub("[^a-z0-9]+", "_", tolower(rez_name))

  out_path  <- file.path(OUT_DIR, sprintf("top%d_%s.png", TOP_N, safe_name))
  ggsave(out_path, p, width = 10, height = 14, dpi = 150)
  message(sprintf("    Saved: %s", out_path))

  cw_idx <- match(names(means), cw$stem)
  tr_idx <- if (!is.null(traits)) match(names(means), traits$stem) else NULL
  tr_col <- function(col) if (is.null(traits)) NA_character_ else traits[[col]][tr_idx]
  all_df <- data.frame(
    species        = pretty_name(names(means)),
    common_name    = cw$common_name[cw_idx],
    scientific_name = cw$scientific_name[cw_idx],
    wlab_taxon_id  = cw$wlab_taxon_id[cw_idx],
    risk_assessed  = cw$risk_assessed[cw_idx],
    wlab_review    = cw$wlab_review[cw_idx],
    migration         = tr_col("migration"),
    garnett_movement  = tr_col("garnett_movement"),
    feeding_guild     = tr_col("feeding_guild"),
    trophic_level     = tr_col("trophic_level"),
    primary_lifestyle = tr_col("primary_lifestyle"),
    foraging_stratum  = tr_col("foraging_stratum"),
    habitat_breadth   = tr_col("habitat_breadth"),
    iucn_status       = tr_col("iucn_status"),
    epbc_status       = tr_col("epbc_status"),
    nsw_status        = tr_col("nsw_status"),
    abundance      = means,
    model_se       = mean_se,
    spatial_cv     = sds / means,
    checklist_freq = freq_mat[names(means), rez_name],
    mean_count     = cnt_mat[names(means), rez_name],
    stringsAsFactors = FALSE
  )
  all_df$ci_lo <- all_df$abundance * exp(-z * all_df$model_se)
  all_df$ci_hi <- all_df$abundance * exp( z * all_df$model_se)
  all_df <- all_df[order(all_df$abundance, decreasing = TRUE), ]

  csv_path <- file.path(OUT_DIR, sprintf("abundance_%s.csv", safe_name))
  write.csv(all_df, csv_path, row.names = FALSE)
  message(sprintf("    Saved: %s", csv_path))
}

message("\nDone. Plots written to: ", OUT_DIR)
