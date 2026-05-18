#!/usr/bin/env Rscript
# CV experiment: does adding PALSAR HV gamma0 as a habitat covariate improve
# held-out prediction for species that depend on sub-2 m woody vegetation
# (which CHMv2 cannot resolve)?
#
# Pattern follows analysis/test_log_transforms.R:
#   - one call to fit_species_model() per species to get the shared $data frame
#   - two formulas: base (current production set) and +s(palsar_hv, k=4)
#   - 5-fold CV on each formula via evaluate_model_cv()
#
# Output: analysis/palsar_hv_cv_results.csv + a comparison plot.

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(data.table)
})

source("analysis/palsar_helpers.R")

# ---- Config ---------------------------------------------------------------
EBD_FILES <- c(
  "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt",
  "ebirdabund/raw_data/ebd_AU-VIC_unv_smp_relMar-2026/ebd_AU-VIC_unv_smp_relMar-2026.txt",
  "ebirdabund/raw_data/ebd_AU-QLD_unv_smp_relMar-2026/ebd_AU-QLD_unv_smp_relMar-2026.txt",
  "ebirdabund/raw_data/ebd_AU-SA_unv_smp_relMar-2026/ebd_AU-SA_unv_smp_relMar-2026.txt",
  "ebirdabund/raw_data/ebd_AU-ACT_unv_smp_relMar-2026/ebd_AU-ACT_unv_smp_relMar-2026.txt"
)
SAMP_FILES <- sub("\\.txt$", "_sampling.txt", EBD_FILES)
CACHE       <- "ebirdabund_cache"
ZEROFILL    <- "ebirdabund_cache_nsw_buffer"
PALSAR_YEAR <- 2020L

# Expect HELP (sub-2 m shrub specialists)
SHRUB_SPECIES <- c(
  "White-fronted Chat",          # chenopod / saltmarsh, low (<1 m) shrubs
  "Crimson Chat",                # chenopod, arid low shrub
  "Orange Chat",                 # saltbush / samphire
  "Spotted Quail-thrush",        # leaf litter + shrub understorey
  "Chestnut-rumped Heathwren",   # dense low heath
  "Shy Heathwren",               # dense low heath / mallee understorey
  "Southern Scrub-Robin"         # low mallee shrub
)
# Expect NO HELP (controls — open-country grass or tall-canopy birds)
CONTROL_SPECIES <- c(
  "White-throated Treecreeper",  # wet-forest canopy specialist
  "Brown Treecreeper",           # open woodland, tree-trunk forager
  "Australian Pipit"             # open grassland / bare ground
)
SPECIES <- c(SHRUB_SPECIES, CONTROL_SPECIES)

# ---- Build cov stack ------------------------------------------------------
message("Loading NSW + 100 km buffer polygon...")
aus <- geodata::gadm("AUS", 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
nsw_buf <- sf::st_buffer(sf::st_transform(nsw, 3577), 100000) |>
  sf::st_transform(4326)

message("Loading standard covariate stack (cached)...")
cov <- prepare_covariates(nsw_buf, cache_dir = CACHE)
message("  Layers: ", paste(names(cov), collapse = ", "))

# Pull PALSAR HV for the buffered bbox, aligned to cov.
bb_buf <- sf::st_bbox(nsw_buf)
ext_buf <- ext(bb_buf[1], bb_buf[3], bb_buf[2], bb_buf[4])
message(sprintf("Loading PALSAR HV (%d) for bbox [%.1f, %.1f, %.1f, %.1f]...",
                PALSAR_YEAR, bb_buf[1], bb_buf[2], bb_buf[3], bb_buf[4]))
palsar_hv <- load_palsar_hv(
  bb       = as.numeric(bb_buf),
  ext      = ext_buf,
  template = cov[[1]],
  cache_dir = CACHE,
  year     = PALSAR_YEAR
)
names(palsar_hv) <- "palsar_hv"

# Stash a summary so we can sanity-check the values
v <- values(palsar_hv); v <- v[!is.na(v)]
message(sprintf("  HV gamma0: median %.1f, p10 %.1f, p90 %.1f dB (n=%s)",
                median(v), quantile(v, 0.1), quantile(v, 0.9),
                format(length(v), big.mark = ",")))

# Append to standard stack.
cov_aug <- c(cov, palsar_hv)
message("Augmented stack layers: ", paste(names(cov_aug), collapse = ", "))

# ---- Per-species CV -------------------------------------------------------
results <- list()

for (sp in SPECIES) {
  message("\n==============================================")
  message("Species: ", sp)
  message("==============================================")

  fit <- tryCatch(
    fit_species_model(
      polygon       = nsw_buf,
      ebird_zip     = EBD_FILES,
      sampling_txt  = SAMP_FILES,
      species       = sp,
      cov_stack     = cov_aug,
      cache_dir     = ZEROFILL
    ),
    error = function(e) {
      message("  fit_species_model failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(fit)) next

  df <- fit$data
  n_orig <- nrow(df)
  # fit_species_model uses a sampling-master fast path that pre-attaches the
  # production covariates and skips extraction when those are already present.
  # palsar_hv isn't in that fast path, so we attach it post-hoc by sampling at
  # the same lat/lon. Bilinear matches the production extract_covariates path.
  if (!"palsar_hv" %in% names(df)) {
    pts <- terra::vect(
      data.frame(x = df$longitude, y = df$latitude),
      geom = c("x", "y"), crs = "EPSG:4326"
    )
    hv_vals <- terra::extract(palsar_hv, pts, method = "bilinear")[, 2]
    df$palsar_hv <- hv_vals
  }
  # Drop checklists missing palsar_hv (mostly coastal/water) so the base and
  # +HV models are compared on the same observations.
  n_na <- sum(is.na(df$palsar_hv))
  df   <- df[!is.na(df$palsar_hv), ]
  n_det <- sum(df$observation_count > 0)
  message(sprintf("  n checklists: %d (dropped %d with NA HV) | detections: %d",
                  nrow(df), n_na, n_det))
  if (n_det < 50L) {
    message("  Skipping (too few detections)")
    next
  }

  hab_cols_base <- grep("^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height)",
                        names(df), value = TRUE)
  hab_cols_base <- setdiff(hab_cols_base, "lc_shrubs")
  hab_cols_base <- hab_cols_base[vapply(hab_cols_base, function(col) {
    length(unique(stats::na.omit(df[[col]]))) >= 4L
  }, logical(1))]

  formula_base <- build_gam_formula(df, hab_cols_base)

  # Augmented formula: base + s(palsar_hv, k=4).
  formula_aug <- stats::update(
    formula_base,
    sprintf(". ~ . + s(palsar_hv, k = %d)",
            safe_k(df$palsar_hv, 4L))
  )

  message("Running 5-fold CV (base)...")
  cv_base <- tryCatch(
    evaluate_model_cv(df, formula = formula_base, k = 5L),
    error = function(e) { message("  CV failed: ", conditionMessage(e)); NULL }
  )
  message("Running 5-fold CV (+palsar_hv)...")
  cv_aug <- tryCatch(
    evaluate_model_cv(df, formula = formula_aug, k = 5L),
    error = function(e) { message("  CV failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(cv_base) || is.null(cv_aug)) next

  message(sprintf(
    "  Spearman r  base: %.4f  +HV: %.4f  [%+.4f]",
    cv_base$summary["spearman_r"], cv_aug$summary["spearman_r"],
    cv_aug$summary["spearman_r"] - cv_base$summary["spearman_r"]
  ))
  message(sprintf(
    "  dev expl    base: %.4f  +HV: %.4f  [%+.4f]",
    cv_base$summary["holdout_dev_expl"], cv_aug$summary["holdout_dev_expl"],
    cv_aug$summary["holdout_dev_expl"] - cv_base$summary["holdout_dev_expl"]
  ))
  message(sprintf(
    "  RMSE        base: %.4f  +HV: %.4f  [%+.4f]",
    cv_base$summary["rmse"], cv_aug$summary["rmse"],
    cv_aug$summary["rmse"] - cv_base$summary["rmse"]
  ))

  results[[sp]] <- data.frame(
    species         = sp,
    group           = if (sp %in% SHRUB_SPECIES) "shrub" else "control",
    n               = nrow(df),
    n_det           = n_det,
    spearman_base   = unname(cv_base$summary["spearman_r"]),
    spearman_aug    = unname(cv_aug$summary["spearman_r"]),
    spearman_delta  = unname(cv_aug$summary["spearman_r"] -
                             cv_base$summary["spearman_r"]),
    devexpl_base    = unname(cv_base$summary["holdout_dev_expl"]),
    devexpl_aug     = unname(cv_aug$summary["holdout_dev_expl"]),
    devexpl_delta   = unname(cv_aug$summary["holdout_dev_expl"] -
                             cv_base$summary["holdout_dev_expl"]),
    rmse_base       = unname(cv_base$summary["rmse"]),
    rmse_aug        = unname(cv_aug$summary["rmse"]),
    rmse_delta      = unname(cv_aug$summary["rmse"] -
                             cv_base$summary["rmse"])
  )
}

if (length(results) == 0L) stop("No species succeeded")

df_res <- rbindlist(results)
write.csv(df_res, "analysis/palsar_hv_cv_results.csv", row.names = FALSE)
message("\nResults: analysis/palsar_hv_cv_results.csv")
print(df_res)

# ---- Comparison plot ------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
df_plot <- df_res
df_plot$species <- factor(df_plot$species, levels = df_plot$species[order(df_plot$spearman_delta)])
p <- ggplot(df_plot, aes(x = species, y = spearman_delta, fill = group)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(shrub = "#cc4c02", control = "#41ab5d")) +
  coord_flip() +
  labs(
    title = "Effect of adding s(palsar_hv) to the GAM",
    subtitle = "Δ 5-fold CV Spearman ρ vs base formula. Positive = HV helps.",
    x = NULL, y = "Δ Spearman ρ (CV)"
  ) +
  theme_minimal(base_size = 11)
ggsave("analysis/palsar_hv_cv_plot.png", p, width = 8, height = 5, dpi = 150)
message("Plot:    analysis/palsar_hv_cv_plot.png")
