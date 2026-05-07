# Workflow Comparison: `06_abundance.Rmd` vs. `ebirdabund` Package

## Overview

`06_abundance.Rmd` is a tutorial from the Cornell Lab / eBird Best Practices guide, illustrating abundance modelling for a single species (Wood Thrush) in BCR 27. The `ebirdabund` package generalises that workflow into a reusable R package designed to run on any region, any species, at scale.

---

## 1. Inputs and Data Loading

| Aspect | `06_abundance.Rmd` | `ebirdabund` package |
|---|---|---|
| eBird data | Pre-processed zero-filled CSV (`ebd_woothr_june_bcr27_zf.csv`) | Raw EBD `.txt`/`.zip` + sampling-events `.txt` — zero-fill performed internally |
| Species | Wood Thrush (hard-coded) | Any eBird common name, passed as a parameter |
| Study region | BCR 27 (hard-coded) | Any `sf` polygon |
| Caching | None | Zero-filled data cached as `.rds`; covariate raster cached as GeoTIFF |

The package performs the zero-fill itself via `load_ebird.R`, reading the raw EBD with `data.table::fread` and left-joining complete checklists with species detections. The tutorial assumes this step has already been done in an earlier chapter.

---

## 2. Habitat Covariates

| Aspect | `06_abundance.Rmd` | `ebirdabund` package |
|---|---|---|
| Source | Pre-downloaded MODIS PLAND land-cover + elevation CSV | ESA WorldCover land cover + SRTM elevation + WorldClim (BIO1, BIO12), all downloaded via `geodata` |
| Variables | 4 species-specific PLAND classes (deciduous broadleaf, mixed forest, cropland, urban) + elevation | 6 generic land-cover fractions (`lc_trees`, `lc_grassland`, `lc_shrubs`, `lc_cropland`, `lc_built`, `lc_water`), `elevation`, `precip_annual`, `temp_annual` |
| Selection | Manually chosen based on Wood Thrush ecology | All available covariates used automatically; columns with < 4 unique values are dropped before model fitting |
| Raster package | `raster` | `terra` |

---

## 3. Spatiotemporal Subsampling

Both workflows use the same hexagonal subsampling strategy (one checklist per hex cell per week, detections and non-detections sampled independently) with `dggridR` at ~5 km spacing. The package exposes `hex_spacing_km` as a parameter; the tutorial hard-codes 5 km.

---

## 4. Train / Test Split and Model Selection

| Aspect | `06_abundance.Rmd` | `ebirdabund` package |
|---|---|---|
| Train/test split | 80 / 20 random split | None — all subsampled data used for fitting |
| Distributions compared | Zero-inflated Poisson, Negative Binomial, Tweedie | Negative Binomial only |
| Model selection | Spearman rank correlation and MAD on held-out test set | Not performed; NB chosen a priori |

The tutorial explicitly compares three distributions and selects the best. The package skips this comparison, fixing the negative binomial as the distribution. This trades flexibility for automation.

---

## 5. GAM Fitting

| Aspect | `06_abundance.Rmd` | `ebirdabund` package |
|---|---|---|
| Fitting function | `mgcv::gam()` | `mgcv::bam()` with `fREML` + `discrete=TRUE` |
| Smooth `k` values | Fixed: `k=5` for most covariates, `k=7` for cyclic time | Data-driven via `safe_k()`: at most 5 for effort terms, 4 for habitat; capped at `n_unique - 1` |
| Cyclic time knots | `seq(0, 24, length.out = k_time)` (full spread) | `c(0, 24)` (boundary knots only) |
| Over-fit penalty | None | `gamma = 1.4` (penalises complexity more strongly) |
| Automatic term shrinkage | No | `select = TRUE` in first attempt; falls back to `select = FALSE` then `discrete = FALSE` if non-convergence |

`bam()` is substantially faster and lower-memory than `gam()` for large datasets, making it more suitable for the package's general-purpose, scalable use case.

---

## 6. Standard-Effort Prediction

| Aspect | `06_abundance.Rmd` | `ebirdabund` package |
|---|---|---|
| Peak time estimation | Scans 300 start times; picks time where LCL is maximised | Fixed at **06:00** (no scanning); overridable via `peak_time` parameter |
| Reference date | `{max_lc_year}-06-15` (mid-June, hard-coded for Wood Thrush) | Circular mean of detection DOYs (handles species peaking near year boundary, e.g. Dec–Jan in the southern hemisphere) |
| Standard effort | 1 km, 60 min, 1 observer, Traveling Count | Same |
| Prediction cap | None | Log-scale cap at the 90th percentile of observed non-zero counts, to prevent overflow in data-sparse extrapolation regions |

---

## 7. Prediction Output and Range Masking

| Aspect | `06_abundance.Rmd` | `ebirdabund` package |
|---|---|---|
| Rasterization | `raster::rasterize()` onto a template | `terra::rast()` built from an explicit bounding-box template (prevents NA gaps in concave polygons) |
| Layers output | `abd`, `abd_se` | `abd`, `abd_se` |
| Range masking | None | Optional masking to eBird status & trends range (`ebirdst::load_ranges`) with fallback to BirdLife BOTW polygons |
| Map output | `raster` + base-R `plot()` with `fields::image.plot` legend | `ggplot2` map returned as a list element |

---

## 8. API and Usability

| Aspect | `06_abundance.Rmd` | `ebirdabund` package |
|---|---|---|
| Interface | Inline script — edit parameters by modifying the document | Three exported functions: `prepare_covariates()`, `fit_species_model()`, `predict_species_map()`, and `estimate_abundance()` as a one-call wrapper |
| Multi-species use | Manual copy-paste | `estimate_abundance_batch()` for running multiple species in sequence |
| Progress feedback | Print statements | Numbered step messages with checklist/deviance summary and ranked predictor importance table |

---

## Summary

The `ebirdabund` package keeps the same core statistical approach (spatiotemporal hex subsampling → negative-binomial GAM → standardised-effort prediction) but makes it general-purpose by:

1. Accepting raw eBird files instead of pre-processed CSVs.
2. Downloading covariates automatically from open global sources.
3. Replacing the interactive model-comparison step with a fixed negative-binomial choice and automated convergence fallbacks.
4. Using `bam()` for scalability, with data-driven `k` and `gamma=1.4` to reduce over-fitting.
5. Adding range masking and a prediction cap to improve map quality without user intervention.
