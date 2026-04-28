# ebirdabund

**Estimate relative abundance of bird species from raw eBird data.**

`ebirdabund` is an R package that takes raw eBird EBD files and an arbitrary study-region polygon and produces a gridded map of relative abundance. It follows the methodology from the Cornell Lab / eBird Best Practices guide but is designed to run on any region and any species, at scale.

## Installation

```r
# install.packages("pak")
pak::pak("wcornwell/ebirdabund")  # or install from local source:
# install.packages("ebirdabund", repos = NULL, type = "source")
```

Required packages are listed in `DESCRIPTION`. The optional `ebirdst` package is needed for eBird-based range masking.

## Quick start

```r
library(sf)
library(ebirdabund)

# 1. Define your study region (any sf polygon)
nsw <- sf::st_as_sf(
  geodata::gadm("AUS", level = 1, path = "ebirdabund_cache")[
    geodata::gadm("AUS", level = 1, path = "ebirdabund_cache")$NAME_1 == "New South Wales",
  ]
)

# 2. Download and cache habitat covariates — once per region
cov <- prepare_covariates(nsw)

# 3. Fit the model — once per species
model_fit <- fit_species_model(
  polygon      = nsw,
  ebird_zip    = "raw_data/ebd_AU-NSW.txt",
  sampling_txt = "raw_data/ebd_AU-NSW_sampling.txt",
  species      = "Superb Fairywren",
  cov_stack    = cov
)

# 4. Predict and map
result <- predict_species_map(
  model_fit,
  polygon = nsw,
  species = "Superb Fairywren"
)
result$plot
```

Or run the full pipeline in one call:

```r
result <- estimate_abundance(
  polygon      = nsw,
  ebird_zip    = "raw_data/ebd_AU-NSW.txt",
  sampling_txt = "raw_data/ebd_AU-NSW_sampling.txt",
  species      = "Superb Fairywren",
  cov_stack    = cov
)
```

For multiple species use `estimate_abundance_batch()`.

---

## How this compares to `06_abundance.Rmd`

The package shares the same core statistical approach as chapter 6 of the [eBird Best Practices guide](https://cornelllabofornithology.github.io/ebird-best-practices/) but generalises and automates it. The table below summarises the key differences.

### 1. Inputs and data loading

| | `06_abundance.Rmd` | `ebirdabund` |
|---|---|---|
| eBird data | Pre-processed zero-filled CSV | Raw EBD `.txt`/`.zip` + sampling-events `.txt`; zero-fill done internally |
| Species | Wood Thrush (hard-coded) | Any eBird common name |
| Study region | BCR 27 (hard-coded) | Any `sf` polygon |
| Caching | None | Zero-filled data cached as `.rds`; covariate raster cached as GeoTIFF |

The package performs the zero-fill itself via `load_ebird.R`, reading the raw EBD with `data.table::fread` and left-joining complete checklists with species detections. The tutorial assumes this step was completed in an earlier chapter.

### 2. Habitat covariates

| | `06_abundance.Rmd` | `ebirdabund` |
|---|---|---|
| Source | Pre-downloaded MODIS PLAND + elevation CSV | ESA WorldCover, SRTM elevation, WorldClim BIO1/BIO12 — downloaded via `geodata` |
| Variables | 4 species-specific PLAND classes + elevation | 6 generic land-cover fractions (`lc_trees`, `lc_grassland`, `lc_shrubs`, `lc_cropland`, `lc_built`, `lc_water`), `elevation`, `precip_annual`, `temp_annual` |
| Selection | Manually chosen for Wood Thrush | All available; columns with < 4 unique values dropped automatically |
| Raster package | `raster` | `terra` |

### 3. Spatiotemporal subsampling

Both workflows use the same hexagonal subsampling strategy — one checklist per hex cell per week, detections and non-detections sampled independently — via `dggridR` at ~5 km spacing. The package exposes `hex_spacing_km` as a parameter; the tutorial hard-codes 5 km.

Note: the current subsampling groups by `year × week × cell` without stratifying by detection status, whereas the tutorial samples detections and non-detections independently. See `subsample.R` for details.

### 4. Train/test split and model selection

| | `06_abundance.Rmd` | `ebirdabund` |
|---|---|---|
| Train/test split | 80/20 random split | None — all subsampled data used |
| Distributions compared | Zero-inflated Poisson, Negative Binomial, Tweedie | Negative Binomial only |
| Model selection | Spearman rank correlation and MAD on held-out test set | Not performed; NB chosen a priori |

The tutorial compares three distributions and selects the best. The package skips this comparison, fixing negative binomial. This trades flexibility for automation.

### 5. GAM fitting

| | `06_abundance.Rmd` | `ebirdabund` |
|---|---|---|
| Fitting function | `mgcv::gam()` | `mgcv::bam()` with `fREML` + `discrete=TRUE` |
| Smooth `k` values | Fixed: `k=5` for most, `k=7` for cyclic time | Data-driven via `safe_k()`: capped at `n_unique - 1`, minimum 3 |
| Cyclic time knots | `seq(0, 24, length.out = k_time)` | `c(0, 24)` (boundary knots only) |
| Over-fit penalty | None | `gamma = 1.4` |
| Term shrinkage | No | `select = TRUE`; falls back to `select = FALSE` then `discrete = FALSE` on non-convergence |

`bam()` is faster and lower-memory than `gam()` for large datasets, making it more suitable for the package's general-purpose use case.

### 6. Standard-effort prediction

| | `06_abundance.Rmd` | `ebirdabund` |
|---|---|---|
| Peak time | Scans 300 start times; picks time where lower CI is maximised | Fixed at **06:00**; overridable via `peak_time` |
| Reference date | Mid-June (hard-coded for Wood Thrush) | Circular mean of detection DOYs (handles species peaking near year boundary, e.g. Dec–Jan in the southern hemisphere) |
| Standard effort | 1 km, 60 min, 1 observer, Traveling Count | Same |
| Prediction cap | None | Log-scale cap at the 90th percentile of observed non-zero counts, to prevent overflow in data-sparse regions |

### 7. Prediction output and range masking

| | `06_abundance.Rmd` | `ebirdabund` |
|---|---|---|
| Rasterization | `raster::rasterize()` onto a template | `terra::rast()` built from an explicit bounding-box template (prevents NA gaps in concave polygons) |
| Layers | `abd`, `abd_se` | `abd`, `abd_se` |
| Range masking | None | Optional masking to eBirdST range or BirdLife BOTW polygons |
| Map | Base-R `plot()` + `fields::image.plot` legend | `ggplot2` map returned as a list element |

### 8. API and usability

| | `06_abundance.Rmd` | `ebirdabund` |
|---|---|---|
| Interface | Inline script | Exported functions: `prepare_covariates()`, `fit_species_model()`, `predict_species_map()`, `estimate_abundance()` |
| Multi-species | Manual copy-paste | `estimate_abundance_batch()` |
| Progress | Print statements | Numbered step messages, convergence fallbacks, ranked predictor importance table |

### Summary

The package keeps the same core statistical approach — spatiotemporal hex subsampling → negative-binomial GAM → standardised-effort prediction — but generalises it by:

1. Accepting raw eBird files instead of pre-processed CSVs.
2. Downloading covariates automatically from open global sources.
3. Replacing the interactive model-comparison step with a fixed NB choice and automated convergence fallbacks.
4. Using `bam()` for scalability, with data-driven `k` and `gamma = 1.4` to reduce over-fitting.
5. Adding range masking and a prediction cap to improve map quality without user intervention.
