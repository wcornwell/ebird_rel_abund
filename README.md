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

## Relationship to Strimas-Mackey et al. (2023)

`ebirdabund` implements the same core statistical workflow as Chapter 6 of *Best Practices for Using eBird Data* (Strimas-Mackey et al. 2023): spatiotemporal hex subsampling → negative-binomial GAM → standardised-effort prediction. The sections below document where the package diverges from that reference implementation.

> Strimas-Mackey, M., W.M. Hochachka, V. Ruiz-Gutierrez, O.J. Robinson, E.T. Miller, T. Auer, S. Kelling, D. Fink, A. Johnston. 2023. *Best Practices for Using eBird Data*. Version 2.0. Cornell Lab of Ornithology, Ithaca, New York. <https://doi.org/10.5281/zenodo.3620739>

### 1. Inputs and data loading

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| eBird data | Pre-processed zero-filled CSV | Raw EBD `.txt`/`.zip` + sampling-events `.txt`; zero-fill performed internally |
| Species | Wood Thrush (hard-coded) | Any eBird common name |
| Study region | BCR 27 (hard-coded) | Any `sf` polygon |
| Caching | None | Zero-filled data cached as `.rds`; covariate raster cached as GeoTIFF |
| Effort filter | None mentioned | `duration_minutes` ∈ [10, 300], `effort_distance_km` ≤ 10, `number_observers` ≤ 10 (eBird Best Practices defaults) |
| Mega-flock records | Retained | Checklists with counts > `max_count` excluded (default 200, parameterised on `clean_ebird()` and `load_ebird()`) |

Zero-fill is performed via `load_ebird.R`, reading the raw EBD with `data.table::fread` and left-joining complete checklists with species detections. Checklists with `observation_count > max_count` are excluded: such records typically reflect targeted sampling of a known flock rather than a passive encounter rate, and their inclusion biases abundance estimates upward for gregarious species. The minimum-duration filter was raised from 5 to 10 minutes following the current eBird Best Practices recommendation; very short checklists tend to under-report species and inflate effort-controlled detection rates.

### 2. Habitat covariates

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| Source | Pre-downloaded MODIS PLAND + elevation CSV | ESA WorldCover, SRTM elevation, WorldClim BIO1/BIO12, WorldPop 2020, JRC Global Surface Water — downloaded automatically via `geodata` / VSICURL |
| Land-cover variables | 4 species-specific PLAND classes | `lc_trees`, `lc_grassland`, `lc_cropland`, `lc_built` (4 generic fractions) |
| Climate | Not included | `precip_annual` (WorldClim BIO12), `temp_annual` (WorldClim BIO1) |
| Water | Not included | `water_occ`: JRC Global Surface Water occurrence 0–100 (Pekel et al. 2016, updated 2021), streamed at 30 m via VSICURL |
| Human footprint | Not included | `pop_density`: log₁₀(persons km⁻² + 1), WorldPop 2020 |
| Soil | Not included | `clay`: SoilGrids 250 m clay fraction 0–5 cm (ISRIC); NA at water set to 0 |
| Canopy height | Not included | `tree_height`: Meta / WRI Canopy Height Maps v2 (DINOv3, 2024 Sentinel-2 imagery), streamed at OVERVIEW_LEVEL=4 (~38 m) and p90-aggregated to the ~1 km master template; pixels >60 m set to NA before aggregation to drop tower/turbine artefacts; NA after aggregation set to 0 |
| Covariate selection | Manually chosen for Wood Thrush | Automated; columns with < 4 unique values dropped; `lc_shrubs` always excluded |
| Raster package | `raster` | `terra` |

`water_occ` is preferred over the ESA WorldCover water fraction because it distinguishes permanent from seasonal surface water at finer spatial resolution. `pop_density` is included to account for the observer-density gradient inherent in citizen-science data. `lc_shrubs` is excluded from model fitting because it is collinear with other land-cover terms and degrades convergence for rare species with sparse detections.

### 3. Spatiotemporal subsampling

Both workflows apply hexagonal subsampling via `dggridR` at ~5 km spacing (one checklist per hex cell per week). The package exposes `hex_spacing_km` as a parameter; Strimas-Mackey et al. (2023) hard-code 5 km.

Note: subsampling groups by `year × week × cell` without stratifying by detection status, whereas Strimas-Mackey et al. (2023) sample detections and non-detections independently. See `subsample.R` for details.

### 4. Model evaluation

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| Train/test split | 80/20 random split used for distribution comparison | No train/test split; all subsampled data used for fitting |
| Distribution comparison | Zero-inflated Poisson, Negative Binomial, Tweedie compared on held-out set | Negative Binomial fixed a priori |
| Cross-validation | Not performed | `evaluate_model_cv()`: k-fold CV reporting Spearman ρ, Pearson r, MAE, RMSE, holdout deviance explained |
| Covariate comparison | Not provided | `compare_covariate_models()`: effort-only vs. full habitat model (and arbitrary additional formulas) via the same CV framework |

Fixing Negative Binomial rather than selecting the distribution from data trades flexibility for automation. `evaluate_model_cv()` and `compare_covariate_models()` provide post-hoc tools for assessing fit quality and diagnosing whether habitat covariates add predictive skill beyond effort alone.

### 5. GAM fitting

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| Fitting function | `mgcv::gam()` | `mgcv::bam()` with `fREML` and `discrete = TRUE` |
| `day_of_year` smooth | Thin-plate spline, `k = 5` | Cyclic cubic spline (`bs = "cc"`), `k = 10`, boundary knots `c(0, 365)` |
| `time_observations_started` smooth | Cyclic cubic spline, `k = 7`, interior knots | Cyclic cubic spline (`bs = "cc"`), `k = 4`, boundary knots `c(0, 24)` |
| Other smooth `k` values | Fixed at 5 | Data-driven via `safe_k()`: min(default, n\_unique − 1), floor of 3 |
| Over-fit penalty | None | `gamma = 1.4` |
| Term shrinkage | No | `select = TRUE`; falls back to `select = FALSE`, then `discrete = FALSE`, on non-convergence |

`bam()` is faster and lower-memory than `gam()` for large datasets. A cyclic spline for `day_of_year` enforces continuity at the year boundary (important for southern-hemisphere species whose abundance peak spans December–January). `gamma = 1.4` penalises complexity beyond what the data support, reducing over-fitting in data-sparse regions.

### 6. Standard-effort prediction

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| Reference date | Mid-June (hard-coded for Wood Thrush) | Circular mean of detection day-of-year |
| Peak observation time | Scans 300 candidate start times; selects time maximising lower confidence limit | Circular mean of detection start times |
| Standard effort | 1 km, 60 min, 1 observer, Traveling Count | Same |
| Prediction cap | None | Log-scale cap at the 90th percentile of observed non-zero counts |
| `abd_se` scale | Response scale | Log (link) scale; asymmetric 95% CI: `[abd × exp(−1.96 × abd_se), abd × exp(+1.96 × abd_se)]` |

The circular mean is used for both reference date and reference time so that species with phenological peaks near the year boundary (e.g. December–January breeders in the southern hemisphere) or near midnight are handled correctly. Storing `abd_se` on the log scale allows asymmetric confidence intervals that respect the positivity constraint on abundance.

### 7. Prediction output and range masking

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| Rasterization | `raster::rasterize()` onto a template | `terra::rast()` from an explicit bounding-box template (prevents NA gaps in concave polygons) |
| Output layers | `abd`, `abd_se` | `abd`, `abd_se` (SE on log scale; see §6) |
| Range masking | None | Optional masking to eBirdST seasonal range or BirdLife BOTW polygons |
| Map | Base-R `plot()` with `fields::image.plot` legend | `ggplot2` map returned as a list element |
| Batch zero-abundance check | N/A | Species with negligible predicted polygon abundance (sum ≤ 10⁻⁵) excluded automatically in `estimate_abundance_batch()` |

### 8. Observer expertise

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| Observer skill | Not modelled | Optional Stage 1 calibration → per-species `s(observer_expertise)` smooth |
| Citation | — | Kelling et al. 2015 (*PLoS ONE*); Johnston et al. 2021 (*Diversity & Distributions*) |

The Cornell tutorial includes `number_observers` (party size on the checklist) but no term for *who* the observer is. Observer skill is known to vary substantially across eBird contributors and is included in the production eBird Status & Trends models (Fink et al.). `ebirdabund` follows the published methodology of Kelling et al. (2015) and Johnston et al. (2021): fit a single global negative-binomial GAM regressing species count per complete checklist on effort covariates with a per-observer random intercept, and use the resulting BLUPs as a continuous per-observer expertise score. The score is then included as `s(observer_expertise)` in each per-species GAM. Predictions are made at the *median* expertise across training data, so the standardised-effort surface represents an average-skill observer.

The Stage 1 calibration runs once via `build_observer_expertise.R` and writes `observer_expertise.rds` into the per-region cache directory. If that file is absent, the per-species pipeline behaves exactly as the tutorial does (no observer term). For a single-species sanity check before generalising, fitting `s(observer_id, bs = "re")` directly into the per-species formula for Superb Fairywren (NSW + 100 km buffer, 387k subsampled checklists, 13k unique observers) dropped per-cell mean predicted abundance by a factor of 0.71 — most of the previously unaccounted observer-skill variance had been bleeding into the habitat smooths. The Stage 1 expertise score is the scalable, citable form of this correction.

### 9. API and usability

| | Strimas-Mackey et al. (2023) | `ebirdabund` |
|---|---|---|
| Interface | Inline R script | Exported functions: `prepare_covariates()`, `fit_species_model()`, `predict_species_map()`, `estimate_abundance()` |
| Multi-species | Manual copy-paste | `estimate_abundance_batch()` |
| Model evaluation | None | `evaluate_model_cv()`, `compare_covariate_models()` |
| Progress | Print statements | Numbered step messages, convergence fallbacks, ranked predictor importance table |

### Summary

The package implements the core pipeline of Strimas-Mackey et al. (2023) — spatiotemporal hex subsampling → negative-binomial GAM → standardised-effort prediction — and generalises it by:

---

## Experiments and methodological tests

A running record of the small evaluations behind the choices encoded above, plus exploratory tests that didn't (or haven't yet) made it into the pipeline. All scripts are committed under `analysis/{topic}/`; outputs under `covariate_tests/{variant}/` or alongside the script. Most use a fixed 19-species test set (`analysis/covariate_test_species.csv`) spanning forest, woodland, open-country, water, urban, arid, alpine, and three deliberately-weak (low dev. expl.) species so deltas are comparable across experiments over time.

Statuses: **Adopted** = encoded in the production pipeline. **Tested, not adopted** = no measurable improvement vs baseline (paired Wilcoxon p > 0.05 or deltas at fold-noise scale). **Pending** = scoped or scripted but not yet run at scale. **Diagnostic** = visual/exploratory, no go/no-go outcome.

### Model family and formulation

| Test | Question | Outcome | Script |
|---|---|---|---|
| Tweedie vs Negative Binomial | Which count family fits eBird better? | **Adopted NB** (chosen a priori; see Modelling Decisions in CLAUDE.md) | `analysis/modeling/test_tweedie_vs_nb.R` |
| Log transforms on effort + habitat | Does log/log1p transforming right-skewed predictors (duration, distance, precip, pop_density) improve held-out skill? | **Adopted** (now baked into `build_gam_formula`) | `analysis/modeling/test_log_transforms.R` |

### Habitat covariates

| Test | Question | Outcome | Script |
|---|---|---|---|
| `pop_density` F-stat scan | Does WorldPop pop_density carry signal across ~40 species, or is it redundant with `lc_built`? | **Adopted** (significant in most species) | `analysis/modeling/test_pop_density.R` |
| Tree-height layer choice | CHMv2 vs OzTreeMap vs ETH-Lang vs GLAD/Potapov? | **CHMv2 adopted** (replaces OzTreeMap; stripe-artefact mitigation tracked in CLAUDE.md) | `analysis/tree_height/test_tree_height_alternatives.R`, `test_chmv2_aggregation.R` |
| CHMv2 two-stage p90 aggregation | Does 38 m → 250 m mean → 1 km p90 beat the current single-step p90 (proposed fix for stripe artefacts)? | **Tested, not adopted** (median ΔSpearman ≈ 0, paired Wilcoxon p > 0.5; see `covariate_tests/tree_height_two_stage_p90/summary.txt`) | `analysis/tree_height/test_tree_height_alternatives.R` |
| CHMv2 layer diagnostics | Visual inspection of the canopy-height raster (native ~1 m, log-scale histograms, outlier provenance) | **Diagnostic** | `analysis/tree_height/check_tree_height_chmv2.R`, `diagnose_chmv2_outliers.R`, `test_chmv2_native_point*.R`, `test_chmv2_small.R` |
| PALSAR HV backscatter | Does ALOS-2 PALSAR HV gamma0 add a sub-2 m woody-vegetation signal CHMv2 cannot resolve? | **Pending** (proof-of-concept + helper + smoke test built; full-scale CV pending) | `analysis/palsar/test_palsar_poc.R`, `test_palsar_helper_smoke.R`, `test_palsar_hv_cv.R` |

### Effort terms and filters

| Test | Question | Outcome | Script |
|---|---|---|---|
| Speed (km/h) as separate effort smooth | Does adding `s(log1p(speed_kph), k=4)` alongside duration & distance smooths improve held-out skill? (eBird S&T uses speed as a distinct effort variable) | **Tested — tiny but consistent improvement, adoption tbd** (19-species median ΔSpearman +0.0003, 17/19 wins, paired Wilcoxon p=0.0002; same direction for DevExpl 16/19, p=0.002 and MAE 12/7, p=0.014. Magnitude is at fold-noise scale.) | `analysis/effort_terms/test_speed_effort.R` |
| Minimum duration: 10 → 5 min | Does lowering the duration threshold from 10 min (current) to 5 min add useful training data, or inflate detection noise? | **Pending** | `analysis/effort_terms/test_min_duration.R` |
| High-count cap sweep | Where should the mega-flock cap sit (sweep 50 → ∞ on three species with different count tails)? | **Adopted cap = 200** | `analysis/count/test_count_threshold.R` |
| X-count checklist exclusion | Should checklists with any X-count entry be dropped from the shared sampling pool *before* zero-fill, so every species sees the same denominator? | **Adopted** (`x_count_ids.rds` is now built once and applied across all species in the buffered cache) | `analysis/count/test_count_threshold.R`, `test_count_filter.R` |

### Observer skill

| Test | Question | Outcome | Script |
|---|---|---|---|
| `s(observer_id, bs="re")` on one species | Does an observer random effect substantively change predictions, or is it noise? | **Adopted as Stage 1 expertise score** — Superb Fairywren mean predicted abundance dropped 1.262 → 0.898 (ratio 0.71); ~30 % of "high numbers" was observer-skill variance bleeding into habitat smooths | `analysis/observer/test_observer_re.R` |

### Integration / regional

| Test | Question | Outcome | Script |
|---|---|---|---|
| Full NSW-only batch | End-to-end smoke test on the NSW polygon with NSW-only training data | **Reference** (original pre-buffer baseline; superseded by buffered run) | `analysis/integration/test_nsw.R` |
| Full NSW + 100 km buffer batch | Same pipeline with ACT/VIC/QLD/SA buffer training data | **Adopted as production setup** (~397 k vs ~133 k checklists per species after subsampling) | `analysis/integration/test_nsw_buffer.R` |
| Cross-region map comparison | Visual side-by-side of NSW-only vs buffered predictions for representative species | **Diagnostic** | `analysis/integration/plot_comparison_nsw.R` |

---



1. Accepting raw eBird EBD files rather than pre-processed CSVs.
2. Downloading covariates automatically from open global sources (ESA WorldCover, SRTM, WorldClim, WorldPop, JRC Global Surface Water).
3. Replacing interactive distribution comparison with a fixed Negative Binomial choice and automated convergence fallbacks.
4. Using `bam()` for scalability, with cyclic splines for periodically-bounded predictors, data-driven `k`, and `gamma = 1.4` to reduce over-fitting.
5. Providing `evaluate_model_cv()` and `compare_covariate_models()` for post-hoc assessment of model skill.
6. Adding range masking, a prediction cap, and a log-scale SE output to improve map quality.
7. Optionally including a per-observer expertise score (Kelling et al. 2015; Johnston et al. 2021) so per-species predictions represent an average-skill observer rather than absorbing observer-skill variance into the habitat smooths.
