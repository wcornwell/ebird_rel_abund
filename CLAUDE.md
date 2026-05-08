# ebird_rel_abund — Project Guide for Claude Code

## What This Project Does

Estimates relative abundance for all bird species recorded in NSW, Australia, using eBird observation data. For each species it fits a negative-binomial GAM that controls for observation effort, then predicts a gridded abundance surface across NSW. Outputs feed into a Renewable Energy Zone (REZ) analysis that ranks bird species by predicted abundance within each energy zone.

The core modelling logic lives in an R package (`ebirdabund/`) that is region-agnostic. The NSW-specific orchestration is in top-level scripts. **The medium-term goal is to scale this to broader geographic scope** — keep that in mind when adding features (prefer parameterisation over hard-coding).

---

## Repo Layout

```
ebird_rel_abund/
├── ebirdabund/              # R package — all reusable modelling logic
│   ├── R/
│   │   ├── load_ebird.R     # Read/filter/zero-fill eBird data
│   │   ├── covariates.R     # Download, cache, and extract habitat layers
│   │   ├── subsample.R      # Spatiotemporal hex-cell subsampling
│   │   ├── model.R          # GAM formula builder + fitter (fit_gam)
│   │   ├── estimate_abundance.R  # High-level API: fit_species_model, predict_species_map
│   │   ├── predict.R        # Grid prediction + log-scale SE output
│   │   ├── batch.R          # Parallel multi-species orchestration
│   │   ├── evaluate.R       # k-fold CV: evaluate_model_cv, compare_covariate_models
│   │   └── utils.R          # Plotting, range loading, time/name helpers
│   └── tests/testthat/      # Unit tests per module
│
├── run_batch_nsw.R          # MAIN PIPELINE — full NSW batch run
├── nsw_species_list.R       # Generate nsw_species_list.csv from raw EBD
├── rez_abundance.R          # Post-processing: per-REZ abundance stats + plots
│
├── analysis/                # Exploratory and evaluation scripts
│   ├── test_tweedie_vs_nb.R     # CV comparison: Tweedie vs NB family
│   ├── test_log_transforms.R    # CV comparison: linear vs log-transformed predictors
│   ├── test_pop_density.R       # F-stat analysis of pop_density across 40 species
│   ├── test_count_filter.R      # Visual check of smooths/predictions for 4 species
│   ├── test_nsw_buffer.R        # Integration test for multi-state buffered run
│   ├── test_nsw.R               # Original single-state test (reference only)
│   ├── 06_abundance.Rmd         # Cornell eBird Best Practices tutorial (reference)
│   └── workflow_comparison.md   # Comparison of tutorial vs ebirdabund package
│
├── ebirdabund_cache/        # Auto-generated cache (do not commit)
│   ├── zerofilled_*.rds     # Per-species zero-filled checklists (~490 files)
│   ├── cov_stack_v3_*.tif   # Covariate raster stack (bbox + version keyed)
│   └── gadm/, climate/, elevation/, landuse/, population/  # Downloaded tiles
│
├── species_maps/            # Per-species output (do not commit)
│   ├── 3km/*.tif            # abd + abd_se layers at 3 km
│   ├── 9km/*.tif            # abd + abd_se layers at 9 km
│   ├── models/*.rds         # Fitted GAM objects
│   ├── smooths/*.png        # Partial-effect plots
│   └── nsw_abundance_stack_{3,9}km.tif   # Multi-band stacks (all species)
│
├── rez_plots/               # REZ analysis outputs
│   ├── top50_{rez_name}.png
│   └── abundance_{rez_name}.csv
│
├── nsw_species_list.csv     # 556 species, reporting_rate >= 0.5%
├── nsw_ebird_taxonomy.csv   # common_name → scientific_name for BOTW lookup
├── batch_nsw_log.csv        # Per-species run log — columns: common_name,
│                            #   scientific_name, status, exclusion_reason,
│                            #   run_date, n_checklists, n_positive, dev_expl,
│                            #   model_sum, peak_doy, peak_time, max_obs_count,
│                            #   max_modeled_abd, range_source, spatial_cv,
│                            #   error_message
└── botw_species/BOTW_2025.gpkg  # BirdLife range polygons (9.3 GB, not committed)
```

---

## Data Inputs

| File | Size | Purpose |
|------|------|---------|
| `ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt` | 6.1 GB | NSW eBird observations (EBD) |
| `ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt` | 316 MB | NSW complete checklists |
| `ebirdabund/raw_data/ebd_AU-ACT_unv_smp_relMar-2026/ebd_AU-ACT_unv_smp_relMar-2026.txt` | — | ACT eBird observations |
| `ebirdabund/raw_data/ebd_AU-ACT_unv_smp_relMar-2026/ebd_AU-ACT_unv_smp_relMar-2026_sampling.txt` | — | ACT complete checklists |
| `ebirdabund/raw_data/ebd_AU-VIC_unv_smp_relMar-2026/ebd_AU-VIC_unv_smp_relMar-2026.txt` | 6.9 GB | VIC eBird observations |
| `ebirdabund/raw_data/ebd_AU-VIC_unv_smp_relMar-2026/ebd_AU-VIC_unv_smp_relMar-2026_sampling.txt` | — | VIC complete checklists |
| `ebirdabund/raw_data/ebd_AU-QLD_unv_smp_relMar-2026/ebd_AU-QLD_unv_smp_relMar-2026.txt` | 7.1 GB | QLD eBird observations |
| `ebirdabund/raw_data/ebd_AU-QLD_unv_smp_relMar-2026/ebd_AU-QLD_unv_smp_relMar-2026_sampling.txt` | — | QLD complete checklists |
| `ebirdabund/raw_data/ebd_AU-SA_unv_smp_relMar-2026/ebd_AU-SA_unv_smp_relMar-2026.txt` | — | SA eBird observations |
| `ebirdabund/raw_data/ebd_AU-SA_unv_smp_relMar-2026/ebd_AU-SA_unv_smp_relMar-2026_sampling.txt` | — | SA complete checklists |
| `botw_species/BOTW_2025.gpkg` | 9.3 GB | BirdLife range polygons for masking |
| `nsw_ebird_taxonomy.csv` | 23 KB | Common → scientific name lookup |
| `nsw_species_list.csv` | 37 KB | Species to model (pre-filtered) |
| `RenewableEnergyZones_Spatial/` | — | NSW REZ polygons for downstream analysis |

The EBD is read with integer column indices (cols 6, 11, 35) not names — `validate_ebd_header()` checks these against the expected column names before any `fread` call. If you update to a different EBD release, run a test species first and watch for header-validation errors.

---

## Pipeline Flow

```
1. nsw_species_list.R
   → reads NSW sampling + EBD, computes reporting rates
   → writes nsw_species_list.csv  (run once per EBD update)

2. run_batch_nsw.R  (main run)
   a. POLYGON: NSW state boundary + 100 km buffer (GDA94/Albers, EPSG:3577)
      — includes all of ACT and border regions of VIC/QLD/SA as training data
   b. PRE-CACHE: reads all 5 state sampling files filtered to buffered bbox,
      then scans all 5 EBD files → zerofilled_{species}.rds per species in
      ebirdabund_cache_nsw_buffer/  (prevents parallel workers competing for files)
   c. EBIRDST RANGES: for each species in the ebirdst 2023 catalogue, downloads
      range .gpkg files (smooth + raw, 27 km + 9 km) to the local ebirdst data
      directory (~11 s/species, ~48 min first run). Skips already-downloaded
      species, so subsequent runs complete instantly. Requires EBIRDST_KEY env var
      (set automatically in the script from the stored key).
   d. COVARIATES: download ESA/SRTM/WorldClim/WorldPop/JRC once →
      cov_stack_v4_{bbox}.tif in ebirdabund_cache/ (bbox-keyed, shared)
   e. BATCH: estimate_abundance_batch() — parallel GAM per species
      For each species:
        load_ebird() [uses cache] → extract_covariates() → subsample_hex()
        → fit_gam() → predict_abundance() [two resolutions] → save .tif + .png
        Map PNGs include NSW state border overlay for reference.
   f. STACK: terra::rast() all .tif → nsw_abundance_stack_{res}.tif
   g. LOG: append to batch_nsw_log.csv after each chunk

3. rez_abundance.R  (post-processing, run after batch)
   → loads stack + REZ polygons + zerofilled cache (ebirdabund_cache_nsw_buffer/)
   → per REZ: mean abundance, checklist frequency, mean count
   → writes top50_{rez}.png + abundance_{rez}.csv
```

**Two cache directories:**
- `ebirdabund_cache/` — covariate tiles, GADM boundaries. Shared across runs; do not delete.
- `ebirdabund_cache_nsw_buffer/` — per-species zerofilled `.rds` files, keyed to the NSW+buffer sampling pool. Delete and rebuild if the study polygon or EBD files change.

To re-run from scratch: delete `species_maps/`, `ebirdabund_cache_nsw_buffer/`, and `batch_nsw_log.csv`. Keep `ebirdabund_cache/` — the covariate stack and static downloads don't change.

---

## Key Modelling Decisions

**GAM family**: Negative binomial (`mgcv::nb()`), chosen a priori via `test_tweedie_vs_nb.R`.

**Effort covariates** (control for observation bias):
- `day_of_year` — cyclic cubic spline, k=10, knots c(0, 365)
- `time_observations_started` — cyclic cubic spline, k=10, knots c(0, 24)
- `log(duration_minutes)` — k=4
- `log1p(effort_distance_km)` — k=4 (log1p because stationary counts = 0)
- `number_observers` — k=4
- `protocol_type` — parametric factor

**Habitat covariates** (abundance signal):
- Land cover fractions from ESA WorldCover: `lc_trees`, `lc_grassland`, `lc_cropland`, `lc_built` — k=4 (excludes `lc_shrubs`, collinear)
- `elevation` — k=6
- `log(precip_annual)` — k=6 (log: ~1 order of magnitude range in NSW)
- `temp_annual` — k=6
- `log1p(pop_density)` — k=4 (log1p: spans ~4 orders of magnitude)
- `water_occ` — k=4
- `clay` — k=4 (SoilGrids 0–5 cm clay fraction; NA at ocean/water set to 0)
- `log1p(tree_height)` — k=4 (OzTreeMap 30 m canopy height median, Pucino et al. 2025; NA where no canopy set to 0)

Effort covariates are log-transformed because they are right-skewed and the log scale distributes knots more evenly. The transformation happens inside the GAM formula string so raw column values flow through unchanged — the prediction surface (`duration_minutes = 60`, `effort_distance_km = 1`) is evaluated correctly.

**Fitting**: `mgcv::bam()` with `fREML`, `discrete=TRUE`, `gamma=1.4`, `select=TRUE`. Falls back to `select=FALSE` then `discrete=FALSE` on convergence failure.

**Prediction**: Standardised effort (60 min, 1 km, 1 observer, Traveling Count) at the circular mean DOY and time of detection checklists. Abundance is capped at the 90th percentile of non-zero training counts (log scale) to prevent extrapolation overflow.

**Exclusion criteria**: fewer than 50 positive checklists after hex subsampling, or model sum ≤ 1e-5 across the polygon.

**Range masking**: ebirdst 2023 preferred; falls back to BOTW 2025 if the species is not in the ebirdst catalogue or if the ebirdst mask produces zero abundance over the study region. ~266 of 479 modelled NSW species are in the ebirdst 2023 catalogue. Ranges are pre-downloaded by `run_batch_nsw.R` before the batch (download step is skipped for species already cached locally). The `range_source` column in `batch_nsw_log.csv` records which source was used per species.

---

## Evaluation

Use `evaluate_model_cv(fit$data, k=5)` for held-out metrics (Spearman ρ, Pearson r, RMSE, holdout deviance explained). The evaluation scripts follow the pattern in `test_tweedie_vs_nb.R`:

1. Call `fit_species_model()` to get the subsampled data frame.
2. Build alternative formulas manually if testing a variant.
3. Pass each formula to `evaluate_model_cv()`.
4. Compare the `$summary` vectors.

`compare_covariate_models(fit$data)` runs effort-only vs. full-habitat in one call and prints a table.

---

## Possible Improvements (flagged during review)

### High priority
- **X-count checklists**: Currently excluded in `clean_ebird()` *after* zero-fill. They should be excluded from the sampling pool *before* zero-fill so all species share the same checklist universe. Otherwise, a species' absence rate is computed against a slightly inflated denominator. Requires re-building the zerofilled cache.
- **Config file**: `run_batch_nsw.R` still has hard-coded paths at the top (EBD files, cache dirs, output dirs). Extracting these into a `config.R` or YAML would make it easier to point the same script at a different region.
- **EBD column index brittleness**: ~~Columns 6, 11, 35 are hard-coded.~~ **Done** — `validate_ebd_header()` now checks column names at those positions before any `fread` call.

### Medium priority
- **Generalise `nsw_species_list.R`**: It's NSW-specific (variable names, output file name, plot title). Wrapping it into a function `generate_species_list(region_name, ebd_path, samp_path, ...)` would make it reusable for geographic scaling.
- **Taxonomy file**: `nsw_ebird_taxonomy.csv` is derived from the NSW EBD. For other regions, this needs to be generated from the regional EBD or replaced with the full eBird taxonomy (available from eBird as `eBird_taxonomy_v*.csv`).
- **Covariate cache versioning**: The cache key uses `v4` as a hard-coded version tag. If covariate layers are updated, this needs a manual bump. A hash of the source URLs or a date stamp would be more robust.
- **REZ script parameterisation**: `rez_abundance.R` paths and the `TOP_N = 50` cutoff are hard-coded. Making it accept command-line arguments would make it easier to run for new regions or different analysis polygons.
- **Parallel memory**: `n_cores = detectCores() - 1` doesn't account for terra's per-worker memory usage. For larger regions (more grid cells, higher resolution), this may OOM. Consider a `max_ram_gb` guard that caps workers.

### Lower priority
- **Incremental stack building**: The stack-building step at the end of `run_batch_nsw.R` re-reads all TIFs. For very large species lists, this is slow. Could build incrementally or defer to a separate script.
- **Smooth plotting**: `plot_gam_smooths()` in `utils.R` saves a static PNG. Interactive HTML (e.g. via `plotly`) would make diagnosing smooth shapes much faster.
- **Minimum duration filter**: Currently 5 minutes. eBird best practices suggest 10 minutes is more reliable. Worth testing via CV before changing.

---

## Geographic Scaling Plan

The package (`ebirdabund/`) is region-agnostic — `fit_species_model()` and `estimate_abundance_batch()` accept any `sf` polygon and any number of EBD/sampling files as character vectors.

**Current setup:** NSW + 100 km buffer as the study polygon, with EBD training data from NSW, ACT, VIC, QLD, and SA. The buffer captures ~397k checklists per species after hex subsampling (vs ~133k with NSW-only data).

The work needed to add a new region:

1. Download regional EBD + sampling files.
2. Run `generate_species_list()` (once generalised) for the new region.
3. Generate a taxonomy CSV for the region (or use the full eBird taxonomy).
4. Create a `run_batch_{region}.R` pointing at the new files — or, better, a single `run_batch.R` that reads a config file.
5. The covariate downloads (ESA, WorldClim, etc.) are global and will re-use cached tiles if the bounding box overlaps.
6. Range masking is already global: BOTW covers all species, and ebirdst covers whichever species it models. Run the ebirdst pre-download step for the new species list before the batch — the download key (`EBIRDST_KEY`) is already set in `run_batch_nsw.R`.
7. Use a separate `ZEROFILL_CACHE` directory for each region's zerofilled data (the cache is keyed only by species name, not by polygon).

The main bottleneck for larger regions will be RAM and disk: the zerofilled cache scales linearly with species count × checklist count.
