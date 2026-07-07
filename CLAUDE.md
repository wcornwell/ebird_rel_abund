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
├── build_observer_expertise.R   # Stage 1 — observer-expertise BLUPs (run once)
├── config.yaml              # Region paths, polygon, thresholds (see CONFIG.md)
│
├── analysis/                # Exploratory, evaluation, and one-off scripts
│   ├── reference/                     # External reference material
│   │   ├── 06_abundance.Rmd           # Cornell eBird Best Practices tutorial
│   │   └── workflow_comparison.md     # Comparison of tutorial vs ebirdabund package
│   ├── tree_height/                   # Canopy-height layer investigations (CHMv2, ETH-Lang, GLAD)
│   ├── palsar/                        # PALSAR backscatter as a candidate covariate
│   ├── observer/                      # One-species sanity check for the observer term
│   ├── count/                         # Visual checks / threshold tuning for observation counts
│   ├── modeling/                      # Family + transform comparisons (NB vs Tweedie, log vs linear, pop_density F-stats)
│   └── integration/                   # Multi-state buffered run + cross-region plot comparisons
│
├── ebirdabund_cache/        # Auto-generated cache (do not commit)
│   ├── cov_stack_v7_*.tif   # Covariate raster stack (bbox + version keyed;
│   │                        #   v5 = Meta CHMv2 replaces OzTreeMap for tree_height;
│   │                        #   v6 = adds nightlights (Falchi/Cinzano 2015);
│   │                        #   v7 = adds palsar_hv (ALOS-2 PALSAR-2 HV gamma0))
│   ├── meta_chmv2_*.tif     # Cached Meta CHMv2 canopy-height raster (one per bbox)
│   ├── palsar_hv_*.tif      # Cached PALSAR HV gamma0 dB raster (one per year + bbox)
│   ├── nightlights/         # Falchi/Cinzano World Atlas 2015 GeoTIFF (local input)
│   └── gadm/, climate/, elevation/, landuse/, population/, soil/  # Downloaded tiles
│
├── ebirdabund_cache_nsw_buffer/  # Region-specific cache (do not commit)
│   ├── zerofilled_*.rds         # Per-species, NSW+buffer pool
│   ├── sampling_master.rds      # Shared sampling pool + covariates + expertise
│   ├── x_count_ids.rds          # Checklists with any X-count entry
│   └── observer_expertise.rds   # Stage 1 BLUPs (run build_observer_expertise.R)
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
├── nsw_species_list.csv     # 551 species, reporting_rate >= 0.5%
├── nsw_ebird_taxonomy.csv   # common_name → scientific_name for BOTW lookup
├── batch_nsw_log.csv        # Per-species run log — columns: common_name,
│                            #   scientific_name, status, exclusion_reason,
│                            #   run_date, n_checklists, n_positive, dev_expl,
│                            #   model_sum, peak_doy, peak_time, max_obs_count,
│                            #   max_modeled_abd, range_source, spatial_cv,
│                            #   error_message, commit_sha
│                            # commit_sha is the short git hash of HEAD when the
│                            # row was written; -dirty suffix means EBIRD_ALLOW_DIRTY=1
│                            # was used to bypass the clean-tree gate.
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

2. build_observer_expertise.R  (Stage 1 calibration — run once after pre-cache,
   before run_batch_nsw.R; ~overnight on the full buffered region)
   → reads all 5 sampling + EBD files, computes species count per checklist
   → filters to observers with >= 20 complete checklists (drops rare contributors
     whose BLUPs would be heavily shrunk anyway)
   → fits a single global NB GAM: n_species ~ effort smooths + protocol_type +
     s(observer_id, bs="re"), via mgcv::bam
   → extracts per-observer random-effect BLUPs as the expertise score
   → writes ebirdabund_cache_nsw_buffer/observer_expertise.rds
   After this completes, delete ebirdabund_cache_nsw_buffer/sampling_master.rds
   so it rebuilds with observer_expertise attached on the next batch run.

3. run_batch_nsw.R  (main run)
   a. GIT GATE: aborts if any tracked file has uncommitted changes (override
      with EBIRD_ALLOW_DIRTY=1). Captures `git rev-parse --short HEAD` as
      commit_sha and stamps it on every row of batch_nsw_log.csv so a result
      can be traced back to the exact source state that produced it.
   b. POLYGON: NSW state boundary + 100 km buffer (GDA94/Albers, EPSG:3577)
      — includes all of ACT and border regions of VIC/QLD/SA as training data
   c. PRE-CACHE: reads all 5 state sampling files filtered to buffered bbox,
      then scans all 5 EBD files once. From that scan it (i) collects the
      union of checklist IDs with any X-count entry → x_count_ids.rds,
      (ii) drops those checklists from the shared sampling pool, then
      (iii) writes zerofilled_{species}.rds per species in
      ebirdabund_cache_nsw_buffer/ (every species now shares an identical
      denominator). Prevents parallel workers competing for files.
   d. EBIRDST RANGES: for each species in the ebirdst 2023 catalogue, downloads
      range .gpkg files (smooth + raw, 27 km + 9 km) to the local ebirdst data
      directory (~11 s/species, ~48 min first run). Skips already-downloaded
      species, so subsequent runs complete instantly. Requires EBIRDST_KEY env var
      (set automatically in the script from the stored key).
   e. COVARIATES: download ESA/SRTM/WorldClim/WorldPop/JRC/SoilGrids and stream
      Meta CHMv2 canopy height once → cov_stack_v7_{bbox}.tif in
      ebirdabund_cache/ (bbox + version keyed, shared across runs). v7 is the
      current key (v5 = CHMv2 tree_height, v6 = nightlights, v7 = palsar_hv;
      see the version history in Repo Layout).
   f. SAMPLING MASTER: one-time covariate extraction + observer_expertise join
      → sampling_master.rds. Built on first call to prepare_sampling_master().
   g. BATCH: estimate_abundance_batch() — parallel GAM per species
      For each species:
        load_ebird() [uses cache + expertise] → subsample_hex()
        → fit_gam() → predict_abundance() [two resolutions] → save .tif + .png
        Map PNGs include NSW state border overlay for reference.
   h. STACK: terra::rast() all .tif → nsw_abundance_stack_{res}.tif
   i. LOG: append to batch_nsw_log.csv after each chunk (every row carries
      the commit_sha captured in step a).

4. rez_abundance.R  (post-processing, run after batch)
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

**GAM family**: Negative binomial (`mgcv::nb()`), chosen a priori via `analysis/modeling/test_tweedie_vs_nb.R`.

**Effort covariates** (control for observation bias):
- `day_of_year` — cyclic cubic spline, k=10, knots c(0, 365)
- `time_observations_started` — cyclic cubic spline, k=10, knots c(0, 24)
- `log(duration_minutes)` — k=4
- `log1p(effort_distance_km)` — k=4 (log1p because stationary counts = 0)
- `number_observers` — k=4
- `protocol_type` — parametric factor
- `observer_expertise` — k=5, added only when `observer_expertise.rds` exists
  in the cache (see Observer Expertise below). Predicted at the median score
  across training data → values represent the population-average observer.

**Habitat covariates** (abundance signal):
- Land cover fractions from ESA WorldCover: `lc_trees`, `lc_grassland`, `lc_cropland`, `lc_built` — k=4 (excludes `lc_shrubs`, collinear)
- `elevation` — k=6
- `log(precip_annual)` — k=6 (log: ~1 order of magnitude range in NSW)
- `temp_annual` — k=6
- `log1p(pop_density)` — k=4 (log1p: spans ~4 orders of magnitude)
- `water_occ` — k=4
- `clay` — k=4 (SoilGrids 0–5 cm clay fraction; NA at ocean/water set to 0)
- `log1p(nightlights)` — k=4 (Falchi/Cinzano World Atlas of Artificial Night Sky Brightness 2015, ~1 km global GeoTIFF, read from a local copy at `ebirdabund_cache/nightlights/World_Atlas_2015.tif`. Float32 radiance proxy in mcd/m²; NoData values over oceans and the polar gap are set to 0 at extraction. log1p because the distribution is heavily right-skewed with most of NSW near 0; chosen as a complement to `pop_density` because population captures *where people live* while nightlights captures *where energy is used at night* — they are correlated but not identical, with major roads, industrial sites, regional towns, and ports lit beyond their resident populations.)
- `log1p(tree_height)` — k=4 (Meta/WRI Canopy Height Maps v2 — DINOv3, 2024 imagery; native ~1 m, streamed via VSICURL from `s3://dataforgood-fb-data/forests/v2/global/dinov3_global_chm_v2_ml3/chm/` at OVERVIEW_LEVEL=4 ≈ 38 m, the closest pre-built overview to our ~30 m cache target. Pixels >60 m are set to NA (model misclassifies tall narrow non-vegetation — wind turbines on ridge crests, transmission towers, silos — as canopy). The cache is then **p90-aggregated** to the master ~1 km template: each output cell takes the 90th percentile of its source pixels, capturing emergent canopy structure (tall trees as hollow/perch habitat) while excluding residual artefacts that sit in the top ~0.1%. NA values at extraction set to 0. Replaces the earlier OzTreeMap layer, which had WCS grid artefacts.)
- `palsar_hv` — k=4, no transform (ALOS-2 PALSAR-2 yearly mosaic gamma0, HV polarization, year 2020 by default; L-band ~24 cm wavelength penetrates canopy and responds to woody biomass/volume scattering, so it picks up sub-canopy structure that CHMv2 height cannot resolve. Streamed via VSICURL from Microsoft Planetary Computer's COG mirror of the JAXA ALOS-2 PALSAR-2 collection at `OVERVIEW_LEVEL=3` ≈ 200 m, then aggregated to the master ~1 km template **in linear power space** before converting back to dB — dB is a log scale, so arithmetic-mean of dB would underweight bright pixels. STAC search → SAS-signed per-asset URLs → VSICURL; no MPC subscription key required. JAXA calibration: gamma0(dB) = 20·log10(DN) − 83; DN==0 is no-data. NA at extraction is median-imputed because `water_occ` already encodes the wet/dry signal and ocean-pixel NAs shouldn't pull predictions to an extreme. No log transform applied because dB is already on a log scale. Adopted from `covariate_tests/palsar_hv_2020/`: forest-interior species are the biggest winners — Superb Lyrebird +0.0039 dev expl, Eastern Whipbird +0.0020, Rose Robin +0.0027; Pacific Black Duck and Welcome Swallow correctly show no gain.)

Effort covariates are log-transformed because they are right-skewed and the log scale distributes knots more evenly. The transformation happens inside the GAM formula string so raw column values flow through unchanged — the prediction surface (`duration_minutes = 60`, `effort_distance_km = 1`) is evaluated correctly.

**Fitting**: `mgcv::bam()` with `fREML`, `discrete=TRUE`, `select=TRUE`, and a BIC-like complexity penalty `gamma = log(nrow(df)) / 2` (≈ 6.4 at the current NSW+buffer scope). Falls back to `select=FALSE` then `discrete=FALSE` on convergence failure. The `gamma` choice was validated against the conventional `gamma=1.4` in a 19-species paired-CV sweep (`gamma_sweep/full/`): at BIC the partial-effect smooths drop from edf ~70 to ~50–65 and the wiggles that lacked biological justification flatten out. Rank-skill metrics (Spearman, holdout deviance explained) degrade by ~0.001 on most species, while RMSE/MAE improve on 14–16/19 species — particularly for low-detection species where over-flexible tails were the main failure mode. Smoother, more interpretable fits were judged worth the marginal predictive cost.

**Prediction**: Standardised effort (60 min, 1 km, 1 observer, Traveling Count) at the **model-based peak DOY and time-of-day** for the species — i.e. the conditions under which detection is highest, following the Cornell eBird Best Practices approach. Peak DOY is the argmax of the link-scale **lower 95% CI** (`fit - 1.96 * se`) from a daily sweep over `1:365` at mean habitat and standard effort; peak time is then the argmax of the same lower-CI sweep over `seq(0, 24, length.out = 300)` at the chosen peak DOY. Using the lower CI rather than the mean keeps the peak inside the well-supported region of the smooth — without it, sparsely-sampled times (e.g. pre-dawn) can win the argmax on point estimate alone. Predictions therefore represent "expected count in a 1-hour, 1-km, single-observer Traveling Count at the species' peak season and time of day", which is systematically higher than an annual/diurnal average. When the model includes the observer-expertise smooth, prediction is at the **median** training-data expertise (population-average observer), not the top of the distribution. Abundance is capped at the 90th percentile of non-zero training counts (on the log scale, before exponentiation) to prevent extrapolation overflow.

**Exclusion criteria**: fewer than 50 positive checklists after hex subsampling, or model sum ≤ 1e-5 across the polygon.

**Range masking**: ebirdst 2023 preferred; falls back to BOTW 2025 if the species is not in the ebirdst catalogue or if the ebirdst mask produces zero abundance over the study region. In the current NSW run (394 modelled species), 218 were masked with ebirdst 2023 ranges, 171 with BOTW 2025, and 41 left unmasked. Ranges are pre-downloaded by `run_batch_nsw.R` before the batch (download step is skipped for species already cached locally). The `range_source` column in `batch_nsw_log.csv` records which source was used per species.

**BOTW name aliases** (`botw_name_aliases.csv` at repo root): eBird/Clements and BirdLife/HBW disagree on a number of recent genus reshuffles, gender-agreement renames, and splits (e.g. Plum-headed Finch is `Emblema modestum` in eBird but `Neochmia modesta` in BOTW). Without a translation layer, the BOTW fallback returns 0 rows for these species and the prediction is left unmasked. `load_range_botw()` consults the alias CSV (passed in by `run_batch_nsw.R`) and rewrites its query to the BOTW name on a hit. Columns: `ebird_sci_name`, `botw_sci_name`, `common_name`, `note`. An empty `botw_sci_name` records "no BOTW equivalent exists" and short-circuits the lookup to avoid a guaranteed-zero round-trip. The file is region-agnostic — add entries as new taxonomy splits surface (audit by filtering `batch_nsw_log.csv` for `range_source == "unmasked"` after each batch and checking BOTW for synonyms).

Aliases only fix the **name-mismatch** failure mode. A separate, larger group of species was previously left `range_source == "unmasked"` (or truncated mid-range) because `load_range_botw`'s old conservative filter (`seasonal IN (1,2,3)`, i.e. resident/breeding/non-breeding only) excluded their qualifying polygons — passage migrants (jaegers, sand-plovers, Wandering Tattler, Aleutian Tern), pelagic seabirds (albatrosses, shearwaters), and nomads whose coastal/eastern occurrence sits in a `seasonal=4` (passage) or `seasonal=5` (occurrence uncertain) polygon. **Red-necked Avocet** was the clearest case: its resident polygon stops at 147°E, so the old filter clipped the predicted abundance along a hard vertical line down the middle of NSW, masking out the whole coast where it actually occurs. The filter now keeps **`seasonal IN (1,2,3,4,5)`** (all season codes; `presence`/`origin` filtering is unchanged), which recovers these ranges. Re-run the batch to refresh the affected TIFs (the change is in code, not in already-written outputs).

---

## Observer Expertise (Stage 1)

eBird checklists are submitted by observers of varying skill. Without controlling for this, per-observer variance leaks into the effort and habitat smooths, biasing predictions toward what frequently-active (often more skilled) observers report. Following **Kelling et al. (2015, *PLoS ONE*)** and **Johnston et al. (2021, *Diversity & Distributions*)**, we partition observer skill out of the per-species models via a two-stage approach.

**Stage 1 — global calibration** (`build_observer_expertise.R`, run once):

1. Read all 5 sampling files (NSW + ACT + VIC + QLD + SA) and clip to the buffered polygon.
2. Scan all 5 EBD files once and compute `n_species` = number of distinct species reported per complete checklist (all rows, including X-counts — they still indicate detection).
3. Apply the same effort/time filters as `clean_ebird()` (duration 10–300 min, distance ≤ 10 km, ≤ 10 observers).
4. Drop observers with fewer than 20 qualifying complete checklists (`MIN_CHECKLISTS_PER_OBSERVER`); their BLUPs would be heavily shrunk anyway and the level cardinality makes `bs="re"` expensive.
5. Apply the same hex × week subsampling used per-species (`spacing_km = 5`).
6. Fit a single global negative-binomial GAM:

   ```r
   n_species ~ s(log(duration_minutes), k = 5) +
               s(log1p(effort_distance_km), k = 5) +
               s(time_observations_started, bs = "cc", k = 7) +
               s(day_of_year, bs = "cc", k = 7) +
               protocol_type +
               s(observer_id, bs = "re")
   ```

7. Extract the observer random-effect BLUPs and save as `ebirdabund_cache_nsw_buffer/observer_expertise.rds` (a small data frame: `observer_id`, `expertise`, plus attributes recording the fit metadata).

The fit is expensive — on a single Superb Fairywren species (387k × 13k observers) `bs="re"` took ~15 h. The Stage 1 fit operates on a similar row count but with more observer levels (~20k after the activity filter), so expect overnight runtime. Re-run only when the EBD or polygon changes.

**Stage 2 — use in per-species models** (automatic):

`attach_observer_expertise()` in `load_ebird.R` joins the score onto each species' training data via `observer_id`. Observers not in the calibration set (rare contributors, or first-time appearances) are assigned the median expertise so they contribute at the population-average level. `observer_id` is dropped after the join. `build_gam_formula()` then adds `s(observer_expertise, k=5)` to the per-species formula. `predict_abundance()` sets `observer_expertise` to the median across training data so the standardised-effort prediction is for an average-skill observer (cf. eBird Status & Trends, which predicts at the *top* of the skill distribution — an "expert eBirder").

**One-species sanity check before Stage 1 was operationalised** (`analysis/observer/test_observer_re.R`): for Superb Fairywren, adding `s(observer_id, bs="re")` directly to the per-species formula dropped the per-cell mean predicted abundance from 1.262 to 0.898 (ratio 0.71). That is, ~30 % of the "high numbers" feeling was observer-skill variance bleeding into the other smooths. The expertise-score approach in Stage 1 + Stage 2 is the citable, scalable form of the same correction.

**If `observer_expertise.rds` is absent**, the pipeline runs as before — no observer term, no expertise join — so Stage 1 is a strict add-on.

**Cache rebuild order after Stage 1:** the sampling master is built once and reused for every species; if it predates Stage 1 it won't carry the expertise column. Delete `ebirdabund_cache_nsw_buffer/sampling_master.rds` after running `build_observer_expertise.R` so the next batch rebuilds the master with expertise attached.

---

## Evaluation

Use `evaluate_model_cv(fit$data, k=5)` for held-out metrics (Spearman ρ, Pearson r, RMSE, holdout deviance explained). The evaluation scripts follow the pattern in `analysis/modeling/test_tweedie_vs_nb.R`:

1. Call `fit_species_model()` to get the subsampled data frame.
2. Build alternative formulas manually if testing a variant.
3. Pass each formula to `evaluate_model_cv()`.
4. Compare the `$summary` vectors.

`compare_covariate_models(fit$data)` runs effort-only vs. full-habitat in one call and prints a table.

---

## Covariate Variant Testing

For evaluating a candidate covariate change — swapping one layer for another (e.g. CHMv2 → ETH-Lang for `tree_height`) or adding a new layer — without invalidating the main cache. Implemented in `ebirdabund/R/test_covariate.R` (`test_covariate_variant()`). Builds on `evaluate_model_cv()` / `compare_covariate_models()`; does not replace them.

**Design constraints:**
- Reuses `sampling_master.rds` and `zerofilled_*.rds` as-is — never rewrites them. The variant column is joined in-memory at the start of each test run.
- Variant covariate values are extracted once from the alternative raster (or extraction function) and cached separately at `ebirdabund_cache_nsw_buffer/variant_extracts/{variant_name}.rds`, keyed by checklist ID. Repeat tests are instant.
- Paired k-fold CV per species: identical fold assignments for baseline and variant (`evaluate_model_cv()` takes an optional `fold_ids` argument for this) → matched-pair comparison with much higher statistical power than independent CVs.

**API:**

```r
test_covariate_variant(
  variant_name   = "tree_height_eth_lang",
  variant_source = "path/to/eth_global_2020.tif",   # or function(bbox, template) -> SpatRaster
  change         = list(swap = "tree_height"),       # or list(add = "tree_height_eth")
  species        = NULL,                              # default: test set CSV
  k              = 5
)
```

**Test species set** lives at `analysis/covariate_test_species.csv` — a fixed list hand-picked to span habitat guilds (forest, open-country, water, urban, arid, alpine). Same set on every variant test so deltas are directly comparable across variants over time. Sized to complete in ~10–20 min.

**Outputs** under `covariate_tests/{variant_name}/`:
- `per_species_metrics.csv` — species × {baseline, variant, Δ} for Spearman ρ, Pearson r, MAE, RMSE, holdout deviance explained.
- `summary.txt` — paired Wilcoxon across species, win/tie/loss counts per metric.
- `smooth_comparison.pdf` — baseline vs variant partial-effect plot per species for the swapped/added covariate.
- `metrics_full.rds` — raw CV objects + full-data fits for re-analysis.

First concrete use: `analysis/tree_height/test_tree_height_alternatives.R` runs the CHMv2 two-stage aggregation variant and a GLAD/Potapov 2019 variant against the current CHMv2 baseline, then merges per-species deltas into a three-way comparison CSV.

---

## Possible Improvements (flagged during review)

### Medium priority

- **CHMv2 stripe artefacts**: The Meta/WRI CHMv2 raster has visible near-vertical stripes inherited from Sentinel-2 orbit/swath boundaries (different image counts, sun angles, and seasonality between adjacent orbital tracks bake discontinuities into DINOv3's height estimate). Visible in `tree_height_chmv2.png`. Same issue is reported in Moudrý et al. 2024 (Ecosphere) for Lang 2023 and Tolan 2024. Three fix paths, in order of effort:
    - **(C, easiest, ~10 LOC)** Two-stage aggregation in `load_meta_chmv2` (covariates.R): 38 m → ~250 m **mean** → 1 km **p90**. Current single-step p90 over a 26×26 window preferentially picks the lit side of stripes; an intermediate mean smooths within-stripe before the p90 captures emergent canopy. Also try anisotropic pre-smoothing (`terra::focal()` 1×5 along the stripe direction) and `OVERVIEW_LEVEL=5` (~76 m) which is already averaged from native 1 m.
    - **(A, medium, ~1 day)** Post-hoc destriping of the cached `meta_chmv2_*.tif` before the p90 step, via combined wavelet–FFT filtering (Münch et al. 2009; reference impl: [DHI-GRAS/rmstripes](https://github.com/DHI-GRAS/rmstripes), Python — call via `reticulate` or port to R with `wavethresh` + `fft`). FFT isolates the periodic vertical-frequency band; wavelet damps only the affected detail bands so real edges are preserved.
    - **(B, biggest payoff for geographic scaling)** Ensemble with an independent product. Stripes are deterministic in a single product but uncorrelated across products. Median of CHMv2 + ETH-Lang 2020 (Sentinel-2 + GEDI, 10 m, free on GEE / Zenodo) cancels them. Optionally blend in GLAD/Potapov 2021 (30 m, Landsat-based, no S2 stripes) or GEDI L4B (1 km) as a long-wavelength anchor for cells where products disagree by >X m.
    - Recommended sequencing: do C first (cheap, may be sufficient), then A on the cached raster if stripes persist, then B if also rolling out to other regions.
- **Generalise `nsw_species_list.R`**: It's NSW-specific (variable names, output file name, plot title). Wrapping it into a function `generate_species_list(region_name, ebd_path, samp_path, ...)` would make it reusable for geographic scaling.
- **Taxonomy file**: `nsw_ebird_taxonomy.csv` is derived from the NSW EBD. For other regions, this needs to be generated from the regional EBD or replaced with the full eBird taxonomy (available from eBird as `eBird_taxonomy_v*.csv`).
- **Covariate cache versioning**: The cache key uses `v7` as a hard-coded version tag. If covariate layers are updated, this needs a manual bump. A hash of the source URLs or a date stamp would be more robust.
- **REZ script parameterisation**: `rez_abundance.R` paths and the `TOP_N = 50` cutoff are hard-coded. Making it accept command-line arguments would make it easier to run for new regions or different analysis polygons.
- **Parallel memory**: `n_cores` is read from `config.yaml` (default `6`). Each worker peaks at ~2 GB RSS (sampling master + zerofill + mgcv::bam state). Three mitigations are now built into `estimate_abundance_batch()`:
    - `terra::terraOptions(memfrac = 0.05)` is set on every worker, so the per-process terra raster cache is capped at ~5% of system RAM. The default of ~60% is per-process and so would, with N workers, claim N×60% — slow drift to jetsam pressure on long runs.
    - `run_species()` does `rm(model_fit, pred); gc(full = TRUE)` immediately before returning, so the bam fit (often 200–500 MB) is released and a generational GC runs before the worker picks up its next species.
    - The worker lambda in the `parLapplyLB` call does the same after `run_species` returns, catching anything that leaked via the error path.
    
    Originally 9 workers on a 32 GB M-series Mac drove the macOS compressor past 13 GB twice in a 4 h 42 m run, killing two chunks with silent worker deaths surfacing as `"error reading from connection"`. Six workers + the cleanup above ran the same workload at pressure level 1 (Normal). For larger regions or higher resolutions consider a `max_ram_gb` guard that caps workers from total RAM, or periodic worker recycling (rebuild the cluster every N species) to fully reset the malloc heap.

### Lower priority

- **Incremental stack building**: The stack-building step at the end of `run_batch_nsw.R` re-reads all TIFs. For very large species lists, this is slow. Could build incrementally or defer to a separate script.
- **Smooth plotting**: `plot_gam_smooths()` in `utils.R` saves a static PNG. Interactive HTML (e.g. via `plotly`) would make diagnosing smooth shapes much faster.
- **GEDI canopy structure (FHD / PAVD) — deferred, conditional on geographic scaling**: Recent bird SDM literature (Burns 2020 ERL; Schulte to Bühne et al. 2024 RSE) shows canopy structure metrics *beyond* height — foliage height diversity, plant area volume density at vertical strata — add signal for forest-interior species. The [Burns/Hakkenberg/Goetz 2024 gridded GEDI product](https://www.nature.com/articles/s41597-024-03668-4) ([GEE: `LARSE/GEDI/GEDI_GRIDDEDVEG_002_COUNTS_V1_1KM`](https://developers.google.com/earth-engine/datasets/catalog/LARSE_GEDI_GRIDDEDVEG_002_COUNTS_V1_1KM)) packages 36 such metrics at 1 km globally. Two caveats make this low priority for NSW:
    - **Coverage is shot-sample, not wall-to-wall.** GEDI is a sampling lidar on the ISS (8 tracks, ~600 m apart, 25 m footprints). Gridded products recommend ≥10 shots per 1 km cell for reliability; cells below threshold are NA. In NSW the arid/saltbush west and Riverina have low shot density over low/bare cover, so a large fraction of the study polygon would be NA. The batch pipeline expects no-NA covariates, so we'd need a hybrid imputation (zero-prior in sparse-canopy cells, bbox median otherwise) plus a binary "GEDI-observed" mask as a separate predictor — non-trivial engineering.
    - **Expected gain is concentrated in a minority of species.** Forest-interior species along the GDR/coast (treecreepers, sittellas, scrubwrens, lyrebird, large hollow-nesters) are where vertical structure beyond canopy height matters; for these, GEDI shot density is adequate. For arid/open-country species (most of NSW by area) the layer would mostly be imputed and uninformative, and CHMv2 already covers the height gradient.
    - **Reconsider when**: scaling to Tasmania, wet-tropics QLD, or NZ — regions with a higher proportion of forest-interior species and denser GEDI coverage. At that point the cost/benefit flips. Scoping notes in chat history (2026-05-18).

---

## Geographic Scaling Plan

The package (`ebirdabund/`) is region-agnostic — `fit_species_model()` and `estimate_abundance_batch()` accept any `sf` polygon and any number of EBD/sampling files as character vectors.

**Current setup:** NSW + 100 km buffer as the study polygon, with EBD training data from NSW, ACT, VIC, QLD, and SA. The buffer captures ~334k checklists per species after hex subsampling (spacing 5 km).

The work needed to add a new region:

1. Download regional EBD + sampling files.
2. Run a per-region adaptation of `nsw_species_list.R` (still needs generalising — see Improvements).
3. Generate a taxonomy CSV for the region (or use the full eBird taxonomy).
4. Copy `config.yaml` to `config_{region}.yaml`, edit paths/polygon/species-list, and run `EBIRD_CONFIG=config_{region}.yaml Rscript run_batch_nsw.R`. `config_example_victoria.yaml` is a worked example; `CONFIG.md` documents every key.
5. The covariate downloads (ESA, WorldClim, etc.) are global and will re-use cached tiles if the bounding box overlaps.
6. Range masking is already global: BOTW covers all species, and ebirdst covers whichever species it models. The pre-download step in `run_batch_nsw.R` runs automatically for the new species list before the batch — `EBIRDST_KEY` is set in the script.
7. Use a separate `zerofill_cache` directory per region in the config (the cache is keyed only by species name, not by polygon, so different regions' caches must not share a directory).

The main bottleneck for larger regions will be RAM and disk: the zerofilled cache scales linearly with species count × checklist count.

---

## REZ Risk-Assessment Layer (seasonality, taxonomy, traits)

Downstream products built on top of the abundance stacks, for a bird collision/
displacement **risk assessment of NSW Renewable Energy Zones (REZs)**, focused on
onshore wind. Three deliverables plus a documented data gap.

### 1. Per-REZ seasonality (`rez_seasonality_batch.R`, `ebirdabund/R/seasonality.R`)

Answers "how seasonal is each taxon in each REZ, and if so what's the window?"
Distinct from the production model's single global `s(day_of_year)` (which is
spatially constant). Per species, **one** shared NB GAM reuses the production
effort+habitat structure and adds an **ordered-factor by-REZ cyclic seasonal
smooth**: reference level `statewide` (rest-of-state) carries the global
`s(day_of_year)`; each target REZ gets a *difference* smooth, so its curve =
global + deviation, and sparse REZs shrink toward statewide (`select=TRUE`).
Regions: `statewide` + 3 target REZs (`central_west` = Central-West Orana,
`new_england`, `south_west` = Hay). Metrics via **posterior simulation** from
`vcov(fit)`: `seasonality_index` = 1 − trough/peak ∈ [0,1]; `is_seasonal` =
(≥20 positive detections) AND P(peak ≥ 2× trough) ≥ 0.95; window = contiguous
DOY where the standardised curve ≥ its annual mean (cyclic-aware).

- **Output**: `rez_seasonality/seasonality_all.csv` (one row per species × region;
  551 × 4 = 2204 rows) + `seasonality_all_metadata.md` (full data dictionary +
  sources + methods — keep it in sync when columns change).
- **Fit chain is discrete-only**: the non-discrete `bam` fallback took 30+ min
  per species on ~334k rows and caused apparent "hangs" — it was removed. Chain:
  full discrete BIC → full discrete γ=1.4 → discrete no-select → simplified
  formula → fail fast. **Do not re-add a non-discrete fallback.**
- **Workers**: run at **≤4** (heavier than the abundance fit; 5 hit jetsam
  deaths). `subsample_hex` draws randomly, so a borderline species can fail on
  one draw and succeed on a re-run (the script resumes on species not yet in the
  output CSV).
- **Caveat**: the model controls for observation *effort* but not for a species'
  own seasonal *detectability* (breeding-season song/display). Cross-check
  `is_seasonal` against the migration columns; residents can read as seasonal.

### 2. Taxonomy / WLAB crosswalk (`ebirdabund/R/taxonomy.R`)

`build_taxonomy_crosswalk()` maps modelled species to the Reid/Baker NSW bird
working list (WLAB, `taxonomy/reidbaker_bird_risk_NSW_species.csv`). Join order:
binomial → common name → `taxonomy/reidbaker_name_aliases.csv` (genus reshuffles,
spelling variants, eBird-split-of-WLAB-subspecies). Output
`taxonomy/species_wlab_crosswalk.csv`; wired into `rez_abundance.R` and the
seasonality output.

- Column is **`risk_assessed`** (renamed from `risk_listed`) = "present in the
  WLAB working list", **not** "threatened" — WLAB includes common species.
- `wlab_review` / `wlab_review_reason = multiple_nsw_subspecies` flags species
  that lump ≥2 NSW-plausible subspecies (see `DEFAULT_NON_NSW_PATTERNS`); we flag,
  don't split.
- **reidbaker is a 493-taxon SUBSET, not a full NSW list.** It genuinely omits
  many common natives (White-faced Heron, Tree Martin, Brown Falcon, Variegated
  Fairywren, …). No aliasing recovers them — they are not in the file. ~415/551
  matched; the rest are introduced species (correctly excluded), rare pelagics
  (unimportant), or these genuine WLAB gaps. For comprehensive native +
  conservation coverage, use **Garnett** (below), not reidbaker.

### 3. Species traits (`extract_traits.R` → `taxonomy/species_traits.csv`)

Per-species functional + conservation traits from four external datasets, joined
by eBird scientific name with per-source alias files. Large source files live in
`taxonomy/traits/` and are **gitignored** (regenerable by the extract scripts):
- **AVONET** (Tobias 2022, eBird taxonomy) → `migration` (Sedentary/Partial/
  Migratory), `feeding_guild` (Trophic.Niche), `trophic_level`,
  `primary_lifestyle`. Alias file: `taxonomy/avonet_name_aliases.csv` (post-2021
  genus reshuffles). AVONET is **conservative for Australian nomads** (Pallid
  Cuckoo coded Sedentary despite spring migration).
- **EltonTraits 1.0** (Wilman 2014) → `foraging_stratum` (dominant vertical
  stratum, incl. `pelagic`).
- **Garnett Australian Bird Data v1.0** (2015, subspecies-level, aggregated to
  species) → `epbc_status`, `nsw_status` (most-threatened across the species'
  taxa), `garnett_movement` (local dispersal / partial / total migrant / nomadic
  / irruptive — catches nomadism AVONET flattens). **Comprehensive Australian
  coverage** (2056 taxa) — fills the reidbaker gaps.
- **BIRDBASE** (Şekercioğlu 2025, 2024 taxonomy) → `iucn_status`,
  `habitat_breadth`.

`extract_migratory.R` writes the intermediate `taxonomy/migratory_status.csv`
(BOTW seasonal-code proxy; superseded by AVONET/Garnett but kept as an audit
column `botw_seasonal_codes`). Traits are joined into `seasonality_all.csv` via
`seasonality_add_traits.R` (idempotent) and into the abundance CSVs by
`rez_abundance.R`.

### 4. The flight-height gap (do not re-chase)

Wind collision risk = exposure (abundance × seasonal presence, which we have) ×
vulnerability. The dominant vulnerability term is **% time flying at rotor-blade
height** — and **no comprehensive species-level flight-height dataset exists**
for terrestrial Australian birds. What exists: good seabird distributions
(Johnston 2014 / BTO-Cook, UK marine only); fragmented per-species GPS tracks
(Movebank, raptors/soaring birds); radar (can't assign to species). Conclusion
after searching: treat as a structural gap — use an **ecology-based proxy**
(rotor-height overlap from `foraging_stratum`/`primary_lifestyle`/`feeding_guild`)
+ real tracking only for the few priority soaring/threatened species.

**AVISTEP** (BirdLife Avian Sensitivity Tool for Energy Planning, Australia): the
download (`avistep/Australia/Australia.gdb`, gitignored) is a **spatial 5 km
sensitivity map** (`avistep_aus_onshore`, category 1–4), **not** per-species
scores — those are internal to AVISTEP and not shipped. Rejected as too
black-box for per-species work, but usable as a REZ-level benchmark: all NSW REZs
rate High–Very High (New England highest ≈3.84, South West lowest ≈3.24; see the
overlay in `rez_plots/avistep_rez_sensitivity.csv`).

**Env note**: R was upgraded 4.5→4.6; run under plain `Rscript` (4.6 has the
needed packages). Do **not** put the 4.5-arm64 library on the path — compiled
packages are ABI-incompatible across the minor bump.
