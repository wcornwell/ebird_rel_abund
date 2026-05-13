# Configuration Guide

The `run_batch_nsw.R` script now uses a YAML configuration file instead of hard-coded paths, making it easy to run the same pipeline for different regions.

## Quick Start

**Run the default NSW pipeline:**
```bash
Rscript run_batch_nsw.R
```
This uses `config.yaml` (NSW + 100 km buffer).

**Run for a different region:**
```bash
export EBIRD_CONFIG=config_victoria.yaml
Rscript run_batch_nsw.R
```

Or in a single command:
```bash
EBIRD_CONFIG=config_victoria.yaml Rscript run_batch_nsw.R
```

## Configuration File Structure

Each YAML file specifies:

| Field | Purpose | Example |
|-------|---------|---------|
| `region` | Region name (for logging) | `NSW`, `Victoria` |
| `description` | Human-readable description | `NSW + 100 km buffer (ACT, VIC, QLD, SA for training data)` |
| `raw_data` | Directory containing EBD/sampling files | `ebirdabund/raw_data` |
| `ebd_files` | List of EBD observation files (relative to `raw_data`) | `[ebd_AU-NSW_*.txt, ebd_AU-ACT_*.txt, ...]` |
| `sampling_files` | List of sampling files (must match `ebd_files` 1:1) | `[ebd_AU-NSW_*_sampling.txt, ...]` |
| `species_list_file` | Species list CSV with reporting rates | `nsw_species_list.csv` |
| `taxonomy_file` | Common → scientific name lookup | `nsw_ebird_taxonomy.csv` |
| `covariate_cache` | Shared cache (GADM, climate tiles) | `ebirdabund_cache` |
| `zerofill_cache` | Region-specific zero-filled cache | `ebirdabund_cache_nsw_buffer` |
| `output_dir` | Output directory for TIFs, PNGs, models | `species_maps` |
| `log_file` | Run log CSV | `batch_nsw_log.csv` |
| `botw_path` | BirdLife range polygons (global) | `botw_species/BOTW_2025.gpkg` |
| `grid_resolutions_km` | Prediction grid resolutions | `[3, 9]` |
| `reporting_rate_threshold` | Species inclusion threshold | `0.005` |
| `study_polygon` | Region boundary and buffer settings | See below |

### Study Polygon

Define the boundary and projection used for the study:

```yaml
study_polygon:
  country: AUS          # GADM country code
  region: New South Wales  # GADM region name (NAME_1)
  buffer_metres: 100000    # Buffer distance for training data
  proj_processing: 3577    # Processing projection (GDA94/Albers)
  proj_output: 4326        # Output projection (WGS84)
```

## Creating a New Region Config

1. Copy `config.yaml` to `config_yourregion.yaml`:
   ```bash
   cp config.yaml config_yourregion.yaml
   ```

2. Update the fields:
   - **EBD/sampling files**: Point to the new region's files
   - **Species list & taxonomy**: Generate these for the new region (see `nsw_species_list.R`)
   - **Cache directories**: Use region-specific names so they don't mix
   - **Output directories**: Region-specific to avoid overwriting
   - **Study polygon**: Region name and buffer settings

3. Run with the new config:
   ```bash
   EBIRD_CONFIG=config_yourregion.yaml Rscript run_batch_nsw.R
   ```

## Example: Victoria

A template example is provided in `config_example_victoria.yaml`. To use it:

1. Ensure you have Victoria EBD/sampling files in `ebirdabund/raw_data/`
2. Generate `vic_species_list.csv` and `vic_ebird_taxonomy.csv` (or copy from NSW as a starting point)
3. Run:
   ```bash
   EBIRD_CONFIG=config_example_victoria.yaml Rscript run_batch_nsw.R
   ```

## Cache and Output Organization

With per-region configs, it's recommended to organize outputs by region:

```
species_maps/          # NSW outputs
species_maps_vic/      # Victoria outputs
species_maps_qld/      # Queensland outputs

ebirdabund_cache_nsw_buffer/     # NSW zero-filled cache (region-specific)
ebirdabund_cache_vic_buffer/     # Victoria zero-filled cache
ebirdabund_cache_qld_buffer/     # Queensland zero-filled cache

ebirdabund_cache/      # Shared covariate cache (GADM, climate, etc.)
```

## Notes

- The `covariate_cache` is shared across all regions (contains global GADM, climate tiles, etc.).
- The `zerofill_cache` should be region-specific to avoid mixing checklist pools.
- EBD/sampling file paths are relative to `raw_data/` — update them for new regions.
- The `EBIRD_CONFIG` environment variable defaults to `config.yaml` if not set.
