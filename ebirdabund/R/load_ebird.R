# Expected EBD column positions used by the package. Validated by
# validate_ebd_header() before any fread call relies on integer indices.
EBD_EXPECTED_COLS <- c(
  "6"  = "COMMON NAME",
  "11" = "OBSERVATION COUNT",
  "35" = "SAMPLING EVENT IDENTIFIER"
)

# Validate that an EBD file has the expected column names at the expected
# positions. eBird occasionally adds columns; any drift would silently
# corrupt the integer-column reads in read_ebd_observations() and the
# pre-cache block of run_batch_*.R.
#
# Reads only the first line of the file (cheap, no full scan).
# Returns invisibly TRUE on success; stops with a diagnostic on failure.
validate_ebd_header <- function(ebd_path,
                                expected = EBD_EXPECTED_COLS) {
  con <- file(ebd_path, "r")
  on.exit(close(con), add = TRUE)
  header_line <- readLines(con, n = 1L)
  if (length(header_line) == 0L) {
    stop("EBD file appears empty or unreadable: ", ebd_path)
  }
  cols <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]

  for (idx_chr in names(expected)) {
    i        <- as.integer(idx_chr)
    expected_name <- expected[[idx_chr]]
    actual_name   <- if (i > length(cols)) NA_character_ else cols[i]
    if (is.na(actual_name) || actual_name != expected_name) {
      stop(sprintf(
        paste0(
          "EBD header check failed for %s\n",
          "  Expected col %d = '%s'; got '%s'.\n",
          "  Column layout differs from the reference EBD release.\n",
          "  Update EBD_EXPECTED_COLS in R/load_ebird.R after confirming new positions."
        ),
        ebd_path, i, expected_name,
        if (is.na(actual_name)) "<missing>" else actual_name
      ))
    }
  }
  invisible(TRUE)
}

# Read sampling events from one or more eBird sampling files, filter to
# complete checklists within `bbox`, optionally cap by `date_cutoff`, and
# return a single combined data.frame with renamed columns.
#
# sampling_txt : character vector (length >= 1) of sampling .txt paths.
# bbox         : c(xmin, ymin, xmax, ymax) in WGS84.
# date_cutoff  : optional Date or character "YYYY-MM-DD". Rows with
#                observation_date strictly later than this are dropped.
#
# Strategy: read all columns (avoids fread column-order ambiguity), then
# filter and select with standard R [[ ]] — no data.table NSE needed.
read_sampling <- function(sampling_txt, bbox, date_cutoff = NULL) {
  if (!is.character(sampling_txt) || length(sampling_txt) == 0L) {
    stop("`sampling_txt` must be a non-empty character vector of file paths.")
  }
  if (!is.null(date_cutoff)) {
    date_cutoff <- as.Date(date_cutoff)
    if (is.na(date_cutoff)) {
      stop("`date_cutoff` could not be parsed as a Date.")
    }
  }

  read_one <- function(p) {
    if (!file.exists(p)) stop("Sampling events file not found: ", p)
    df <- as.data.frame(data.table::fread(
      p, sep = "\t", quote = "", showProgress = FALSE, na.strings = ""
    ))

    # Standard R subsetting — no NSE, no backtick issues
    df <- df[df[["ALL SPECIES REPORTED"]] == 1 &
             !is.na(df[["LATITUDE"]]) &
             df[["LATITUDE"]]  >= bbox[2] & df[["LATITUDE"]]  <= bbox[4] &
             df[["LONGITUDE"]] >= bbox[1] & df[["LONGITUDE"]] <= bbox[3], ]

    raw <- c("SAMPLING EVENT IDENTIFIER", "LATITUDE", "LONGITUDE",
             "OBSERVATION DATE", "TIME OBSERVATIONS STARTED",
             "DURATION MINUTES", "EFFORT DISTANCE KM",
             "NUMBER OBSERVERS", "PROTOCOL NAME", "OBSERVER ID")
    df  <- df[, raw, drop = FALSE]
    names(df) <- c("checklist_id", "latitude", "longitude",
                   "observation_date", "time_observations_started",
                   "duration_minutes", "effort_distance_km",
                   "number_observers", "protocol_type", "observer_id")
    df
  }

  message(sprintf(
    "Reading %d sampling event file(s)...", length(sampling_txt)
  ))
  parts <- lapply(sampling_txt, read_one)
  combined <- if (length(parts) == 1L) parts[[1L]] else do.call(rbind, parts)

  if (!is.null(date_cutoff)) {
    obs_date <- as.Date(combined$observation_date)
    keep     <- !is.na(obs_date) & obs_date <= date_cutoff
    n_drop   <- sum(!keep)
    combined <- combined[keep, ]
    if (n_drop > 0L) {
      message(sprintf(
        "  Dropped %d sampling rows after %s (date cutoff)",
        n_drop, format(date_cutoff, "%Y-%m-%d")
      ))
    }
  }

  combined
}

# Read raw observations (common_name, observation_count, checklist_id) from
# one or more EBD .txt files in a single concatenated pass. The column-index
# integer reads (cols 6/11/35) are validated against each file's header
# before fread is called.
#
# ebd_txt              : character vector (length >= 1) of EBD .txt paths.
# species_set          : optional character vector. If non-NULL, only rows
#                        whose `common_name` is in this set are kept.
# valid_checklist_ids  : optional character vector. If non-NULL, only rows
#                        whose `checklist_id` is in this set are kept.
# return_x_count_ids   : if TRUE, the unique set of checklist_ids that had
#                        any observation_count == "X" is captured BEFORE the
#                        species/checklist filters and returned alongside the
#                        observation frame. Used by run_batch_*.R to filter
#                        the shared sampling pool so every species has an
#                        identical denominator.
#
# Returns a data.frame with columns common_name, observation_count, checklist_id,
# or — when return_x_count_ids = TRUE — a list with elements `obs` and
# `x_count_ids`.
read_ebd_observations <- function(ebd_txt,
                                  species_set         = NULL,
                                  valid_checklist_ids = NULL,
                                  return_x_count_ids  = FALSE) {
  if (!is.character(ebd_txt) || length(ebd_txt) == 0L) {
    stop("`ebd_txt` must be a non-empty character vector of file paths.")
  }

  paths <- vapply(ebd_txt, resolve_ebird_path, character(1))

  # Validate headers cheaply (reads one line per file) before the big scan.
  invisible(lapply(paths, validate_ebd_header))

  if (length(paths) == 1L) {
    message("  Scanning EBD: ", basename(paths))
    cmd <- NULL
    src <- paths
  } else {
    # Single-pass concatenation via awk: print header from file 1 only,
    # then data rows from all files. Much faster than 5 separate freads.
    message(sprintf("  Scanning %d EBD files in one pass...", length(paths)))
    cmd <- paste("awk 'FNR>1 || NR==1'", paste(shQuote(paths), collapse = " "))
    src <- NULL
  }

  df <- as.data.frame(data.table::fread(
    if (is.null(cmd)) src else cmd,
    select       = c(6L, 11L, 35L),
    sep          = "\t",
    quote        = "",
    showProgress = TRUE,
    na.strings   = ""
  ))
  names(df) <- c("common_name", "observation_count", "checklist_id")

  # Capture X-count checklist IDs from the unfiltered data; this is the only
  # place they're cheaply available across all species.
  x_count_ids <- if (return_x_count_ids) {
    unique(df[["checklist_id"]][df[["observation_count"]] == "X"])
  } else NULL

  if (!is.null(species_set)) {
    df <- df[df[["common_name"]] %in% species_set, , drop = FALSE]
  }
  if (!is.null(valid_checklist_ids)) {
    df <- df[df[["checklist_id"]] %in% valid_checklist_ids, , drop = FALSE]
  }

  if (return_x_count_ids) list(obs = df, x_count_ids = x_count_ids) else df
}

# Scan one or more EBD files for checklist IDs containing any observation_count
# == "X" entry. Streams via awk + sort -u so memory stays bounded regardless of
# EBD size. Use this when the per-species pre-cache scan is being skipped
# (i.e. all per-species zerofills already exist) but a shared X-count filter
# still needs to be applied.
#
# Returns a character vector of unique checklist IDs.
find_x_count_checklists <- function(ebd_txt) {
  if (!is.character(ebd_txt) || length(ebd_txt) == 0L) {
    stop("`ebd_txt` must be a non-empty character vector of file paths.")
  }
  paths <- vapply(ebd_txt, resolve_ebird_path, character(1))
  invisible(lapply(paths, validate_ebd_header))

  message(sprintf("  Scanning %d EBD file(s) for X-count checklists...",
                  length(paths)))
  cmd <- paste(
    "awk -F'\\t' 'FNR>1 && $11==\"X\"{print $35}'",
    paste(shQuote(paths), collapse = " "),
    "| sort -u"
  )
  res <- system(cmd, intern = TRUE)
  message(sprintf("    Found %d unique X-count checklists.", length(res)))
  res
}

# Single-species wrapper kept for backward compatibility with load_ebird()'s
# fallback path. Returns checklist_id + observation_count only.
read_ebd_species <- function(ebd_txt, species) {
  message("Scanning EBD file for '", species, "' (large file, please wait)...")
  df <- read_ebd_observations(ebd_txt, species_set = species)
  df[, c("checklist_id", "observation_count"), drop = FALSE]
}

# Left-join all complete checklists with species observations.
# Checklists where the species was not detected receive observation_count = 0.
zero_fill <- function(sampling_df, ebd_df) {
  merged <- merge(sampling_df, ebd_df, by = "checklist_id", all.x = TRUE)
  merged$observation_count[is.na(merged$observation_count)] <- "0"
  merged$species_observed <- merged$observation_count != "0"
  merged
}

# Discard ecologically unreliable rows and standardise column types.
#
# The `observation_count != "X"` filter is retained as a safety net for code
# paths that bypass the shared-pool X-count exclusion (e.g. single-species
# slow path in load_ebird() when no x_count_ids.rds is present). When the
# shared pool has already excluded X-count checklists from sampling, this
# filter matches nothing.
clean_ebird <- function(zf, max_count = 200L, min_duration_minutes = 10L) {
  zf |>
    dplyr::filter(
      .data$observation_count != "X",
      .data$duration_minutes >= min_duration_minutes,
      .data$duration_minutes <= 300,
      is.na(.data$effort_distance_km) | .data$effort_distance_km <= 10,
      .data$number_observers <= 10,
      !is.na(time_to_decimal(.data$time_observations_started))
    ) |>
    dplyr::mutate(
      observation_count         = as.integer(.data$observation_count),
      observation_date          = as.Date(.data$observation_date),
      day_of_year               = lubridate::yday(.data$observation_date),
      year                      = lubridate::year(.data$observation_date),
      week                      = lubridate::week(.data$observation_date),
      time_observations_started = time_to_decimal(
        .data$time_observations_started
      ),
      effort_distance_km        = dplyr::if_else(
        is.na(.data$effort_distance_km), 0, .data$effort_distance_km
      ),
      protocol_type             = factor(.data$protocol_type)
    ) |>
    # Exclude mega-flock checklists: very high counts reflect targeted non-
    # random sampling (observer sought out the flock) rather than encounter-
    # rate signal. Cap configurable via `max_count`; default 200.
    dplyr::filter(.data$observation_count <= max_count)
}

# Attach the per-observer expertise score from observer_expertise.rds onto a
# data frame containing observer_id, then drop observer_id. Observers absent
# from the calibration set are assigned the median score so they contribute
# at the population-average level. No-op if the score file isn't present.
attach_observer_expertise <- function(df, cache_dir) {
  if (!"observer_id" %in% names(df)) return(df)
  exp_f <- file.path(cache_dir, "observer_expertise.rds")
  if (!file.exists(exp_f)) {
    df$observer_id <- NULL
    return(df)
  }
  exp     <- readRDS(exp_f)
  exp_med <- stats::median(exp$expertise, na.rm = TRUE)
  df      <- merge(df, exp, by = "observer_id", all.x = TRUE)
  n_na    <- sum(is.na(df$expertise))
  df$expertise[is.na(df$expertise)] <- exp_med
  names(df)[names(df) == "expertise"] <- "observer_expertise"
  df$observer_id <- NULL
  message(sprintf(
    "  Observer expertise attached (%d / %d rows fell back to median %.3f).",
    n_na, nrow(df), exp_med
  ))
  df
}

# Build a sampling master: raw sampling rows clipped to polygon with covariates
# pre-extracted. Saved as {cache_dir}/sampling_master.rds.
#
# When the master exists, load_ebird() uses a fast path that skips
# extract_covariates(), st_filter(), and drop_na() per species — replacing
# 475 terra::extract calls with one.
#
# Call this after prepare_covariates() and the zerofill pre-cache, before
# estimate_abundance_batch().
#' @export
prepare_sampling_master <- function(sampling_txt, polygon, cov_stack, cache_dir) {
  master_f <- file.path(cache_dir, "sampling_master.rds")
  if (file.exists(master_f)) {
    message("Sampling master already exists: ", master_f)
    return(invisible(readRDS(master_f)))
  }

  polygon_wgs84 <- sf::st_transform(polygon, 4326)
  bbox          <- as.numeric(sf::st_bbox(polygon_wgs84))

  message("Building sampling master (one-time covariate extraction for all species)...")
  samp <- read_sampling(sampling_txt, bbox)

  # Apply shared-pool X-count filter: drop any checklist that had X for any
  # species in the EBD, so every species sees an identical sampling pool.
  # The IDs are written by the pre-cache step in run_batch_*.R.
  x_ids_f <- file.path(cache_dir, "x_count_ids.rds")
  if (file.exists(x_ids_f)) {
    x_ids    <- readRDS(x_ids_f)
    n_before <- nrow(samp)
    samp     <- samp[!samp$checklist_id %in% x_ids, ]
    if (nrow(samp) < n_before)
      message(sprintf("  Dropped %d X-count checklists (shared-pool filter).",
                      n_before - nrow(samp)))
  }

  # Spatial clip (keeps only checklists inside the polygon, not just the bbox)
  samp_sf <- sf::st_as_sf(samp, coords = c("longitude", "latitude"),
                           crs = 4326, remove = FALSE)
  samp    <- sf::st_drop_geometry(sf::st_filter(samp_sf, polygon_wgs84))

  # Attach observer expertise scores (and drop observer_id) if Stage 1
  # calibration has been run (observer_expertise.rds in cache_dir).
  samp <- attach_observer_expertise(samp, cache_dir)

  # Extract covariates once for all checklist locations
  samp <- extract_covariates(samp, cov_stack)

  # Drop rows where covariates are missing (ocean edges, raster gaps)
  cov_cols <- grep(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height|nightlights|palsar_hv)",
    names(samp), value = TRUE
  )
  n_before <- nrow(samp)
  samp     <- tidyr::drop_na(samp, dplyr::all_of(cov_cols))
  if (nrow(samp) < n_before)
    message(sprintf("  Dropped %d checklists with missing covariates (%d remaining).",
                    n_before - nrow(samp), nrow(samp)))

  saveRDS(samp, master_f)
  message(sprintf("Sampling master saved: %d checklists → %s", nrow(samp), master_f))
  invisible(samp)
}

# Main loader: read, zero-fill, clean, clip to polygon.
# If a sampling_master.rds exists in cache_dir, uses a fast path that skips
# extract_covariates(), st_filter(), and drop_na(cov_cols) — the expensive
# steps that are identical across all species.
# Returns a flat data.frame with one row per checklist.
load_ebird <- function(polygon, ebird_zip, sampling_txt, species, cache_dir,
                       max_count = 200L) {
  polygon_wgs84 <- sf::st_transform(polygon, 4326)
  bb <- as.numeric(sf::st_bbox(polygon_wgs84))

  spp      <- safe_name(species)
  cache_f  <- file.path(cache_dir, sprintf("zerofilled_%s.rds", spp))
  master_f <- file.path(cache_dir, "sampling_master.rds")
  x_ids_f  <- file.path(cache_dir, "x_count_ids.rds")

  # ── Fast path: sampling master available ─────────────────────────────────
  if (file.exists(master_f) && file.exists(cache_f)) {
    message("Loading sampling master + per-species obs for '", species, "'.")
    master  <- readRDS(master_f)
    obs     <- readRDS(cache_f)[, c("checklist_id", "observation_count"),
                                drop = FALSE]
    zf      <- merge(master, obs, by = "checklist_id", all.x = TRUE)
    zf$observation_count[is.na(zf$observation_count)] <- "0"
    zf$species_observed <- zf$observation_count != "0"
    zf <- clean_ebird(zf, max_count = max_count)
    message(sprintf("Loaded %d checklists (%d with detections) inside polygon.",
                    nrow(zf), sum(zf$species_observed)))
    return(zf)
  }

  # ── Slow path: build from scratch ────────────────────────────────────────
  if (file.exists(cache_f)) {
    message("Loading cached zero-filled data for '", species, "'.")
    zf <- readRDS(cache_f)
  } else {
    ebd_txt  <- vapply(ebird_zip, resolve_ebird_path, character(1))
    sampling <- read_sampling(sampling_txt, bb)
    # Shared-pool X-count filter (if pre-cache has populated it)
    if (file.exists(x_ids_f)) {
      x_ids    <- readRDS(x_ids_f)
      n_before <- nrow(sampling)
      sampling <- sampling[!sampling$checklist_id %in% x_ids, ]
      if (nrow(sampling) < n_before)
        message(sprintf("  Dropped %d X-count checklists (shared-pool filter).",
                        n_before - nrow(sampling)))
    }
    ebd_spp  <- read_ebd_species(ebd_txt, species)
    zf       <- zero_fill(sampling, ebd_spp)
    saveRDS(zf, cache_f)
  }

  zf <- clean_ebird(zf, max_count = max_count)
  zf <- attach_observer_expertise(zf, cache_dir)

  if (nrow(zf) == 0) {
    stop(
      "No usable checklists found for '", species,
      "' within the polygon after filtering."
    )
  }

  zf_sf  <- sf::st_as_sf(
    zf, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
  )
  zf_sf  <- sf::st_filter(zf_sf, polygon_wgs84)
  zf_out <- sf::st_drop_geometry(zf_sf)

  message(sprintf(
    "Loaded %d checklists (%d with detections) inside polygon.",
    nrow(zf_out), sum(zf_out$species_observed)
  ))

  zf_out
}

# Resolve the EBD path — accepts either a .txt or a .zip.
# If a .zip is given, looks for the extracted .txt alongside it or in a
# same-named subdirectory (the layout 'unzip' produces).
resolve_ebird_path <- function(ebird_path) {
  if (grepl("\\.txt$", ebird_path, ignore.case = TRUE)) {
    if (!file.exists(ebird_path)) stop("EBD file not found: ", ebird_path)
    return(ebird_path)
  }

  zip_dir <- dirname(ebird_path)
  base    <- tools::file_path_sans_ext(basename(ebird_path))

  candidates <- c(
    file.path(zip_dir, paste0(base, ".txt")),
    file.path(zip_dir, base, paste0(base, ".txt"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) return(found[1])

  stop(
    "Could not find the extracted EBD .txt alongside the zip.\n",
    "Please extract manually, e.g.:\n",
    "  unzip ", ebird_path, " -d ", zip_dir, "\n",
    "Or pass the .txt path directly to ebird_zip."
  )
}
