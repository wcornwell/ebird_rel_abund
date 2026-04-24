# Read the sampling events file (316 MB), filter to complete checklists within
# the bounding box, and return a plain data.frame with renamed columns.
#
# Strategy: read all columns (avoids fread column-order ambiguity), then
# filter and select with standard R [[ ]] — no data.table NSE needed.
read_sampling <- function(sampling_txt, bbox) {
  message("Reading sampling events file...")
  df <- as.data.frame(data.table::fread(
    sampling_txt, sep = "\t", quote = "", showProgress = FALSE, na.strings = ""
  ))

  # Standard R subsetting — no NSE, no backtick issues
  df <- df[df[["ALL SPECIES REPORTED"]] == 1 &
           !is.na(df[["LATITUDE"]]) &
           df[["LATITUDE"]]  >= bbox[2] & df[["LATITUDE"]]  <= bbox[4] &
           df[["LONGITUDE"]] >= bbox[1] & df[["LONGITUDE"]] <= bbox[3], ]

  raw <- c("SAMPLING EVENT IDENTIFIER", "LATITUDE", "LONGITUDE",
           "OBSERVATION DATE", "TIME OBSERVATIONS STARTED",
           "DURATION MINUTES", "EFFORT DISTANCE KM",
           "NUMBER OBSERVERS", "PROTOCOL NAME")
  df  <- df[, raw, drop = FALSE]
  names(df) <- c("checklist_id", "latitude", "longitude",
                 "observation_date", "time_observations_started",
                 "duration_minutes", "effort_distance_km",
                 "number_observers", "protocol_type")
  df
}

# Scan the EBD observations file for a single species.
# Uses integer column indices (not names) to avoid fread's select-vs-file-order
# ambiguity. Positions verified against the Feb-2026 EBD header:
#   col 6  = COMMON NAME
#   col 11 = OBSERVATION COUNT
#   col 35 = SAMPLING EVENT IDENTIFIER
read_ebd_species <- function(ebd_txt, species) {
  message("Scanning EBD file for '", species, "' (large file, please wait)...")
  df <- as.data.frame(data.table::fread(
    ebd_txt,
    select       = c(6L, 11L, 35L),
    sep          = "\t",
    quote        = "",
    showProgress = TRUE,
    na.strings   = ""
  ))
  # fread returns columns in file order: COMMON NAME, OBSERVATION COUNT,
  # SAMPLING EVENT IDENTIFIER
  names(df) <- c("common_name", "observation_count", "checklist_id")
  df <- df[df[["common_name"]] == species, c("checklist_id", "observation_count"),
           drop = FALSE]
  df
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
clean_ebird <- function(zf) {
  zf |>
    dplyr::filter(
      .data$observation_count != "X",
      .data$duration_minutes >= 5,
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
    )
}

# Main loader: read, zero-fill, clean, clip to polygon.
# Returns a flat data.frame with one row per checklist.
load_ebird <- function(polygon, ebird_zip, sampling_txt, species, cache_dir) {
  polygon_wgs84 <- sf::st_transform(polygon, 4326)
  bb <- as.numeric(sf::st_bbox(polygon_wgs84))  # xmin ymin xmax ymax

  spp     <- safe_name(species)
  cache_f <- file.path(cache_dir, sprintf("zerofilled_%s.rds", spp))

  if (file.exists(cache_f)) {
    message("Loading cached zero-filled data for '", species, "'.")
    zf <- readRDS(cache_f)
  } else {
    ebd_txt  <- resolve_ebird_path(ebird_zip)
    sampling <- read_sampling(sampling_txt, bb)
    ebd_spp  <- read_ebd_species(ebd_txt, species)
    zf       <- zero_fill(sampling, ebd_spp)
    saveRDS(zf, cache_f)
  }

  zf <- clean_ebird(zf)

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
