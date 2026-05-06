test_that("resolve_ebird_path returns an existing .txt path unchanged", {
  tmp <- tempfile(fileext = ".txt")
  writeLines("header\trow", tmp)
  expect_equal(ebirdabund:::resolve_ebird_path(tmp), tmp)
  unlink(tmp)
})

test_that("resolve_ebird_path errors on missing .txt file", {
  expect_error(
    ebirdabund:::resolve_ebird_path("/nonexistent/path/file.txt"),
    "EBD file not found"
  )
})

test_that("resolve_ebird_path errors when extracted .txt not found for .zip", {
  tmp_zip <- tempfile(fileext = ".zip")
  writeLines("", tmp_zip)
  expect_error(
    ebirdabund:::resolve_ebird_path(tmp_zip),
    "Could not find"
  )
  unlink(tmp_zip)
})

test_that("resolve_ebird_path finds .txt alongside .zip", {
  tmp_dir <- tempdir()
  base    <- "ebd_test_rel2025"
  zip_f   <- file.path(tmp_dir, paste0(base, ".zip"))
  txt_f   <- file.path(tmp_dir, paste0(base, ".txt"))
  writeLines("", zip_f)
  writeLines("header", txt_f)
  expect_equal(ebirdabund:::resolve_ebird_path(zip_f), txt_f)
  unlink(c(zip_f, txt_f))
})

test_that("zero_fill merges and assigns zeros for undetected checklists", {
  sampling_df <- data.frame(
    checklist_id = c("S1", "S2", "S3"),
    stringsAsFactors = FALSE
  )
  ebd_df <- data.frame(
    checklist_id      = "S1",
    observation_count = "3",
    stringsAsFactors  = FALSE
  )
  result <- ebirdabund:::zero_fill(sampling_df, ebd_df)

  expect_equal(nrow(result), 3)
  expect_equal(
    result$observation_count[result$checklist_id == "S1"], "3"
  )
  expect_equal(
    result$observation_count[result$checklist_id == "S2"], "0"
  )
  expect_true(result$species_observed[result$checklist_id == "S1"])
  expect_false(result$species_observed[result$checklist_id == "S2"])
})

test_that("zero_fill returns all sampling rows even with no detections", {
  sampling_df <- data.frame(checklist_id = c("S1", "S2"))
  ebd_df      <- data.frame(checklist_id = character(0),
                             observation_count = character(0))
  result <- ebirdabund:::zero_fill(sampling_df, ebd_df)
  expect_equal(nrow(result), 2)
  expect_true(all(result$observation_count == "0"))
})

test_that("clean_ebird removes X counts, short/long durations, far-ranging efforts", {
  df <- data.frame(
    observation_count         = c("2", "0", "X", "1",  "1",   "1"),
    duration_minutes          = c(60,  60,  60,  3,    400,   60),
    effort_distance_km        = c(1,   1,   1,   1,    1,     15),
    number_observers          = rep(1, 6),
    protocol_type             = "Traveling Count",
    observation_date          = as.Date("2020-06-15"),
    time_observations_started = "07:00",
    stringsAsFactors          = FALSE
  )
  # rows 3 (X), 4 (< 5 min), 5 (> 300 min), 6 (> 10 km) should be dropped
  result <- ebirdabund:::clean_ebird(df)
  expect_equal(nrow(result), 2)
})

test_that("clean_ebird converts column types correctly", {
  df <- data.frame(
    observation_count         = "5",
    duration_minutes          = 60,
    effort_distance_km        = 1,
    number_observers          = 1,
    protocol_type             = "Traveling Count",
    observation_date          = as.Date("2020-06-15"),
    time_observations_started = "07:30",
    stringsAsFactors          = FALSE
  )
  result <- ebirdabund:::clean_ebird(df)
  expect_true(is.integer(result$observation_count))
  expect_true(inherits(result$observation_date, "Date"))
  expect_true(is.factor(result$protocol_type))
  expect_true(is.numeric(result$time_observations_started))
  expect_true("day_of_year" %in% names(result))
})

test_that("clean_ebird allows NA effort_distance_km (stationary counts)", {
  df <- data.frame(
    observation_count         = "2",
    duration_minutes          = 60,
    effort_distance_km        = NA_real_,
    number_observers          = 1,
    protocol_type             = "Stationary Count",
    observation_date          = as.Date("2020-06-15"),
    time_observations_started = "06:00",
    stringsAsFactors          = FALSE
  )
  expect_equal(nrow(ebirdabund:::clean_ebird(df)), 1)
})

# ── validate_ebd_header ─────────────────────────────────────────────────────

# Build a header line that matches the eBird EBD layout up to col 35 so the
# validator's index check finds the expected names. The exact contents of
# unused columns don't matter — only the positions of cols 6, 11, and 35.
make_valid_ebd_header <- function() {
  cols <- character(35L)
  cols[1:5]   <- c("GLOBAL UNIQUE IDENTIFIER", "LAST EDITED DATE",
                   "TAXONOMIC ORDER", "CATEGORY", "TAXON CONCEPT ID")
  cols[6]     <- "COMMON NAME"
  cols[7:10]  <- c("SCIENTIFIC NAME", "SUBSPECIES COMMON NAME",
                   "SUBSPECIES SCIENTIFIC NAME", "EXOTIC CODE")
  cols[11]    <- "OBSERVATION COUNT"
  cols[12:34] <- paste0("FILLER_", 12:34)
  cols[35]    <- "SAMPLING EVENT IDENTIFIER"
  paste(cols, collapse = "\t")
}

write_tmp_ebd <- function(header_line, body_rows = character(0)) {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(header_line, body_rows), tmp)
  tmp
}

test_that("validate_ebd_header passes on a correct header", {
  tmp <- write_tmp_ebd(make_valid_ebd_header())
  on.exit(unlink(tmp), add = TRUE)
  expect_true(ebirdabund:::validate_ebd_header(tmp))
})

test_that("validate_ebd_header errors when col 6 is wrong", {
  cols <- strsplit(make_valid_ebd_header(), "\t", fixed = TRUE)[[1L]]
  cols[6] <- "WRONG NAME"
  tmp <- write_tmp_ebd(paste(cols, collapse = "\t"))
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    ebirdabund:::validate_ebd_header(tmp),
    "Expected col 6 = 'COMMON NAME'; got 'WRONG NAME'"
  )
})

test_that("validate_ebd_header errors when col 35 is missing entirely", {
  cols <- strsplit(make_valid_ebd_header(), "\t", fixed = TRUE)[[1L]]
  short <- cols[1:30]   # truncate before col 35
  tmp <- write_tmp_ebd(paste(short, collapse = "\t"))
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    ebirdabund:::validate_ebd_header(tmp),
    "Expected col 35.*<missing>"
  )
})

test_that("validate_ebd_header errors on empty file", {
  tmp <- tempfile(fileext = ".txt")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    ebirdabund:::validate_ebd_header(tmp),
    "empty or unreadable"
  )
})

# ── read_sampling: multi-state + date cutoff ────────────────────────────────

# Build a tab-separated sampling file with the columns read_sampling expects.
make_tmp_sampling <- function(rows) {
  cols <- c("SAMPLING EVENT IDENTIFIER", "ALL SPECIES REPORTED",
            "LATITUDE", "LONGITUDE", "OBSERVATION DATE",
            "TIME OBSERVATIONS STARTED", "DURATION MINUTES",
            "EFFORT DISTANCE KM", "NUMBER OBSERVERS", "PROTOCOL NAME")
  body <- vapply(rows, function(r) paste(r[cols], collapse = "\t"),
                 character(1L))
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(paste(cols, collapse = "\t"), body), tmp)
  tmp
}

# Helper: one sampling row with sensible defaults — caller can override fields.
samp_row <- function(...) {
  defaults <- list(
    `SAMPLING EVENT IDENTIFIER` = "S0",
    `ALL SPECIES REPORTED`      = "1",
    `LATITUDE`                  = "-35.0",
    `LONGITUDE`                 = "149.5",
    `OBSERVATION DATE`          = "2024-06-01",
    `TIME OBSERVATIONS STARTED` = "07:00:00",
    `DURATION MINUTES`          = "60",
    `EFFORT DISTANCE KM`        = "1.0",
    `NUMBER OBSERVERS`          = "1",
    `PROTOCOL NAME`             = "Traveling"
  )
  modifyList(defaults, list(...))
}

test_that("read_sampling errors on an empty path vector", {
  expect_error(
    ebirdabund:::read_sampling(character(0), c(140, -40, 155, -30)),
    "non-empty character vector"
  )
})

test_that("read_sampling concatenates multiple sampling files", {
  tmp1 <- make_tmp_sampling(list(
    samp_row(`SAMPLING EVENT IDENTIFIER` = "S1"),
    samp_row(`SAMPLING EVENT IDENTIFIER` = "S2")
  ))
  tmp2 <- make_tmp_sampling(list(
    samp_row(`SAMPLING EVENT IDENTIFIER` = "S3", `LATITUDE` = "-34.0",
             `LONGITUDE` = "148.0")
  ))
  on.exit(unlink(c(tmp1, tmp2)), add = TRUE)

  out <- ebirdabund:::read_sampling(c(tmp1, tmp2),
                                    bbox = c(140, -40, 155, -30))

  expect_equal(nrow(out), 3L)
  expect_setequal(out$checklist_id, c("S1", "S2", "S3"))
  expect_equal(
    sort(names(out)),
    sort(c("checklist_id", "latitude", "longitude", "observation_date",
           "time_observations_started", "duration_minutes",
           "effort_distance_km", "number_observers", "protocol_type"))
  )
})

test_that("read_sampling drops checklists outside bbox", {
  tmp <- make_tmp_sampling(list(
    samp_row(`SAMPLING EVENT IDENTIFIER` = "INSIDE",
             `LATITUDE` = "-35.0", `LONGITUDE` = "149.5"),
    samp_row(`SAMPLING EVENT IDENTIFIER` = "OUTSIDE",
             `LATITUDE` = "-10.0", `LONGITUDE` = "100.0")
  ))
  on.exit(unlink(tmp), add = TRUE)

  out <- ebirdabund:::read_sampling(tmp, bbox = c(140, -40, 155, -30))
  expect_equal(out$checklist_id, "INSIDE")
})

test_that("read_sampling drops rows with ALL SPECIES REPORTED == 0", {
  tmp <- make_tmp_sampling(list(
    samp_row(`SAMPLING EVENT IDENTIFIER` = "OK"),
    samp_row(`SAMPLING EVENT IDENTIFIER` = "INCOMPLETE",
             `ALL SPECIES REPORTED` = "0")
  ))
  on.exit(unlink(tmp), add = TRUE)

  out <- ebirdabund:::read_sampling(tmp, bbox = c(140, -40, 155, -30))
  expect_equal(out$checklist_id, "OK")
})

test_that("read_sampling applies date_cutoff inclusively", {
  tmp <- make_tmp_sampling(list(
    samp_row(`SAMPLING EVENT IDENTIFIER` = "BEFORE",
             `OBSERVATION DATE` = "2024-01-15"),
    samp_row(`SAMPLING EVENT IDENTIFIER` = "ON_CUTOFF",
             `OBSERVATION DATE` = "2024-06-30"),
    samp_row(`SAMPLING EVENT IDENTIFIER` = "AFTER",
             `OBSERVATION DATE` = "2024-07-01")
  ))
  on.exit(unlink(tmp), add = TRUE)

  out <- ebirdabund:::read_sampling(
    tmp, bbox = c(140, -40, 155, -30),
    date_cutoff = as.Date("2024-06-30")
  )

  expect_setequal(out$checklist_id, c("BEFORE", "ON_CUTOFF"))
})

test_that("read_sampling accepts character date_cutoff", {
  tmp <- make_tmp_sampling(list(
    samp_row(`SAMPLING EVENT IDENTIFIER` = "OK",
             `OBSERVATION DATE` = "2024-01-15"),
    samp_row(`SAMPLING EVENT IDENTIFIER` = "FUTURE",
             `OBSERVATION DATE` = "2030-01-01")
  ))
  on.exit(unlink(tmp), add = TRUE)

  out <- ebirdabund:::read_sampling(
    tmp, bbox = c(140, -40, 155, -30),
    date_cutoff = "2024-12-31"
  )
  expect_equal(out$checklist_id, "OK")
})

test_that("read_sampling errors on unparseable date_cutoff", {
  tmp <- make_tmp_sampling(list(samp_row()))
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    ebirdabund:::read_sampling(tmp, bbox = c(140, -40, 155, -30),
                               date_cutoff = "not-a-date"),
    "Date"
  )
})

# ── read_ebd_observations ───────────────────────────────────────────────────

# Build a tab-separated EBD file with the 35-column header layout. Each row
# is fully synthetic (filler values for unused columns) so cols 6/11/35 can
# be tested without depending on the real EBD.
make_tmp_ebd <- function(rows) {
  header <- make_valid_ebd_header()
  body <- vapply(rows, function(r) {
    cols <- character(35L)
    cols[]   <- ""              # filler defaults
    cols[6]  <- r[["common_name"]]
    cols[11] <- r[["observation_count"]]
    cols[35] <- r[["checklist_id"]]
    paste(cols, collapse = "\t")
  }, character(1L))
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(header, body), tmp)
  tmp
}

test_that("read_ebd_observations validates header before reading", {
  cols <- strsplit(make_valid_ebd_header(), "\t", fixed = TRUE)[[1L]]
  cols[6] <- "WRONG NAME"
  tmp <- write_tmp_ebd(paste(cols, collapse = "\t"))
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    ebirdabund:::read_ebd_observations(tmp),
    "EBD header check failed"
  )
})

test_that("read_ebd_observations concatenates multiple EBD files", {
  tmp1 <- make_tmp_ebd(list(
    list(common_name = "Sp A", observation_count = "2", checklist_id = "S1"),
    list(common_name = "Sp B", observation_count = "1", checklist_id = "S1")
  ))
  tmp2 <- make_tmp_ebd(list(
    list(common_name = "Sp A", observation_count = "5", checklist_id = "S2")
  ))
  on.exit(unlink(c(tmp1, tmp2)), add = TRUE)

  out <- ebirdabund:::read_ebd_observations(c(tmp1, tmp2))

  expect_equal(nrow(out), 3L)
  expect_setequal(out$common_name, c("Sp A", "Sp B"))
  expect_setequal(out$checklist_id, c("S1", "S2"))
})

test_that("read_ebd_observations filters by species_set", {
  tmp <- make_tmp_ebd(list(
    list(common_name = "Sp A", observation_count = "2", checklist_id = "S1"),
    list(common_name = "Sp B", observation_count = "1", checklist_id = "S2"),
    list(common_name = "Sp C", observation_count = "3", checklist_id = "S3")
  ))
  on.exit(unlink(tmp), add = TRUE)

  out <- ebirdabund:::read_ebd_observations(tmp,
                                            species_set = c("Sp A", "Sp C"))
  expect_equal(nrow(out), 2L)
  expect_setequal(out$common_name, c("Sp A", "Sp C"))
})

test_that("read_ebd_observations filters by valid_checklist_ids", {
  tmp <- make_tmp_ebd(list(
    list(common_name = "Sp A", observation_count = "2", checklist_id = "S1"),
    list(common_name = "Sp A", observation_count = "5", checklist_id = "S2")
  ))
  on.exit(unlink(tmp), add = TRUE)

  out <- ebirdabund:::read_ebd_observations(tmp,
                                            valid_checklist_ids = "S2")
  expect_equal(out$checklist_id, "S2")
  expect_equal(out$observation_count, "5")
})

# ── Real-data smoke test: ACT (smallest state) ──────────────────────────────
# Skipped in CI / portable runs where the raw EBD isn't present.

test_that("read_sampling and validate_ebd_header work on the real ACT dump", {
  act_dir   <- "ebirdabund/raw_data/ebd_AU-ACT_unv_smp_relMar-2026"
  act_dir   <- file.path(rprojroot::find_package_root_file(), "..", act_dir)
  act_samp  <- file.path(act_dir, "ebd_AU-ACT_unv_smp_relMar-2026_sampling.txt")
  act_ebd   <- file.path(act_dir, "ebd_AU-ACT_unv_smp_relMar-2026.txt")

  skip_if_not(file.exists(act_samp), "ACT sampling file not present")
  skip_if_not(file.exists(act_ebd),  "ACT EBD file not present")

  expect_true(ebirdabund:::validate_ebd_header(act_ebd))

  # ACT bbox roughly: lon 148.7-149.4, lat -35.93 to -35.12
  out <- ebirdabund:::read_sampling(act_samp,
                                    bbox = c(148.5, -36.5, 149.5, -35.0))
  expect_gt(nrow(out), 0L)
  expect_true(all(c("checklist_id", "latitude", "longitude",
                    "observation_date") %in% names(out)))
  expect_true(all(out$latitude  >= -36.5 & out$latitude  <= -35.0))
  expect_true(all(out$longitude >= 148.5 & out$longitude <= 149.5))
})
