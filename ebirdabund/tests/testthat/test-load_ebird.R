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
