test_that("subsample_hex reduces checklist count", {
  df <- make_test_checklists(n = 200)
  result <- suppressMessages(ebirdabund:::subsample_hex(df, spacing_km = 50))
  expect_lt(nrow(result), nrow(df))
  expect_gt(nrow(result), 0)
})

test_that("subsample_hex does not leave a hex_cell column", {
  df <- make_test_checklists(n = 100)
  result <- suppressMessages(ebirdabund:::subsample_hex(df, spacing_km = 50))
  expect_false("hex_cell" %in% names(result))
})

test_that("subsample_hex preserves all original non-hex columns", {
  df     <- make_test_checklists(n = 100)
  result <- suppressMessages(ebirdabund:::subsample_hex(df, spacing_km = 50))
  expect_true(all(names(df) %in% names(result)))
})

test_that("subsample_hex with very large spacing retains at most one row per year/week/detection combo", {
  df <- make_test_checklists(n = 200)
  result <- suppressMessages(ebirdabund:::subsample_hex(df, spacing_km = 500))
  # With a single global hex cell, each (year, week, species_observed) triplet
  # should have at most one row
  combos <- paste(result$year, result$week, result$species_observed)
  expect_true(all(!duplicated(combos)))
})
