test_that("time_to_decimal converts HH:MM correctly", {
  expect_equal(ebirdabund:::time_to_decimal("06:30"), 6.5)
  expect_equal(ebirdabund:::time_to_decimal("12:00"), 12.0)
  expect_equal(ebirdabund:::time_to_decimal("00:00"), 0.0)
})

test_that("time_to_decimal handles HH:MM:SS", {
  expect_equal(
    ebirdabund:::time_to_decimal("06:30:30"),
    6 + 30 / 60 + 30 / 3600
  )
})

test_that("time_to_decimal handles NA and empty string", {
  expect_true(is.na(ebirdabund:::time_to_decimal(NA_character_)))
  expect_true(is.na(ebirdabund:::time_to_decimal("")))
})

test_that("time_to_decimal is vectorised and preserves length", {
  result <- ebirdabund:::time_to_decimal(c("06:00", "12:30", NA))
  expect_length(result, 3)
  expect_equal(result[1], 6.0)
  expect_equal(result[2], 12.5)
  expect_true(is.na(result[3]))
})

test_that("safe_name produces lowercase alphanumeric-and-underscore strings", {
  expect_equal(ebirdabund:::safe_name("Superb Fairywren"), "superb_fairywren")
  expect_equal(ebirdabund:::safe_name("Red-tailed Black-Cockatoo"),
               "red_tailed_black_cockatoo")
  expect_match(ebirdabund:::safe_name("Any Species Name!"), "^[a-z0-9_]+$")
})

test_that("safe_name trims leading/trailing whitespace", {
  expect_equal(ebirdabund:::safe_name("  Superb Fairywren  "),
               "superb_fairywren")
})
