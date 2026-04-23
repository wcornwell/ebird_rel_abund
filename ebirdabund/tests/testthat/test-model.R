test_that("safe_k returns an integer", {
  expect_type(ebirdabund:::safe_k(1:10, 5), "integer")
})

test_that("safe_k respects default_k upper bound", {
  expect_equal(ebirdabund:::safe_k(1:20, 5), 5L)
})

test_that("safe_k enforces minimum of 3", {
  expect_equal(ebirdabund:::safe_k(1:2, 5), 3L)   # n_uniq=2, n_uniq-1=1 → clamped to 3
  expect_equal(ebirdabund:::safe_k(1:3, 5), 3L)   # n_uniq=3, n_uniq-1=2 → clamped to 3
})

test_that("safe_k uses n_unique - 1 when between 3 and default_k", {
  expect_equal(ebirdabund:::safe_k(1:5, 10), 4L)  # n_uniq=5, n_uniq-1=4
})

test_that("fit_gam returns an mgcv gam object", {
  df <- make_test_checklists()
  model <- make_test_gam(df)
  expect_true(inherits(model, "gam"))
})

test_that("fit_gam errors when no habitat covariate columns present", {
  df <- make_test_checklists()
  df <- df[, !grepl("^(lc_|elevation)", names(df))]
  expect_error(
    suppressMessages(ebirdabund:::fit_gam(df)),
    "No habitat covariate columns"
  )
})

test_that("fit_gam deviance explained is between 0 and 1", {
  model <- make_test_gam()
  dev_expl <- summary(model)$dev.expl
  expect_gte(dev_expl, 0)
  expect_lte(dev_expl, 1)
})

test_that("fit_gam drops hab columns with fewer than 4 unique values", {
  df <- make_test_checklists()
  df$lc_trees <- rep(c(0.1, 0.2, 0.3), length.out = nrow(df))  # only 3 unique
  model <- make_test_gam(df)
  expect_false(any(grepl("lc_trees", as.character(formula(model)))))
})
