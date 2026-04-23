# Input validation — these tests do not need real eBird files.

test_that("fit_species_model rejects non-sf polygon", {
  expect_error(
    fit_species_model("not_a_polygon", "x.txt", "y.txt", "Species"),
    "must be an sf or sfc object"
  )
})

test_that("fit_species_model errors on missing eBird file", {
  poly <- make_test_polygon()
  expect_error(
    fit_species_model(poly, "/nonexistent/ebd.txt", "y.txt", "Species"),
    "eBird data file not found"
  )
})

test_that("fit_species_model errors on missing sampling file", {
  poly    <- make_test_polygon()
  tmp_ebd <- tempfile(fileext = ".txt")
  writeLines("header", tmp_ebd)
  on.exit(unlink(tmp_ebd))
  expect_error(
    fit_species_model(poly, tmp_ebd, "/nonexistent/sampling.txt", "Species"),
    "Sampling events file not found"
  )
})

test_that("fit_species_model errors on empty species name", {
  poly    <- make_test_polygon()
  tmp_ebd <- tempfile(fileext = ".txt")
  tmp_smp <- tempfile(fileext = ".txt")
  writeLines("header", tmp_ebd)
  writeLines("header", tmp_smp)
  on.exit(unlink(c(tmp_ebd, tmp_smp)))
  expect_error(
    fit_species_model(poly, tmp_ebd, tmp_smp, ""),
    "non-empty character string"
  )
})

# predict_species_map — use synthetic model_fit to avoid file I/O.

test_that("predict_species_map returns predictions and plot", {
  polygon   <- make_test_polygon()
  cov_stack <- make_test_cov_stack()
  df        <- make_test_checklists()
  model     <- make_test_gam(df)
  cov_cols  <- grep("^(lc_|elevation)", names(df), value = TRUE)

  model_fit <- list(model     = model,
                    data      = df,
                    cov_stack = cov_stack,
                    cov_cols  = cov_cols)

  result <- suppressMessages(suppressWarnings(
    predict_species_map(model_fit, polygon,
                        species     = "Test Species",
                        grid_res_km = 10,
                        use_range   = FALSE)
  ))

  expect_named(result, c("predictions", "plot"))
  expect_true(inherits(result$predictions, "SpatRaster"))
  expect_true(inherits(result$plot, "ggplot"))
})

test_that("predict_species_map with use_range=FALSE skips ebirdst lookup", {
  polygon   <- make_test_polygon()
  cov_stack <- make_test_cov_stack()
  df        <- make_test_checklists()
  model     <- make_test_gam(df)
  cov_cols  <- grep("^(lc_|elevation)", names(df), value = TRUE)

  model_fit <- list(model     = model,
                    data      = df,
                    cov_stack = cov_stack,
                    cov_cols  = cov_cols)

  # Should not error even with a bogus species name because use_range = FALSE
  expect_no_error(suppressMessages(
    predict_species_map(model_fit, polygon,
                        species     = "Definitely Not A Real Species 99",
                        grid_res_km = 10,
                        use_range   = FALSE)
  ))
})

# estimate_abundance — same validation surface as fit_species_model.

test_that("estimate_abundance rejects non-sf polygon", {
  expect_error(
    estimate_abundance("not_sf", "x.txt", "y.txt", "Species"),
    "must be an sf or sfc object"
  )
})
