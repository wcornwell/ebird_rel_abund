test_that("make_prediction_surface returns a data.frame with longitude/latitude", {
  polygon <- make_test_polygon()
  surface <- suppressMessages(
    ebirdabund:::make_prediction_surface(polygon, grid_res_km = 10)
  )
  expect_true(is.data.frame(surface))
  expect_true(all(c("longitude", "latitude") %in% names(surface)))
})

test_that("make_prediction_surface points fall inside the polygon", {
  polygon <- make_test_polygon()
  surface <- suppressMessages(
    ebirdabund:::make_prediction_surface(polygon, grid_res_km = 10)
  )
  expect_true(all(surface$longitude >= 150 & surface$longitude <= 151))
  expect_true(all(surface$latitude  >= -34 & surface$latitude  <= -33))
  expect_gt(nrow(surface), 0)
})

test_that("make_prediction_surface finer grid produces more points", {
  polygon  <- make_test_polygon()
  coarse   <- suppressMessages(ebirdabund:::make_prediction_surface(polygon, 20))
  fine     <- suppressMessages(ebirdabund:::make_prediction_surface(polygon, 5))
  expect_gt(nrow(fine), nrow(coarse))
})

test_that("extract_covariates appends covariate columns without altering row count", {
  cov_stack <- make_test_cov_stack()
  df <- data.frame(longitude = c(150.1, 150.5, 150.9),
                   latitude  = c(-33.2, -33.5, -33.8))
  result <- ebirdabund:::extract_covariates(df, cov_stack)
  expect_equal(nrow(result), nrow(df))
  expect_true(all(c("lc_trees", "lc_grassland", "elevation") %in% names(result)))
  expect_true(all(c("longitude", "latitude") %in% names(result)))
})
