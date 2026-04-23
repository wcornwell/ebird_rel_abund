test_that("find_peak_time returns a value in [0, 24]", {
  df    <- make_test_checklists()
  model <- make_test_gam(df)

  template <- df[1, ]
  template$day_of_year        <- 180
  template$duration_minutes   <- 60
  template$effort_distance_km <- 1
  template$number_observers   <- 1

  peak <- suppressMessages(ebirdabund:::find_peak_time(model, template))
  expect_true(is.numeric(peak))
  expect_length(peak, 1)
  expect_gte(peak, 0)
  expect_lte(peak, 24)
})

test_that("ref_protocol returns 'Traveling Count' when present", {
  model <- make_test_gam()
  expect_equal(ebirdabund:::ref_protocol(model), "Traveling Count")
})

test_that("predict_abundance returns a SpatRaster with abd and abd_se layers", {
  polygon   <- make_test_polygon()
  cov_stack <- make_test_cov_stack()
  df        <- make_test_checklists()
  model     <- make_test_gam(df)

  surface  <- suppressMessages(
    ebirdabund:::make_prediction_surface(polygon, 10)
  )
  surface  <- ebirdabund:::extract_covariates(surface, cov_stack)
  cov_cols <- grep("^(lc_|elevation)", names(surface), value = TRUE)
  surface  <- tidyr::drop_na(surface, dplyr::all_of(cov_cols))

  r <- suppressMessages(
    ebirdabund:::predict_abundance(model, surface, polygon,
                                   grid_res_km = 10, peak_time = 7)
  )
  expect_true(inherits(r, "SpatRaster"))
  expect_equal(terra::nlyr(r), 2)
  expect_true(all(c("abd", "abd_se") %in% names(r)))
})

test_that("predict_abundance abd values are non-negative", {
  polygon   <- make_test_polygon()
  cov_stack <- make_test_cov_stack()
  model     <- make_test_gam()

  surface  <- suppressMessages(
    ebirdabund:::make_prediction_surface(polygon, 10)
  )
  surface  <- ebirdabund:::extract_covariates(surface, cov_stack)
  cov_cols <- grep("^(lc_|elevation)", names(surface), value = TRUE)
  surface  <- tidyr::drop_na(surface, dplyr::all_of(cov_cols))

  r    <- suppressMessages(
    ebirdabund:::predict_abundance(model, surface, polygon,
                                   grid_res_km = 10, peak_time = 7)
  )
  vals <- terra::values(r[["abd"]], na.rm = TRUE)
  expect_true(all(vals >= 0))
})

test_that("predict_abundance output has correct CRS and extent", {
  polygon   <- make_test_polygon()
  cov_stack <- make_test_cov_stack()
  model     <- make_test_gam()

  surface  <- suppressMessages(
    ebirdabund:::make_prediction_surface(polygon, 10)
  )
  surface  <- ebirdabund:::extract_covariates(surface, cov_stack)
  cov_cols <- grep("^(lc_|elevation)", names(surface), value = TRUE)
  surface  <- tidyr::drop_na(surface, dplyr::all_of(cov_cols))

  r <- suppressMessages(
    ebirdabund:::predict_abundance(model, surface, polygon,
                                   grid_res_km = 10, peak_time = 7)
  )
  expect_equal(terra::crs(r, describe = TRUE)$code, "4326")

  # Raster extent should not exceed the (buffered) polygon bbox
  poly_bbox <- as.numeric(sf::st_bbox(sf::st_transform(polygon, 4326)))
  r_ext     <- as.vector(terra::ext(r))
  expect_gte(r_ext[1], poly_bbox[1] - 0.2)  # xmin
  expect_lte(r_ext[2], poly_bbox[3] + 0.2)  # xmax
})
