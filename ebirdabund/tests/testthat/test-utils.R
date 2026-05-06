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

test_that("make_fit_polygon returns sf in WGS84 and grows the bbox", {
  poly <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = 149, ymin = -36, xmax = 150, ymax = -35), crs = 4326)
  ) |> sf::st_as_sf()

  buf <- ebirdabund:::make_fit_polygon(poly, buffer_km = 50)

  expect_s3_class(buf, c("sf", "sfc"))
  expect_equal(sf::st_crs(buf)$epsg, 4326L)

  # 50 km ≈ 0.45° latitude — buffered bbox should extend the original by ~0.4°.
  bb_in  <- sf::st_bbox(poly)
  bb_out <- sf::st_bbox(buf)
  expect_lt(bb_out[["xmin"]], bb_in[["xmin"]] - 0.3)
  expect_gt(bb_out[["xmax"]], bb_in[["xmax"]] + 0.3)
  expect_lt(bb_out[["ymin"]], bb_in[["ymin"]] - 0.3)
  expect_gt(bb_out[["ymax"]], bb_in[["ymax"]] + 0.3)
})

test_that("make_fit_polygon with buffer_km = 0 leaves the area effectively unchanged", {
  poly <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = 149, ymin = -36, xmax = 150, ymax = -35), crs = 4326)
  ) |> sf::st_as_sf()

  buf <- ebirdabund:::make_fit_polygon(poly, buffer_km = 0)

  a_in  <- as.numeric(sf::st_area(sf::st_transform(poly, "EPSG:3577")))
  a_out <- as.numeric(sf::st_area(sf::st_transform(buf,  "EPSG:3577")))
  # Round-trip through 3577 + small reprojection noise; tolerate < 0.5 %.
  expect_lt(abs(a_out - a_in) / a_in, 0.005)
})

test_that("make_fit_polygon clips to land_polygon when supplied", {
  # Polygon partially over open ocean (east of NSW coast):
  poly <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = 152, ymin = -34, xmax = 154, ymax = -33), crs = 4326)
  ) |> sf::st_as_sf()
  # Synthetic "land" — only the western strip is land:
  land <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = 150, ymin = -35, xmax = 153, ymax = -32), crs = 4326)
  ) |> sf::st_as_sf()

  buf <- ebirdabund:::make_fit_polygon(poly, buffer_km = 10, land_polygon = land)

  expect_s3_class(buf, c("sf", "sfc"))
  bb <- sf::st_bbox(buf)
  # Eastern edge should be clipped to the land polygon's eastern bound (153).
  expect_lt(bb[["xmax"]], 153.05)
})

test_that("make_fit_polygon validates inputs", {
  expect_error(ebirdabund:::make_fit_polygon("not_sf", 10), "must be an sf")
  poly <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = 149, ymin = -36, xmax = 150, ymax = -35), crs = 4326)
  ) |> sf::st_as_sf()
  expect_error(ebirdabund:::make_fit_polygon(poly, -5), "non-negative")
  expect_error(ebirdabund:::make_fit_polygon(poly, c(1, 2)), "non-negative")
  expect_error(ebirdabund:::make_fit_polygon(poly, NA_real_), "non-negative")
})
