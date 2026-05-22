make_test_botw <- function() {
  # Build a minimal GeoPackage with the same schema columns load_range_botw
  # queries against: sci_name, presence, origin, seasonal, geom.
  geom <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
    sf::st_polygon(list(rbind(c(2, 0), c(3, 0), c(3, 1), c(2, 1), c(2, 0)))),
    crs = 4326
  )
  rows <- sf::st_sf(
    sci_name = c("Neochmia modesta", "Charadrius mongolus"),
    presence = c(1L, 1L),
    origin   = c(1L, 1L),
    seasonal = c(1L, 1L),
    geometry = geom
  )
  path <- tempfile(fileext = ".gpkg")
  sf::st_write(rows, path, layer = "all_species", quiet = TRUE)
  path
}

test_that("load_range_botw with NULL aliases preserves prior behaviour", {
  botw <- make_test_botw()
  on.exit(unlink(botw), add = TRUE)

  # Exact-name lookup still works.
  r1 <- ebirdabund:::load_range_botw("Neochmia modesta", botw, aliases = NULL)
  expect_s3_class(r1, "sf")
  expect_equal(nrow(r1), 1L)

  # Unknown name returns NULL (no aliases consulted).
  r2 <- ebirdabund:::load_range_botw("Emblema modestum", botw, aliases = NULL)
  expect_null(r2)
})

test_that("load_range_botw rewrites the query on alias hit", {
  botw <- make_test_botw()
  on.exit(unlink(botw), add = TRUE)

  aliases <- data.frame(
    ebird_sci_name = c("Emblema modestum", "Anarhynchus mongolus"),
    botw_sci_name  = c("Neochmia modesta",  "Charadrius mongolus"),
    note           = c("genus rename",      "genus rename"),
    stringsAsFactors = FALSE
  )

  # Aliased lookup resolves to the BOTW name and returns its polygon.
  # load_range_botw selects geom only, so we verify the correct row by its
  # bbox — the two fixture polygons sit at x in [0,1] vs [2,3].
  expect_message(
    r1 <- ebirdabund:::load_range_botw("Emblema modestum", botw, aliases),
    "aliasing 'Emblema modestum' -> 'Neochmia modesta'",
    fixed = TRUE
  )
  expect_s3_class(r1, "sf")
  expect_equal(nrow(r1), 1L)
  expect_equal(unname(sf::st_bbox(r1)[c("xmin", "xmax")]), c(0, 1))

  # Second alias should resolve to the polygon at x in [2,3].
  r2 <- suppressMessages(
    ebirdabund:::load_range_botw("Anarhynchus mongolus", botw, aliases)
  )
  expect_equal(nrow(r2), 1L)
  expect_equal(unname(sf::st_bbox(r2)[c("xmin", "xmax")]), c(2, 3))
})

test_that("load_range_botw passes through names not in the alias map", {
  botw <- make_test_botw()
  on.exit(unlink(botw), add = TRUE)

  aliases <- data.frame(
    ebird_sci_name = "Emblema modestum",
    botw_sci_name  = "Neochmia modesta",
    note           = "genus rename",
    stringsAsFactors = FALSE
  )

  # Name not in alias map → exact-name BOTW lookup.
  r <- ebirdabund:::load_range_botw("Neochmia modesta", botw, aliases)
  expect_s3_class(r, "sf")
  expect_equal(nrow(r), 1L)
})

test_that("load_range_botw short-circuits when botw_sci_name is empty/NA", {
  botw <- make_test_botw()
  on.exit(unlink(botw), add = TRUE)

  aliases <- data.frame(
    ebird_sci_name = c("Ptiloris paradiseus",   "Some species"),
    botw_sci_name  = c("",                       NA_character_),
    note           = c("no BOTW equivalent",     "no BOTW equivalent"),
    stringsAsFactors = FALSE
  )

  expect_message(
    r1 <- ebirdabund:::load_range_botw("Ptiloris paradiseus", botw, aliases),
    "no BOTW equivalent",
    fixed = TRUE
  )
  expect_null(r1)

  expect_message(
    r2 <- ebirdabund:::load_range_botw("Some species", botw, aliases),
    "no BOTW equivalent",
    fixed = TRUE
  )
  expect_null(r2)
})

test_that("load_range_botw handles empty alias data frames gracefully", {
  botw <- make_test_botw()
  on.exit(unlink(botw), add = TRUE)

  empty_aliases <- data.frame(
    ebird_sci_name = character(0),
    botw_sci_name  = character(0),
    stringsAsFactors = FALSE
  )

  r <- ebirdabund:::load_range_botw("Neochmia modesta", botw, empty_aliases)
  expect_equal(nrow(r), 1L)
})
