# Shared test fixtures used across the test suite.
# Functions prefixed with make_ build synthetic objects on demand so each test
# gets a fresh copy without global state.

make_test_polygon <- function() {
  sf::st_as_sfc(
    sf::st_bbox(c(xmin = 150, ymin = -34, xmax = 151, ymax = -33),
                crs = 4326)
  ) |> sf::st_as_sf()
}

# Synthetic checklist dataframe matching the schema returned by load_ebird().
make_test_checklists <- function(n = 200, n_detections = 60, seed = 42) {
  withr::with_seed(seed, {
    counts <- c(
      sample(1:10, n_detections, replace = TRUE),
      rep(0L, n - n_detections)
    )
    data.frame(
      checklist_id              = paste0("S", seq_len(n)),
      latitude                  = runif(n, -34, -33),
      longitude                 = runif(n, 150, 151),
      observation_date          = as.Date("2020-01-01") + sample(0:364, n, replace = TRUE),
      day_of_year               = sample(1:365, n, replace = TRUE),
      year                      = sample(2018:2022, n, replace = TRUE),
      week                      = sample(1:52, n, replace = TRUE),
      time_observations_started = runif(n, 5, 10),
      duration_minutes          = sample(c(30, 60, 90, 120), n, replace = TRUE),
      effort_distance_km        = runif(n, 0.1, 5),
      number_observers          = sample(1:4, n, replace = TRUE),
      protocol_type             = factor(
        sample(c("Traveling Count", "Stationary Count"), n, replace = TRUE)
      ),
      observation_count         = counts,
      species_observed          = counts > 0,
      lc_trees                  = runif(n, 0, 1),
      lc_grassland              = runif(n, 0, 1),
      lc_shrubs                 = runif(n, 0, 1),
      lc_cropland               = runif(n, 0, 1),
      lc_built                  = runif(n, 0, 1),
      lc_water                  = runif(n, 0, 1),
      elevation                 = runif(n, 0, 500),
      stringsAsFactors          = FALSE
    )
  })
}

# Small synthetic SpatRaster covering the test polygon extent.
make_test_cov_stack <- function(seed = 42) {
  layer_names <- c("lc_trees", "lc_grassland", "lc_shrubs",
                   "lc_cropland", "lc_built", "lc_water", "elevation")
  withr::with_seed(seed, {
    layers <- lapply(layer_names, function(nm) {
      r <- terra::rast(
        nrows = 20, ncols = 20,
        xmin = 149.5, xmax = 151.5,
        ymin = -34.5, ymax = -32.5,
        crs  = "EPSG:4326"
      )
      terra::values(r) <- runif(terra::ncell(r))
      r
    })
  })
  stack <- terra::rast(layers)
  names(stack) <- layer_names
  stack
}

# Fit a minimal GAM on synthetic data for use in prediction tests.
# Cached in the test environment to avoid refitting for every test.
make_test_gam <- function(df = NULL) {
  if (is.null(df)) df <- make_test_checklists()
  suppressMessages(ebirdabund:::fit_gam(df))
}
