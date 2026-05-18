# Tests for covariate-variant testing helpers (R/test_covariate.R) and the
# fold_ids extension on evaluate_model_cv().

# ── validate_change ─────────────────────────────────────────────────────────
test_that("validate_change accepts swap and add", {
  expect_silent(ebirdabund:::validate_change(list(swap = "tree_height")))
  expect_silent(ebirdabund:::validate_change(list(add  = "new_col")))
})

test_that("validate_change rejects bad shapes", {
  expect_error(ebirdabund:::validate_change(list()),
               "length-1 list")
  expect_error(ebirdabund:::validate_change(list(replace = "x")),
               "named 'swap' or 'add'")
  expect_error(ebirdabund:::validate_change(list(swap = c("a", "b"))),
               "single non-empty column name")
  expect_error(ebirdabund:::validate_change(list(swap = "")),
               "single non-empty column name")
})


# ── extract_variant_covariate ───────────────────────────────────────────────
test_that("extract_variant_covariate writes a cache and returns aligned rows", {
  tmp_cache <- withr::local_tempdir()
  cov_stack <- make_test_cov_stack()
  # Build a tiny variant raster aligned with the test polygon extent
  variant_r <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 149.5, xmax = 151.5,
    ymin = -34.5, ymax = -32.5,
    crs = "EPSG:4326"
  )
  terra::values(variant_r) <- runif(terra::ncell(variant_r), 0, 30)
  names(variant_r) <- "tree_height_alt"
  variant_path <- file.path(tmp_cache, "alt.tif")
  terra::writeRaster(variant_r, variant_path)

  master <- data.frame(
    checklist_id = paste0("S", 1:5),
    longitude    = c(150.1, 150.4, 150.6, 150.8, 150.95),
    latitude     = c(-33.2, -33.4, -33.6, -33.8, -33.9),
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(ebirdabund:::extract_variant_covariate(
    variant_source    = variant_path,
    variant_name      = "alt_height",
    master            = master,
    variant_cache_dir = tmp_cache
  ))

  expect_equal(nrow(out), nrow(master))
  expect_setequal(out$checklist_id, master$checklist_id)
  expect_true("variant_value" %in% names(out))
  expect_equal(attr(out, "variant_layer_name"), "tree_height_alt")
  expect_true(file.exists(file.path(tmp_cache, "alt_height.rds")))
})

test_that("extract_variant_covariate uses cache on a second call", {
  tmp_cache <- withr::local_tempdir()
  variant_r <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 149.5, xmax = 151.5, ymin = -34.5, ymax = -32.5,
    crs  = "EPSG:4326"
  )
  terra::values(variant_r) <- 1
  names(variant_r) <- "h"
  variant_path <- file.path(tmp_cache, "v.tif")
  terra::writeRaster(variant_r, variant_path)

  master <- data.frame(
    checklist_id = paste0("S", 1:3),
    longitude    = c(150.1, 150.5, 150.9),
    latitude     = c(-33.2, -33.5, -33.8),
    stringsAsFactors = FALSE
  )

  first <- suppressMessages(ebirdabund:::extract_variant_covariate(
    variant_source    = variant_path,
    variant_name      = "v",
    master            = master,
    variant_cache_dir = tmp_cache
  ))
  msgs <- capture.output(
    second <- ebirdabund:::extract_variant_covariate(
      variant_source    = variant_path,
      variant_name      = "v",
      master            = master,
      variant_cache_dir = tmp_cache
    ),
    type = "message"
  )
  expect_equal(first, second)
  expect_true(any(grepl("loaded from cache", msgs)))
})

test_that("extract_variant_covariate accepts a function source", {
  tmp_cache <- withr::local_tempdir()
  master <- data.frame(
    checklist_id = paste0("S", 1:4),
    longitude    = c(150.1, 150.4, 150.6, 150.9),
    latitude     = c(-33.2, -33.4, -33.6, -33.9),
    stringsAsFactors = FALSE
  )
  polygon <- make_test_polygon()
  src <- function(bbox, template) {
    r <- terra::rast(
      nrows = 5, ncols = 5,
      xmin = bbox[1], xmax = bbox[3],
      ymin = bbox[2], ymax = bbox[4],
      crs = "EPSG:4326"
    )
    terra::values(r) <- seq_len(terra::ncell(r))
    names(r) <- "derived_col"
    r
  }
  out <- suppressMessages(ebirdabund:::extract_variant_covariate(
    variant_source    = src,
    variant_name      = "func_test",
    master            = master,
    variant_cache_dir = tmp_cache,
    polygon           = polygon
  ))
  expect_true("variant_value" %in% names(out))
  expect_equal(attr(out, "variant_layer_name"), "derived_col")
  expect_equal(nrow(out), nrow(master))
})


# ── apply_variant ───────────────────────────────────────────────────────────
test_that("apply_variant swap replaces the target column", {
  df <- make_test_checklists(n = 20, n_detections = 8)
  df$checklist_id <- paste0("S", seq_len(nrow(df)))
  df$tree_height  <- runif(nrow(df), 0, 10)
  variant <- data.frame(
    checklist_id  = df$checklist_id,
    variant_value = runif(nrow(df), 100, 110),  # very different range
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(ebirdabund:::apply_variant(
    df, variant, list(swap = "tree_height")
  ))
  expect_true(all(out$tree_height >= 100))
  expect_false("variant_value" %in% names(out))
  expect_equal(nrow(out), nrow(df))
})

test_that("apply_variant swap survives same-named variant raster layer", {
  # Regression for the bug where variant layer name == swap target caused
  # merge() to produce target.x / target.y and the new column to vanish.
  df <- make_test_checklists(n = 12, n_detections = 4)
  df$checklist_id <- paste0("S", seq_len(nrow(df)))
  df$tree_height  <- 1   # baseline value
  variant <- data.frame(
    checklist_id  = df$checklist_id,
    variant_value = 99,
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(ebirdabund:::apply_variant(
    df, variant, list(swap = "tree_height")
  ))
  expect_true(all(out$tree_height == 99))
  expect_equal(nrow(out), nrow(df))
})

test_that("apply_variant add appends a new column without dropping rows", {
  df <- make_test_checklists(n = 15, n_detections = 6)
  df$checklist_id <- paste0("S", seq_len(nrow(df)))
  variant <- data.frame(
    checklist_id  = df$checklist_id,
    variant_value = runif(nrow(df)),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(ebirdabund:::apply_variant(
    df, variant, list(add = "new_layer")
  ))
  expect_true("new_layer" %in% names(out))
  expect_equal(nrow(out), nrow(df))
})

test_that("apply_variant add refuses to clobber an existing column", {
  df <- make_test_checklists(n = 10, n_detections = 4)
  df$checklist_id <- paste0("S", seq_len(nrow(df)))
  df$elevation    <- 100
  variant <- data.frame(
    checklist_id  = df$checklist_id,
    variant_value = runif(nrow(df)),
    stringsAsFactors = FALSE
  )
  expect_error(
    suppressMessages(ebirdabund:::apply_variant(
      df, variant, list(add = "elevation")
    )),
    "already exists"
  )
})

test_that("apply_variant drops rows with NA in the variant column", {
  df <- make_test_checklists(n = 10, n_detections = 4)
  df$checklist_id <- paste0("S", seq_len(nrow(df)))
  df$tree_height  <- runif(nrow(df))
  variant <- data.frame(
    checklist_id  = df$checklist_id,
    variant_value = c(rep(NA_real_, 3), runif(7)),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(ebirdabund:::apply_variant(
    df, variant, list(swap = "tree_height")
  ))
  expect_equal(nrow(out), 7L)
})


# ── assemble_metric_row / summarise_deltas ──────────────────────────────────
test_that("assemble_metric_row produces base/var/delta columns", {
  cv_base <- list(summary = c(spearman_r = 0.5, pearson_r = 0.4,
                              mae = 1.0, rmse = 1.5,
                              holdout_dev_expl = 0.2))
  cv_var  <- list(summary = c(spearman_r = 0.6, pearson_r = 0.5,
                              mae = 0.9, rmse = 1.4,
                              holdout_dev_expl = 0.3))
  row <- ebirdabund:::assemble_metric_row("Foo", cv_base, cv_var, n_pos = 100L)
  expect_equal(row$spearman_r_base,  0.5)
  expect_equal(row$spearman_r_var,   0.6)
  expect_equal(row$spearman_r_delta, 0.1, tolerance = 1e-12)
  expect_equal(row$rmse_delta, -0.1, tolerance = 1e-12)
})

test_that("summarise_deltas reports win/tie/loss with correct direction", {
  per_species <- data.frame(
    species = c("A", "B", "C", "D"),
    spearman_r_base  = c(0.5, 0.5, 0.5, 0.5),
    spearman_r_var   = c(0.6, 0.7, 0.5, 0.4),
    pearson_r_base   = c(0.4, 0.4, 0.4, 0.4),
    pearson_r_var    = c(0.4, 0.4, 0.4, 0.4),
    mae_base         = c(1.0, 1.0, 1.0, 1.0),
    mae_var          = c(0.8, 0.9, 1.0, 1.2),
    rmse_base        = c(1.5, 1.5, 1.5, 1.5),
    rmse_var         = c(1.4, 1.4, 1.4, 1.4),
    holdout_dev_expl_base = c(0.2, 0.2, 0.2, 0.2),
    holdout_dev_expl_var  = c(0.3, 0.3, 0.3, 0.3),
    stringsAsFactors = FALSE
  )
  # Add the *_delta columns the function expects to NOT exist (it computes
  # them itself)
  out <- ebirdabund:::summarise_deltas(per_species)
  sp <- out[out$metric == "spearman_r", ]
  expect_equal(sp$wins,   2)  # B, A
  expect_equal(sp$losses, 1)  # D
  expect_equal(sp$ties,   1)  # C

  mae <- out[out$metric == "mae", ]
  expect_equal(mae$wins,   2)  # A, B (lower is better)
  expect_equal(mae$losses, 1)  # D
  expect_equal(mae$ties,   1)
})


# ── evaluate_model_cv fold_ids extension ────────────────────────────────────
test_that("evaluate_model_cv with explicit fold_ids reproduces seeded result", {
  df <- make_test_checklists(n = 80, n_detections = 30)
  set.seed(42L)
  expected_fold <- sample(rep_len(seq_len(5L), nrow(df)))

  # Reuse the same seed inside the function (default seed=42)
  out_default <- suppressMessages(suppressWarnings(
    ebirdabund:::evaluate_model_cv(df, k = 5L, seed = 42L)
  ))
  out_explicit <- suppressMessages(suppressWarnings(
    ebirdabund:::evaluate_model_cv(df, k = 5L, fold_ids = expected_fold)
  ))
  expect_equal(out_default$summary, out_explicit$summary, tolerance = 1e-8)
})

test_that("evaluate_model_cv rejects wrong-length fold_ids", {
  df <- make_test_checklists(n = 30, n_detections = 10)
  expect_error(
    suppressMessages(ebirdabund:::evaluate_model_cv(
      df, k = 5L, fold_ids = seq_len(5L)
    )),
    "length nrow"
  )
})

test_that("evaluate_model_cv rejects fold_ids outside 1..k", {
  df <- make_test_checklists(n = 30, n_detections = 10)
  bad <- rep(c(1L, 6L), length.out = nrow(df))
  expect_error(
    suppressMessages(ebirdabund:::evaluate_model_cv(
      df, k = 5L, fold_ids = bad
    )),
    "integers in 1..k"
  )
})
