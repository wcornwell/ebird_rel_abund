suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(data.table)
  library(parallel)
})

ebd_path   <- paste0("ebirdabund/raw_data/",
                     "ebd_AU-NSW_unv_smp_relFeb-2026/",
                     "ebd_AU-NSW_unv_smp_relFeb-2026.txt")
samp_path  <- paste0("ebirdabund/raw_data/",
                     "ebd_AU-NSW_unv_smp_relFeb-2026/",
                     "ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt")
cache_dir  <- "ebirdabund_cache"

# ── NSW boundary ──────────────────────────────────────────────────────────────
aus <- geodata::gadm(country = "AUS", level = 1, path = cache_dir)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])

# ── Build new covariate stack (includes pop_density) ─────────────────────────
message("Building covariate stack with population density...")
cov <- prepare_covariates(nsw, cache_dir = cache_dir)
message("Layers: ", paste(names(cov), collapse = ", "))

# ── Stratified species sample ─────────────────────────────────────────────────
# Pick ~40 species spread evenly across reporting-rate deciles so the test
# covers both common and uncommon species.
species_df <- read.csv("nsw_species_list.csv", stringsAsFactors = FALSE)
species_df <- species_df[species_df$reporting_rate >= 0.005, ]

set.seed(42)
species_df$decile <- cut(
  species_df$reporting_rate,
  breaks   = quantile(species_df$reporting_rate, probs = seq(0, 1, 0.1)),
  include.lowest = TRUE,
  labels   = FALSE
)
sample_sp <- do.call(rbind, lapply(1:10, function(d) {
  pool <- species_df[species_df$decile == d, ]
  pool[sample(nrow(pool), min(4, nrow(pool))), ]
}))
message(sprintf(
  "Testing %d species across %d deciles",
  nrow(sample_sp), length(unique(sample_sp$decile))
))

# ── Fit one species, return F-stats for all smooth terms ─────────────────────
fit_one <- function(sp, cov_stack) {
  tryCatch({
    mf  <- fit_species_model(
      polygon        = nsw,
      ebird_zip      = ebd_path,
      sampling_txt   = samp_path,
      species        = sp,
      cov_stack      = cov_stack,
      cache_dir      = cache_dir,
      hex_spacing_km = 5
    )
    s   <- summary(mf$model)
    tbl <- as.data.frame(s$s.table)
    tbl$term     <- rownames(tbl)
    tbl$species  <- sp
    tbl$dev_expl <- s$dev.expl
    tbl
  }, error = function(e) {
    message(sprintf("  [FAILED] %s: %s", sp, conditionMessage(e)))
    NULL
  })
}

# ── Parallel run ──────────────────────────────────────────────────────────────
n_cores <- max(1L, parallel::detectCores() - 1L)
workers <- min(n_cores, nrow(sample_sp))
wrapped <- terra::wrap(cov)
pkg_src <- normalizePath("ebirdabund", mustWork = FALSE)

message(sprintf(
  "\nFitting %d species on %d cores...", nrow(sample_sp), workers
))
cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterExport(
  cl,
  c("pkg_src", "nsw", "ebd_path", "samp_path", "cache_dir",
    "wrapped", "fit_one"),
  envir = environment()
)
parallel::clusterEvalQ(cl, {
  devtools::load_all(pkg_src, quiet = TRUE)
  suppressPackageStartupMessages(library(terra))
})

results <- parallel::parLapplyLB(
  cl,
  sample_sp$common_name,
  function(sp) {
    cov_local <- terra::unwrap(wrapped)
    fit_one(sp, cov_local)
  }
)

df <- rbindlist(Filter(Negate(is.null), results))
names(df)[names(df) == "p-value"] <- "p_value"
names(df)[names(df) == "F"]       <- "f_stat"

# ── Aggregate: habitat covariates only ───────────────────────────────────────
hab <- df[grepl("^s[(](lc_|elevation|precip_|temp_|pop_)", df$term)]

summary_tbl <- hab[, .(
  n_models   = .N,
  median_f   = round(median(f_stat, na.rm = TRUE), 1),
  q10_f      = round(quantile(f_stat, 0.10, na.rm = TRUE), 1),
  q90_f      = round(quantile(f_stat, 0.90, na.rm = TRUE), 1),
  median_edf = round(median(edf, na.rm = TRUE), 2),
  pct_sig    = round(100 * mean(p_value < 0.05, na.rm = TRUE), 1),
  pct_shrunk = round(100 * mean(edf < 1.05, na.rm = TRUE), 1)
), by = term][order(median_f)]

cat(
  "\n\nHabitat covariate importance incl. pop_density",
  sprintf("(%d species)\n", length(unique(df$species)))
)
cat(strrep("=", 85), "\n")
print(summary_tbl, row.names = FALSE)
cat("\npct_shrunk = % of species where smooth was penalised to ~zero\n")

saveRDS(summary_tbl, "pop_density_test_result.rds")
message("\nResult saved to pop_density_test_result.rds")
