suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
  library(geodata)
  if (!"dplyr" %in% .packages()) library(dplyr)
})

EBD      <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt"
SAMP     <- "ebirdabund/raw_data/ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026_sampling.txt"
CACHE    <- "ebirdabund_cache"
OUT_CSV  <- "nsw_species_list.csv"
OUT_PLOT <- "nsw_species_rarity_cutoff.png"

# ── 1. Complete checklists from the sampling file ─────────────────────────────
message("Reading sampling events file (complete checklists only)...")
samp <- as.data.frame(fread(
  SAMP, sep = "\t", quote = "", showProgress = FALSE, na.strings = "",
  select = c("SAMPLING EVENT IDENTIFIER", "ALL SPECIES REPORTED",
             "LATITUDE", "LONGITUDE",
             "DURATION MINUTES", "EFFORT DISTANCE KM", "NUMBER OBSERVERS")
))
samp <- samp[samp[["ALL SPECIES REPORTED"]] == 1, ]

# Apply standard eBird effort filters to match the modelling pipeline
samp <- samp[
  !is.na(samp[["DURATION MINUTES"]]) &
  samp[["DURATION MINUTES"]] >= 5 &
  samp[["DURATION MINUTES"]] <= 300 &
  (is.na(samp[["EFFORT DISTANCE KM"]]) | samp[["EFFORT DISTANCE KM"]] <= 10) &
  !is.na(samp[["NUMBER OBSERVERS"]]) &
  samp[["NUMBER OBSERVERS"]] <= 10, ]

message(sprintf("  %d complete checklists after effort filtering", nrow(samp)))

# ── 1b. Remove offshore checklists ────────────────────────────────────────────
# Checklists outside the NSW land polygon have no covariate coverage and are
# unmodelable (pelagic seabird surveys, offshore ferry routes, etc.).
# Filtering them here ensures reporting_rate reflects land-based detections only.
message("Loading NSW land boundary to filter offshore checklists...")
aus <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
nsw <- sf::st_transform(nsw, 4326)

samp_pts <- sf::st_as_sf(
  samp[!is.na(samp$LATITUDE) & !is.na(samp$LONGITUDE), ],
  coords = c("LONGITUDE", "LATITUDE"),
  crs    = 4326,
  remove = FALSE
)
on_land <- lengths(sf::st_within(samp_pts, nsw)) > 0
samp    <- samp_pts[on_land, ]
sf::st_geometry(samp) <- NULL   # back to plain data.frame

n_offshore <- sum(!on_land)
message(sprintf("  Removed %d offshore/out-of-polygon checklists (%d remaining)",
                n_offshore, nrow(samp)))

valid_ids  <- samp[["SAMPLING EVENT IDENTIFIER"]]
n_complete <- length(valid_ids)
message(sprintf("  %d land-based complete checklists", n_complete))

# ── 2. All species observations in those checklists ───────────────────────────
message("Reading EBD (columns: COMMON NAME, OBSERVATION COUNT, SAMPLING EVENT IDENTIFIER)...")
ebd <- as.data.frame(fread(
  EBD,
  select       = c(6L, 11L, 35L),   # positions stable in Feb-2026 EBD
  sep          = "\t", quote        = "",
  showProgress = TRUE, na.strings   = ""
))
names(ebd) <- c("common_name", "observation_count", "checklist_id")

# Keep only valid complete checklists; drop X counts; exclude non-species:
#   slash taxa  "A/B"
#   spuh        "...sp."           (also matches "...sp. (qualifier)" forms
#                                   like "Accipitrine hawk sp. (former
#                                   Accipiter sp.)" — anchoring to $ misses
#                                   those because of the trailing
#                                   parenthetical)
#   hybrids     "A x B (hybrid)"
is_true_species <- function(name) {
  !grepl("/|\\bsp\\.| x .*\\(hybrid\\)", name, perl = TRUE)
}

ebd <- ebd[
  ebd[["checklist_id"]] %in% valid_ids &
  !is.na(ebd[["observation_count"]]) &
  ebd[["observation_count"]] != "X" &
  is_true_species(ebd[["common_name"]]), ]

ebd[["observation_count"]] <- as.integer(ebd[["observation_count"]])
message(sprintf("  %d observation records retained", nrow(ebd)))

# ── 3. Per-species metrics ────────────────────────────────────────────────────
message("Computing per-species metrics...")
species_stats <- ebd |>
  group_by(common_name) |>
  summarise(
    n_detections     = n(),
    mean_count       = round(mean(observation_count, na.rm = TRUE), 3),
    median_count     = median(observation_count, na.rm = TRUE),
    max_count        = max(observation_count, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    reporting_rate   = round(n_detections / n_complete * 100, 4),
    n_checklists     = n_complete
  ) |>
  arrange(desc(reporting_rate))

message(sprintf("  %d species detected in NSW", nrow(species_stats)))

# ── 4. Evidence-based rarity cutoff analysis ──────────────────────────────────
# Three candidate thresholds commonly used in the literature:
#   T1: 1% reporting rate  — detectable, but relatively rare
#   T2: 0.1% reporting rate — the standard eBird "sparse data" floor
#   T3: natural break at the knee of the cumulative rank-abundance curve

rr <- sort(species_stats$reporting_rate, decreasing = TRUE)

t1_n <- sum(rr >= 1)
t2_n <- sum(rr >= 0.1)

# Knee: maximise distance from the diagonal of the log-transformed rank-abundance curve
log_rr <- log10(rr + 1e-6)
x_norm <- seq(0, 1, length.out = length(log_rr))
y_norm <- (log_rr - min(log_rr)) / diff(range(log_rr))
dist_diag <- abs(y_norm - (1 - x_norm))   # distance from y=1-x line
knee_idx  <- which.max(dist_diag)
t3_n      <- knee_idx
t3_rr     <- rr[knee_idx]

message("\n── Rarity cutoff candidates ─────────────────────────────────────")
message(sprintf("  >= 1%%   reporting rate  →  %d species", t1_n))
message(sprintf("  >= 0.1%% reporting rate  →  %d species", t2_n))
message(sprintf("  Knee of rank-abundance curve (%.3f%%) → %d species", t3_rr, t3_n))

# ── 5. Add rarity classification column ──────────────────────────────────────
species_stats <- species_stats |>
  mutate(
    category = case_when(
      reporting_rate >= 1   ~ "regular (>=1%)",
      reporting_rate >= 0.1 ~ "uncommon (0.1-1%)",
      TRUE                  ~ "rare (<0.1%)"
    )
  )

# ── 6. Write CSV ──────────────────────────────────────────────────────────────
write.csv(species_stats, OUT_CSV, row.names = FALSE)
message(sprintf("\nSpecies list written to: %s", OUT_CSV))

# ── 7. Diagnostic plots ───────────────────────────────────────────────────────
message("Generating rarity cutoff diagnostic plot...")

# Rank-abundance curve (log scale) with candidate thresholds annotated
rank_df <- data.frame(
  rank            = seq_along(rr),
  reporting_rate  = rr
)

p1 <- ggplot(rank_df, aes(rank, reporting_rate)) +
  geom_line(colour = "#2c7bb6", linewidth = 0.7) +
  geom_hline(yintercept = 1,   linetype = "dashed", colour = "#d7191c", linewidth = 0.7) +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "#fdae61", linewidth = 0.7) +
  geom_vline(xintercept = t3_n, linetype = "dotted", colour = "#1a9641", linewidth = 0.8) +
  annotate("text", x = t1_n + 5, y = 1.3,
           label = sprintf("1%% (n=%d)", t1_n),
           colour = "#d7191c", hjust = 0, size = 3) +
  annotate("text", x = t2_n + 5, y = 0.14,
           label = sprintf("0.1%% (n=%d)", t2_n),
           colour = "#a65400", hjust = 0, size = 3) +
  annotate("text", x = t3_n + 5, y = max(rr) * 0.6,
           label = sprintf("Knee %.3f%% (n=%d)", t3_rr, t3_n),
           colour = "#1a9641", hjust = 0, size = 3) +
  scale_y_log10(
    labels = function(x) paste0(x, "%"),
    breaks = c(0.01, 0.1, 1, 10, 100)
  ) +
  labs(
    title    = "NSW Bird Species Rank-Abundance Curve",
    subtitle = sprintf(
      "eBird Feb-2026 release  |  %d complete checklists  |  %d species detected",
      n_complete, nrow(species_stats)
    ),
    x = "Species rank (most to least common)",
    y = "Reporting rate (% of complete checklists)"
  ) +
  theme_minimal(base_size = 12)

# Histogram of reporting rates in the low end (zoom on the rare tail)
p2 <- species_stats |>
  filter(reporting_rate < 5) |>
  ggplot(aes(reporting_rate, fill = category)) +
  geom_histogram(bins = 80, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 1,   linetype = "dashed", colour = "#d7191c") +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "#fdae61") +
  geom_vline(xintercept = t3_rr, linetype = "dotted", colour = "#1a9641") +
  scale_fill_manual(values = c(
    "regular (>=1%)"   = "#2c7bb6",
    "uncommon (0.1-1%)" = "#fdae61",
    "rare (<0.1%)"     = "#d7191c"
  )) +
  labs(
    title = "Distribution of reporting rates (< 5% shown)",
    x = "Reporting rate (%)", y = "Number of species", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

combined <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1.4, 1))
cowplot::save_plot(OUT_PLOT, combined, base_height = 10, base_width = 10)
message(sprintf("Diagnostic plot saved: %s", OUT_PLOT))

# ── 8. Console summary ────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════════════════════════════\n")
cat("NSW SPECIES LIST SUMMARY\n")
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("Total checklists analysed : %d\n", n_complete))
cat(sprintf("Total species detected    : %d\n", nrow(species_stats)))
cat(sprintf("\nCandidate rarity thresholds:\n"))
cat(sprintf("  ≥ 1%%   reporting rate → %d  species  (recommended default)\n", t1_n))
cat(sprintf("  ≥ 0.1%% reporting rate → %d species\n", t2_n))
cat(sprintf("  Knee of curve (≥%.3f%%) → %d species\n", t3_rr, t3_n))
cat("\nTop 20 most-reported species:\n")
print(head(species_stats[, c("common_name", "reporting_rate", "mean_count", "category")], 20),
      row.names = FALSE)
cat(sprintf("\nFull table in: %s\n", OUT_CSV))
cat(sprintf("Plot in      : %s\n", OUT_PLOT))
