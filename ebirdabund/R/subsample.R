# Spatiotemporal subsampling: keep at most one checklist per
# hexagonal cell × year × week × detection-status combination.
# This reduces spatial and temporal clustering bias before model fitting.
#
# spacing_km: approximate diameter of hexagonal cells in km.
subsample_hex <- function(df, spacing_km = 5) {
  n_before <- nrow(df)
  dgg <- dggridR::dgconstruct(spacing = spacing_km)

  df <- df |>
    dplyr::mutate(
      hex_cell = dggridR::dgGEO_to_SEQNUM(dgg, .data$longitude, .data$latitude)$seqnum
    ) |>
    dplyr::group_by(.data$year, .data$week, .data$hex_cell, .data$species_observed) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-"hex_cell")

  message(sprintf(
    "  %d → %d checklists retained (%.0f%% removed by spatiotemporal subsampling).",
    n_before, nrow(df), 100 * (n_before - nrow(df)) / n_before
  ))

  df
}
