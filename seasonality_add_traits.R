# seasonality_add_traits.R
# Additively join species traits (AVONET migration + feeding guild, EltonTraits
# stratum) onto rez_seasonality/seasonality_all.csv, keyed by common_name.
# No rows dropped or filtered — columns only.

SEAS   <- "rez_seasonality/seasonality_all.csv"
TRAITS <- "taxonomy/species_traits.csv"

d <- read.csv(SEAS, stringsAsFactors = FALSE, check.names = FALSE)
t <- read.csv(TRAITS, stringsAsFactors = FALSE)

trait_cols <- c("migration", "feeding_guild", "trophic_level",
                "primary_lifestyle", "foraging_stratum")
d <- d[, setdiff(names(d), trait_cols)]         # idempotent: drop if re-run
i <- match(d$common_name, t$common_name)
add <- t[i, trait_cols]

# Insert the trait block right after the taxonomy columns (after
# wlab_review_reason) if present, else append at the end.
anchor <- match("wlab_review_reason", names(d))
if (is.na(anchor)) {
  out <- cbind(d, add)
} else {
  out <- cbind(d[, 1:anchor, drop = FALSE], add,
               d[, (anchor + 1):ncol(d), drop = FALSE])
}

write.csv(out, SEAS, row.names = FALSE)
message(sprintf("Added %s to %s (%d rows).",
                paste(trait_cols, collapse = ", "), SEAS, nrow(out)))
message(sprintf("  migration populated: %d / %d rows",
                sum(!is.na(out$migration)), nrow(out)))
