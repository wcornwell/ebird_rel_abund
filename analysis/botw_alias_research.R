#!/usr/bin/env Rscript
# Research BOTW scientific names for the 44 species currently logged as
# range_source = "unmasked". For each one, query BOTW directly with the
# eBird name (expected to fail), then with a fuzzy LIKE on the species
# epithet, then with known historical genus names. Output:
# analysis/botw_alias_research.csv
#
# This is investigative scaffolding for botw_name_aliases.csv. The CSV it
# produces is NOT the alias file — it's a working report that a human
# (Claude in this conversation, or anyone later) curates into conservative
# aliases.

suppressPackageStartupMessages({
  library(sf)
})

BOTW <- "botw_species/BOTW_2025.gpkg"
LOG  <- "batch_nsw_log.csv"

stopifnot(file.exists(BOTW), file.exists(LOG))

log <- read.csv(LOG, stringsAsFactors = FALSE)
unmasked <- log[!is.na(log$range_source) & log$range_source == "unmasked",
                c("common_name", "scientific_name")]
unmasked <- unmasked[order(unmasked$common_name), ]
cat(sprintf("Investigating %d unmasked species against BOTW.\n", nrow(unmasked)))

# --- helpers ----------------------------------------------------------------

botw_query <- function(where_clause) {
  q <- sprintf(
    "SELECT DISTINCT sci_name FROM all_species WHERE %s",
    where_clause
  )
  r <- tryCatch(st_read(BOTW, query = q, quiet = TRUE),
                error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(character(0))
  sort(unique(r$sci_name))
}

epithet <- function(sci_name) {
  parts <- strsplit(sci_name, "\\s+")[[1]]
  if (length(parts) < 2) return(NA_character_)
  parts[2]
}

# For each species: exact, epithet wildcard, "ends with epithet"
investigate <- function(ebird_name) {
  ep <- epithet(ebird_name)
  exact <- botw_query(sprintf("sci_name = '%s'",
                              gsub("'", "''", ebird_name)))
  # Epithet match: any sci_name whose second word matches the epithet OR a
  # close feminine/masculine variant (-um <-> -a, -us <-> -a, etc.).
  # Use a single LIKE pattern so we also catch trinomials.
  variants <- unique(c(
    ep,
    sub("um$", "a",  ep),
    sub("a$",  "um", ep),
    sub("us$", "a",  ep),
    sub("a$",  "us", ep),
    sub("is$", "e",  ep),
    sub("e$",  "is", ep)
  ))
  variants <- variants[!is.na(variants) & nchar(variants) > 2]
  like_clauses <- sprintf("sci_name LIKE '%% %s' OR sci_name LIKE '%% %s %%'",
                          variants, variants)
  epithet_hits <- if (length(variants))
    botw_query(paste(like_clauses, collapse = " OR ")) else character(0)
  # Drop the eBird name itself if it happens to be in there (shouldn't, since
  # these are unmasked).
  epithet_hits <- setdiff(epithet_hits, ebird_name)
  list(
    ebird_name   = ebird_name,
    epithet      = ep,
    exact_hit    = if (length(exact)) exact[1] else NA_character_,
    epithet_hits = paste(epithet_hits, collapse = " | ")
  )
}

# --- run --------------------------------------------------------------------

results <- lapply(unmasked$scientific_name, investigate)
df <- data.frame(
  common_name  = unmasked$common_name,
  ebird_sci    = unmasked$scientific_name,
  epithet      = vapply(results, `[[`, character(1), "epithet"),
  exact_in_botw = vapply(results, `[[`, character(1), "exact_hit"),
  botw_candidates = vapply(results, `[[`, character(1), "epithet_hits"),
  stringsAsFactors = FALSE
)

out <- "analysis/botw_alias_research.csv"
write.csv(df, out, row.names = FALSE)
cat(sprintf("\nWrote %s (%d rows)\n", out, nrow(df)))
cat("\nSummary:\n")
cat(sprintf("  exact match in BOTW : %d\n",
            sum(!is.na(df$exact_in_botw))))
cat(sprintf("  no candidates       : %d\n",
            sum(df$botw_candidates == "" & is.na(df$exact_in_botw))))
cat(sprintf("  >=1 candidate       : %d\n",
            sum(df$botw_candidates != "")))
