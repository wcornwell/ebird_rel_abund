# Taxonomy crosswalk: map modeled eBird species to the Reid/Baker NSW bird-risk
# working list (WLAB). The WLAB list is subspecies-level (many trinomials), so
# the join is many-to-one: several WLAB taxa can collapse onto one eBird species.
#
# Join order per species: (1) binomial (genus + species) — the reliable key;
# (2) exact common-name; (3) manual alias override. A species with no match is
# simply not risk-listed. `risk_listed` is TRUE iff >= 1 WLAB taxon matched.

# Region qualifiers in a WLAB subspecies common name that place it clearly
# OUTSIDE NSW. Used only for the wlab_review flag (count NSW-plausible ssp).
# Border-region qualifiers (Southern Queensland, Gippsland, Murray Mallee,
# Lake Eyre Basin, Southern Victoria) are deliberately NOT excluded — they may
# reach into NSW, so we keep them NSW-plausible and err toward flagging. Edit
# this list to tune the flag.
DEFAULT_NON_NSW_PATTERNS <- c(
  "\\bWestern\\b", "South Australian", "\\bTasmanian\\b", "Kangaroo Island",
  "Nullarbor", "Kimberley", "Northern Territory", "New Zealand", "Norfolk",
  "Lord Howe", "Cape York", "Arnhem", "Pilbara", "Top End",
  "Central Queensland Coast", "Western Queensland"
)

# First two whitespace-separated tokens of a scientific name.
binomial <- function(x) {
  vapply(strsplit(trimws(x), "\\s+"), function(p)
    paste(utils::head(p, 2L), collapse = " "), character(1))
}

#' Build the eBird -> WLAB taxonomy crosswalk
#'
#' @param species_common Character vector of modeled species (eBird common names).
#' @param ebird_taxonomy data.frame with `common_name`, `scientific_name`.
#' @param reidbaker data.frame with `WLAB_ScientificName`, `WLAB_CommonName`,
#'   `WLAB_Taxon_ID`.
#' @param aliases Optional data.frame with `ebird_common_name` and
#'   `wlab_scientificname` (the WLAB name to force-join for taxonomy
#'   mismatches). One row per override.
#' @param non_nsw_patterns Character vector of case-insensitive regex tokens in
#'   a WLAB common name that place a subspecies clearly *outside* NSW. Used only
#'   to decide the `wlab_review` flag: a species is flagged when >= 2 of its
#'   WLAB subspecies are NSW-plausible (i.e. match none of these tokens).
#' @return data.frame, one row per modeled species: `common_name`,
#'   `scientific_name`, `wlab_taxon_id`, `wlab_scientificname`,
#'   `n_wlab_subspecies`, `n_nsw_subspecies`, `risk_listed`, `match_type`,
#'   `wlab_review`, `wlab_review_reason`.
#' @export
build_taxonomy_crosswalk <- function(species_common, ebird_taxonomy, reidbaker,
                                     aliases = NULL,
                                     non_nsw_patterns = DEFAULT_NON_NSW_PATTERNS) {
  reidbaker$binom <- binomial(reidbaker$WLAB_ScientificName)

  sci <- ebird_taxonomy$scientific_name[
    match(species_common, ebird_taxonomy$common_name)]
  bn  <- binomial(sci)

  alias_binom <- rep(NA_character_, length(species_common))
  if (!is.null(aliases) && nrow(aliases) > 0) {
    ai <- match(species_common, aliases$ebird_common_name)
    alias_binom <- binomial(aliases$wlab_scientificname[ai])
  }

  collapse <- function(hit) list(
    wlab_taxon_id       = paste(hit$WLAB_Taxon_ID,       collapse = ";"),
    wlab_scientificname = paste(hit$WLAB_ScientificName, collapse = ";"),
    n_wlab_subspecies   = nrow(hit))

  rows <- lapply(seq_along(species_common), function(i) {
    mt <- "none"; hit <- reidbaker[0, ]
    if (!is.na(alias_binom[i])) {                       # 3. alias override
      hit <- reidbaker[reidbaker$binom == alias_binom[i], ]
      if (nrow(hit)) mt <- "alias"
    }
    if (!nrow(hit) && !is.na(bn[i])) {                  # 1. binomial
      hit <- reidbaker[reidbaker$binom == bn[i], ]
      if (nrow(hit)) mt <- "binomial"
    }
    if (!nrow(hit)) {                                   # 2. common name
      hit <- reidbaker[tolower(reidbaker$WLAB_CommonName) ==
                         tolower(species_common[i]), ]
      if (nrow(hit)) mt <- "common"
    }
    cc <- collapse(hit)

    # Flag only genuine lumping where >= 2 subspecies plausibly occur in NSW.
    # A single (NSW) subspecies is treated as the species itself and not
    # flagged. Subspecies whose common name places them clearly outside NSW are
    # not counted. We only FLAG for review — we never model subspecies apart.
    n_ssp <- nrow(hit)
    n_nsw <- if (n_ssp)
      sum(!grepl(paste(non_nsw_patterns, collapse = "|"),
                 hit$WLAB_CommonName, ignore.case = TRUE)) else 0L
    review_reason <- if (n_nsw >= 2L) "multiple_nsw_subspecies" else NA_character_

    data.frame(
      common_name         = species_common[i],
      scientific_name     = sci[i],
      wlab_taxon_id       = if (nrow(hit)) cc$wlab_taxon_id else NA_character_,
      wlab_scientificname = if (nrow(hit)) cc$wlab_scientificname else NA_character_,
      n_wlab_subspecies   = if (nrow(hit)) cc$n_wlab_subspecies else 0L,
      n_nsw_subspecies    = as.integer(n_nsw),
      risk_listed         = nrow(hit) > 0,
      match_type          = mt,
      wlab_review         = !is.na(review_reason),
      wlab_review_reason  = review_reason,
      stringsAsFactors    = FALSE)
  })
  do.call(rbind, rows)
}
