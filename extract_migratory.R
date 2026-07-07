# extract_migratory.R
# Derive a per-species migratory flag from the BirdLife BOTW geopackage.
#
# BOTW `all_species` codes each range polygon's `seasonal`:
#   1 resident | 2 breeding | 3 non-breeding | 4 passage | 5 uncertain.
# A species mapped with any breeding/non-breeding/passage range (2/3/4) is
# treated as migratory; resident-only (codes in {1,5}) is non-migratory.
# Same presence/origin filter as load_range_botw (extant, established).
#
# CAVEAT: BirdLife maps many intra-Australian partial migrants / nomads as
# resident (seasonal=1), so this UNDER-detects Australian partial migrants
# (e.g. Flame Robin). Conservative: "migrant" is high-precision, "resident"
# may include partial migrants. Output feeds the seasonality report as a
# filter/annotation, not a hard truth.

suppressPackageStartupMessages({ library(sf); library(yaml) })

cfg   <- yaml::read_yaml("config.yaml")
BOTW  <- cfg$botw_path
TAX   <- cfg$taxonomy_file
SPL   <- cfg$species_list_file
ALIAS <- "botw_name_aliases.csv"
OUT   <- "taxonomy/migratory_status.csv"

tax <- read.csv(TAX, stringsAsFactors = FALSE)
spl <- read.csv(SPL, stringsAsFactors = FALSE)
ali <- if (file.exists(ALIAS)) read.csv(ALIAS, stringsAsFactors = FALSE) else
  data.frame(ebird_sci_name = character(), botw_sci_name = character())

# One attribute-only query: seasonal codes for all established populations.
message("Querying BOTW seasonal codes (attribute-only)...")
q <- "SELECT sci_name, seasonal FROM all_species
      WHERE presence IN (1,2) AND origin IN (1,2,3,6)"
codes <- sf::st_read(BOTW, query = q, quiet = TRUE)
codes <- sf::st_drop_geometry(codes)
message(sprintf("  %d range rows across %d BOTW names.",
                nrow(codes), length(unique(codes$sci_name))))

# sci_name -> sorted distinct seasonal codes
by_name <- tapply(codes$seasonal, codes$sci_name,
                  function(x) sort(unique(as.integer(x))))

classify <- function(sci) {
  # apply eBird->BOTW alias if present
  a <- ali$botw_sci_name[!is.na(ali$ebird_sci_name) & ali$ebird_sci_name == sci]
  qn <- if (length(a) && nzchar(a[1])) a[1] else sci
  cc <- by_name[[qn]]
  if (is.null(cc)) return(c(NA_character_, "unknown"))
  status <- if (any(cc %in% c(2L, 3L, 4L))) "migrant" else "resident"
  c(paste(cc, collapse = ";"), status)
}

sci <- tax$scientific_name[match(spl$common_name, tax$common_name)]
cl  <- t(vapply(sci, classify, character(2)))
out <- data.frame(
  common_name        = spl$common_name,
  scientific_name    = sci,
  botw_seasonal_codes = cl[, 1],
  migratory_status   = cl[, 2],
  stringsAsFactors   = FALSE)

write.csv(out, OUT, row.names = FALSE)
message(sprintf("\nWrote %s (%d species)", OUT, nrow(out)))
print(table(out$migratory_status, useNA = "ifany"))
