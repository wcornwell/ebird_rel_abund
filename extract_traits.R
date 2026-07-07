# extract_traits.R
# Per-species trait table combining:
#   - migration        : AVONET (eBird taxonomy) Migration — Sedentary /
#                        Partial migrant / Migratory. Captures partial &
#                        nomadic movers (e.g. Red Wattlebird = Partial) that
#                        the BirdLife range-polygon seasonal codes miss.
#   - feeding_guild    : AVONET Trophic.Niche (eBird-keyed)
#   - trophic_level, primary_lifestyle : AVONET
#   - foraging_stratum : EltonTraits 1.0 dominant vertical stratum (extra)
#   - botw_seasonal_codes : BirdLife BOTW seasonal codes (audit / cross-check)
# Output: taxonomy/species_traits.csv
#
# AVONET is the authority for migration & guild because it is eBird-keyed and
# expert-scored (aligned with HBW/Birds of the World movement accounts).

suppressPackageStartupMessages({ library(readxl); library(yaml) })

cfg   <- yaml::read_yaml("config.yaml")
AVO   <- "taxonomy/traits/AVONET.xlsx"
ELTON <- "taxonomy/traits/BirdFuncDat.txt"
MIGR  <- "taxonomy/migratory_status.csv"          # BOTW seasonal (audit)
OUT   <- "taxonomy/species_traits.csv"

tax <- read.csv(cfg$taxonomy_file, stringsAsFactors = FALSE)
spl <- read.csv(cfg$species_list_file, stringsAsFactors = FALSE)
sci <- tax$scientific_name[match(spl$common_name, tax$common_name)]

binomial <- function(x) vapply(strsplit(trimws(x), "\\s+"),
                               function(p) paste(utils::head(p, 2L), collapse = " "),
                               character(1))

# ── AVONET (eBird sheet) ──────────────────────────────────────────────────────
avo <- as.data.frame(read_excel(AVO, sheet = "AVONET2_eBird"))
avo$binom <- binomial(avo$Species2)
mig_lab <- c("1" = "Sedentary", "2" = "Partial migrant", "3" = "Migratory")

# Apply eBird->AVONET name aliases for taxonomy reshuffles/splits since AVONET's
# 2021 eBird taxonomy (e.g. Anarhynchus->Charadrius, Heteroscenes->Cacomantis).
ALIAS <- "taxonomy/avonet_name_aliases.csv"
sci_q <- sci
if (file.exists(ALIAS)) {
  av_alias <- read.csv(ALIAS, stringsAsFactors = FALSE)
  am <- match(sci, av_alias$ebird_sci_name)
  sci_q[!is.na(am)] <- av_alias$avonet_species2[am[!is.na(am)]]
}

# join eBird scientific -> Species2, then binomial fallback
ai <- match(sci_q, avo$Species2)
miss <- is.na(ai); ai[miss] <- match(binomial(sci_q)[miss], avo$binom)

migration        <- mig_lab[as.character(avo$Migration[ai])]
feeding_guild    <- avo$Trophic.Niche[ai]
trophic_level    <- avo$Trophic.Level[ai]
primary_lifestyle <- avo$Primary.Lifestyle[ai]

# ── EltonTraits: dominant foraging stratum (extra vertical detail) ─────────────
elton <- read.delim(ELTON, stringsAsFactors = FALSE, fileEncoding = "latin1",
                    check.names = FALSE)
elton$binom <- binomial(elton$Scientific)
strat_cols <- c(`ForStrat-watbelowsurf` = "below-water",
                `ForStrat-wataroundsurf` = "around-water",
                `ForStrat-ground` = "ground", `ForStrat-understory` = "understory",
                `ForStrat-midhigh` = "mid-high", `ForStrat-canopy` = "canopy",
                `ForStrat-aerial` = "aerial")
ei <- match(binomial(sci), elton$binom)
ei[is.na(ei)] <- match(tolower(spl$common_name)[is.na(ei)], tolower(elton$English))
foraging_stratum <- vapply(seq_along(ei), function(k) {
  i <- ei[k]; if (is.na(i)) return(NA_character_)
  if (isTRUE(as.integer(elton$PelagicSpecialist[i]) == 1L)) return("pelagic")
  vals <- suppressWarnings(as.numeric(elton[i, names(strat_cols)]))
  if (all(is.na(vals)) || max(vals, na.rm = TRUE) == 0) return(NA_character_)
  unname(strat_cols[which.max(vals)])
}, character(1))

# ── BIRDBASE (2024 taxonomy): IUCN status + habitat breadth ───────────────────
bb   <- as.data.frame(read_excel("taxonomy/traits/BIRDBASE.xlsx", sheet = 1, skip = 1))
bkey <- "eBird/Clements (V2024b)"
bi <- match(sci, bb[[bkey]])
bmiss <- is.na(bi); bi[bmiss] <- match(binomial(sci)[bmiss], binomial(bb[[bkey]]))
iucn_status     <- bb[["2024 IUCN Red List category"]][bi]
habitat_breadth <- suppressWarnings(as.numeric(bb[["HB"]][bi]))

# ── Garnett Australian Bird Data: EPBC/NSW status + national movement ──────────
# Australia-specific & subspecies-level; aggregated to species (most-threatened
# status, union of movement flags). Fills natives that reidbaker/WLAB omits.
gar  <- read.csv("taxonomy/traits/Australian_Bird_Data_v1.csv",
                 stringsAsFactors = FALSE, check.names = FALSE)
gcol <- function(p) grep(p, names(gar), value = TRUE)[1]
gar$binom <- paste(gar[[gcol("Genus_name")]], gar[[gcol("Species_name")]])
gcommon <- gcol("Taxon_common_name"); gepbc <- gcol("EPBC_Status_July")
gnsw <- gcol("NSW_status_2015")
mvcols <- grep("National_movement", names(gar), value = TRUE)
mvlab  <- c("local dispersal", "partial migrant", "total migrant",
            "nomadic", "irruptive")
sev_ord <- c("EX (presumed)"=5,"EX"=5,"CR"=4,"EN (population)"=3,"EN"=3,
             "VU"=2,"not listed"=1)
worst <- function(x){ x <- x[!is.na(x) & x != ""]; if(!length(x)) return(NA_character_)
  x[which.max(unname(ifelse(is.na(sev_ord[x]), 0, sev_ord[x])))] }

gbin <- binomial(sci)
epbc_status <- nsw_status <- garnett_movement <- rep(NA_character_, length(sci))
for (i in seq_along(sci)) {
  rows <- gar[gar$binom == gbin[i], , drop = FALSE]
  if (!nrow(rows)) rows <- gar[tolower(gar[[gcommon]]) ==
                                 tolower(spl$common_name[i]), , drop = FALSE]
  if (!nrow(rows)) next
  epbc_status[i] <- worst(rows[[gepbc]]); nsw_status[i] <- worst(rows[[gnsw]])
  flags <- vapply(mvcols, function(c) any(rows[[c]] == "1", na.rm = TRUE), logical(1))
  garnett_movement[i] <- if (any(flags)) paste(mvlab[flags], collapse = ";")
    else if (any(!is.na(unlist(rows[mvcols])) & unlist(rows[mvcols]) != "")) "sedentary"
    else NA_character_
}

# ── BOTW seasonal codes (audit column) ────────────────────────────────────────
botw <- if (file.exists(MIGR)) {
  m <- read.csv(MIGR, stringsAsFactors = FALSE)
  m$botw_seasonal_codes[match(spl$common_name, m$common_name)]
} else NA_character_

out <- data.frame(
  common_name       = spl$common_name,
  scientific_name   = sci,
  migration         = unname(migration),
  garnett_movement  = garnett_movement,
  feeding_guild     = feeding_guild,
  trophic_level     = trophic_level,
  primary_lifestyle = primary_lifestyle,
  foraging_stratum  = foraging_stratum,
  habitat_breadth   = habitat_breadth,
  iucn_status       = iucn_status,
  epbc_status       = epbc_status,
  nsw_status        = nsw_status,
  botw_seasonal_codes = botw,
  stringsAsFactors  = FALSE)

write.csv(out, OUT, row.names = FALSE)
message(sprintf("Wrote %s (%d species)", OUT, nrow(out)))
message(sprintf("  AVONET matched: %d/%d | EltonTraits stratum: %d/%d",
                sum(!is.na(migration)), nrow(out),
                sum(!is.na(foraging_stratum)), nrow(out)))
cat("\nmigration:\n");     print(table(out$migration, useNA = "ifany"))
cat("\nfeeding_guild:\n"); print(table(out$feeding_guild, useNA = "ifany"))
