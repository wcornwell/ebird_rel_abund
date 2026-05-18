suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(ggrepel)
})

CACHE      <- "ebirdabund_cache"
STACK_PATH <- "species_maps/nsw_abundance_stack_3km.tif"
TAXONOMY   <- "nsw_ebird_taxonomy.csv"

# NSW boundary
aus <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])

# Species display names: safe-name → common name
taxonomy    <- read.csv(TAXONOMY, stringsAsFactors = FALSE)
sp_names_vec <- setNames(
  taxonomy$common_name,
  ebirdabund:::safe_name(taxonomy$common_name)
)

message("Building pairwise comparison plot for NSW...")
p <- plot_abundance_comparison(
  stack         = STACK_PATH,
  polygon       = nsw,
  cache_dir     = CACHE,
  species_names = sp_names_vec,
  top_label_n   = 15,
  hex_spacing_km = 5
)

out_path <- "abundance_comparison_nsw.png"
ggplot2::ggsave(out_path, p, width = 22, height = 14, dpi = 150)
message("Saved: ", out_path)
