suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(terra)
  library(geodata)
  library(ggplot2)
})

CACHE <- "ebirdabund_cache"

message("Building NSW polygon (+ 100 km buffer)...")
aus     <- geodata::gadm(country = "AUS", level = 1, path = CACHE)
nsw     <- sf::st_as_sf(aus[aus$NAME_1 == "New South Wales", ])
polygon <- sf::st_transform(
  sf::st_buffer(sf::st_transform(nsw, 3577), 100000),
  4326
)

message("Building covariate stack (downloads Meta/WRI CHMv2 first run)...")
cov <- prepare_covariates(polygon, cache_dir = CACHE)

message("Generating verification map...")
th <- cov[["tree_height"]]

df <- as.data.frame(th, xy = TRUE)
names(df) <- c("x", "y", "tree_height")
df <- df[!is.na(df$tree_height), ]

nsw_border <- sf::st_transform(nsw, 4326)

p <- ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = tree_height)) +
  geom_sf(data = nsw_border, fill = NA, colour = "black", linewidth = 0.4) +
  scale_fill_gradientn(
    colours = c("white", "#cde9c0", "#5ab566", "#1a7337", "#003d12"),
    na.value = NA,
    name = "Tree height (m)"
  ) +
  coord_sf(expand = FALSE) +
  labs(title = "Tree height — Meta/WRI CHMv2 (DINOv3), OL=4 ≈ 38 m") +
  theme_minimal(base_size = 11)

ggsave("tree_height_chmv2.png", p, width = 8, height = 7, dpi = 150)
message("Saved tree_height_chmv2.png")
