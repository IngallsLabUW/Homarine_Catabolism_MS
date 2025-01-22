library(tidyverse)
library(rstatix)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggrepel)

# TODO KRH: Use this as starting point for map for homarine experiment

# Load data
hom_fate_file <- "figures/maps/2025_1_14_hom_fate_coords_map.csv"

hom_fate_dat <- read_csv(hom_fate_file, show_col_types = FALSE) %>%
  select(cruise, lat, lon, station) %>%
  arrange(cruise, station) %>%
  mutate(label = paste(cruise, station, sep = " St"))
  # Add year here?

  # Load world map
  world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot the data with a reduced map extent
hom_fate_map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(
    xlim = c(-160, -100), ylim = c(-10, 55)
  ) + # Restrict map to the west of 100°W
  geom_point(
    data = hom_fate_dat, aes(x = lon, y = lat, shape = cruise), size = 3
  ) + # Larger circles
  geom_text_repel(
    data = hom_fate_dat, aes(x = lon, y = lat, label = label),
    size = 3
 ) +
  theme_bw() +
  labs(
      shape = "Cruise", # Add legend label
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = "none",
    legend.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  )
hom_fate_map
# remove everything but p
rm("hom_fate_file", "hom_fate_dat", "world")
