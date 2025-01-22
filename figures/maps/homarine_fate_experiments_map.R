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
  geom_sf(fill = "grey80", color = "grey50") +
  coord_sf(
    xlim = c(-160, -100), ylim = c(-5, 55)
  ) + # Restrict map to the west of 100°W
  geom_point(
    data = hom_fate_dat, aes(x = lon, y = lat, shape = cruise), size = 2
  ) + # Larger circles
  geom_text_repel(
    data = hom_fate_dat, aes(x = lon, y = lat, label = label),
    size = 2.5,
    # make bold
    fontface = "bold",
 ) +
    # Make x axis have fewer breaks
  scale_x_continuous(breaks = c(-160, -130, -100)) +
    scale_y_continuous(breaks = c(-5, 0, 25, 50)) +
  labs(
      shape = "Cruise", # Add legend label
    x = "Longitude",
    y = "Latitude"
  ) +
    theme_bw() +

  theme(
    legend.position = "none",
    legend.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  )
hom_fate_map

# remove everything but p
rm("hom_fate_file", "hom_fate_dat", "world")
