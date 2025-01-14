library(tidyverse)
library(rstatix)
library(sf)
library(rnaturalearthdata)
library(rnaturalearth)

# Load data
g5.isolate.file <- "C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/manuscript/map/2025_1_14_hom_fate_coords_map.csv"

g5.iso.dat <- read_csv(g5.isolate.file) %>%
  select(cruise, lat, lon, station)

all.iso.dat <- g5.iso.dat

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot the data with a reduced map extent
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-160, -100), ylim = c(-10, 55)) + # Restrict map to the west of 100°W
  geom_point(data = all.iso.dat, aes(x = lon, y = lat, fill = cruise), size = 5, shape = 21) + # Larger circles
  geom_text(data = all.iso.dat, aes(x = lon, y = lat, label = station), 
            nudge_x = 3, nudge_y = 0.5, size = 3, color = "black") +
  scale_fill_manual(values = c("#86C6E5", "#9467BD", "#FF9896")) + # Use custom colors
  theme_bw() +
  labs(fill = "Cruise", # Add legend label
       title = "Homarine Fate Experiments",
       x = "Longitude",
       y = "Latitude")
