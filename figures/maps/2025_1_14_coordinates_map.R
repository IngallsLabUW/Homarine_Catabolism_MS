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

# Plot the data with different colors for each cruise and adjusted station labels
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-160, -50), ylim = c(-10, 55)) +
  geom_point(data = all.iso.dat, aes(x = lon, y = lat, fill = cruise), size = 3, shape = 21) +
  geom_text(data = all.iso.dat, aes(x = lon, y = lat, label = station), 
            nudge_x = 2, nudge_y = 0.5, size = 3, color = "black") +
  scale_fill_viridis_d(option = "C") + # Use a colorblind-friendly palette
  theme_bw() +
  labs(fill = "Cruise", # Add legend label
       title = "Isolate Distribution by Cruise",
       x = "Longitude",
       y = "Latitude")
