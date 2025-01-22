library(tidyverse)
library(rstatix)
library(rnaturalearthdata)
library(rnaturalearth)

#TODO KRH: Use this as starting point for map for homarine experiment

# Load data
g5.isolate.file <- "figures/maps/Isolation_Locations_01142025.csv"

g5.iso.dat <- read_csv(g5.isolate.file, show_col_types = FALSE) %>%
  select(isolate, group, Lat, Lon, Depth_m)

all.iso.dat <- g5.iso.dat

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Function to plot semi-circles for overlapping points
create_semi_circles <- function(lon, lat, radius = 1, angle_start = 0, angle_end = 180, fill) {
  tibble(
    x = lon + radius * cos(seq(angle_start, angle_end, length.out = 100) * pi / 180),
    y = lat + radius * sin(seq(angle_start, angle_end, length.out = 100) * pi / 180),
    fill = fill
  )
}

# Identify overlapping points
overlapping_points <- all.iso.dat %>%
  filter(group %in% c("RW", "PS")) %>%
  group_by(Lon, Lat) %>%
  filter(n() > 1) %>%
  ungroup()

# Generate semi-circle data for overlaps
semi_circle_data <- overlapping_points %>%
  rowwise() %>%
  mutate(
    data = list(bind_rows(
      create_semi_circles(Lon, Lat, radius = 1, angle_start = 0, angle_end = 180, fill = "RW"),
      create_semi_circles(Lon, Lat, radius = 1, angle_start = 180, angle_end = 360, fill = "PS")
    ))
  ) %>%
  select(data) %>%
  unnest(data)

# Plot the map with semi-circles
hom_fate_map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-155, -55), ylim = c(-5, 50)) + # Restrict map to the west of 100°W
  geom_point(data = all.iso.dat %>% filter(!group %in% c("RW", "PS")), aes(x = Lon, y = Lat, fill = group), size = 5, shape = 21) +
  geom_polygon(data = semi_circle_data, aes(x = x, y = y, fill = fill), alpha = 0.8, color = "black", size = 0.3) +
  # geom_text(data = all.iso.dat, aes(x = Lon, y = Lat, label = Depth_m),
  #           nudge_x = 5, nudge_y = 0.5, size = 3, color = "black") +
  scale_fill_manual(values = c("#86C6E5", "lightpink", "#9467BD", "#FF9896", "#8DCA8D", "#7FDBFF", "#FFDC00")) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, color = "black", stroke = 0.5))) + # Ensure legend symbols are circular
  theme_bw() +
  labs(fill = "Cruise", # Add legend label
       title = "Isolates grown in Homarine",
       x = "Longitude",
       y = "Latitude")

# remove everything but p
rm(all.iso.dat, g5.iso.dat, g5.isolate.file, overlapping_points, semi_circle_data, world)
