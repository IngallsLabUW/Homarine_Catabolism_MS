library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Load transect data
transect_file <- "C:/SynologyDrive/GITHUB/Homarine/Homarine_Catabolism_MS/figures/maps/Chlorophyll_a_from_AC_S_Instrument.csv" # Replace with your file path
transect_data <- read_csv(transect_file)

# Select unique lat/lon points
transect_unique <- transect_data %>%
  select(lat, lon) %>%
  distinct()

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot the map with transect line
p <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-160, -156), ylim = c(25.5, 44)) + # Adjust limits as needed
  geom_path(data = transect_unique, aes(x = lon, y = lat), color = "#008080", size = 1.2) + # Plot transect line with custom color
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  labs(
    title = "Gradients 3",
    x = "Longitude",
    y = "Latitude"
  )

print(p)
