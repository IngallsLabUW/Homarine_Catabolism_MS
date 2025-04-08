library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork) # For multi-panel figures

# Load transect data for all three files
transect_file1 <- "C:/SynologyDrive/GITHUB/Homarine/Homarine_Catabolism_MS/figures/maps/CTD_Chlorophyll_a_Gradients_1.csv"
transect_file2 <- "C:/SynologyDrive/GITHUB/Homarine/Homarine_Catabolism_MS/figures/maps/chlorophyll_a_from_AC_S_instrument_Gradients_2.csv"
transect_file3 <- "C:/SynologyDrive/GITHUB/Homarine/Homarine_Catabolism_MS/figures/maps/Chlorophyll_a_from_AC_S_Instrument.csv"

transect_data1 <- read_csv(transect_file1) %>%
  select(lat, lon) %>%
  distinct() %>%
  mutate(transect = "Transect 1")

transect_data2 <- read_csv(transect_file2) %>%
  select(lat, lon) %>%
  distinct() %>%
  mutate(transect = "Transect 2")

transect_data3 <- read_csv(transect_file3) %>%
  select(lat, lon) %>%
  distinct() %>%
  mutate(transect = "Transect 3")

# Combine all transects into one data frame
all_transects <- bind_rows(transect_data1, transect_data2, transect_data3)

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create individual plots for each transect
p1 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-160, -156), ylim = c(20, 44)) +
  geom_path(data = filter(all_transects, transect == "Transect 1"), aes(x = lon, y = lat), color = "#008080", size = 1.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  scale_x_continuous(breaks = seq(-160, -156, by = 2)) + # Customize x-axis breaks
  labs(title = "G1", x = " ", y = "Latitude")

p2 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-160, -156), ylim = c(20, 44)) +
  geom_path(data = filter(all_transects, transect == "Transect 2"), aes(x = lon, y = lat), color = "#008080", size = 1.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    axis.title.y = element_blank(), # Remove y-axis title
    axis.text.y = element_blank(),  # Remove y-axis text (degrees)
    axis.ticks.y = element_blank()  # Remove y-axis ticks
  ) +
  scale_x_continuous(breaks = seq(-160, -156, by = 2)) + # Customize x-axis breaks
  labs(title = "G2", x = "Longitude")

p3 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-160, -156), ylim = c(20, 44)) +
  geom_path(data = filter(all_transects, transect == "Transect 3"), aes(x = lon, y = lat), color = "#008080", size = 1.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    axis.title.y = element_blank(), # Remove y-axis title
    axis.text.y = element_blank(),  # Remove y-axis text (degrees)
    axis.ticks.y = element_blank()  # Remove y-axis ticks
  ) +
  scale_x_continuous(breaks = seq(-160, -156, by = 2)) + # Customize x-axis breaks
  labs(title = "G3", x = " ")

# Combine the plots into a 3-panel figure arranged horizontally
combined_plot <- p1 | p2 | p3  # Arrange plots in a horizontal layout

# Print the combined plot
print(combined_plot)

