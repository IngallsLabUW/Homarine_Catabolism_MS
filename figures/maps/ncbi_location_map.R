# This script makes a map of the lat long coordinates of the NCBI metadata with the most common environmental source for each lat long point


# Load required libraries ----
library(here)
library(tidyverse)
library(rnaturalearthdata)
library(rnaturalearth)

# Define a function to return the mode value from a vector of values ------
get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

# Load the metadata file -----
dat <- read_csv(here("data", "intermediate", "ncbi_metadata_clean.csv"),
                     show_col_types = FALSE)

## FILTER DAT dataframe to only those that have homi genes here!!!!!!-----


# Summarize by unique lat and lon, grabbing the most common environmental source for each
dat_plot2 <- dat %>%
    filter(isolation_source_mapped == "environmental") %>% # Only keep environmental sources
    group_by(lat, lon) %>%
    summarise(
        n = n(),
        # TODO Change this to mode of environmental source or prepare for pie graphs
        environmental_source_mapped = get_mode(environmental_source_mapped)
    ) %>%
    ungroup() %>%
    # Make environmental_source_mapped factor
    mutate(environmental_source_mapped = factor(
        environmental_source_mapped,
        levels = 
            c(
            "marine", "freshwater", "sediment", "soil", "other"
        )))

# Plot the lat long on a map!
world <- ne_countries(scale = "medium", returnclass = "sf")

ncbi_map <- ggplot(data = world) +
    geom_sf(fill = "grey80", color = "grey80") +
    coord_sf() + 
    geom_point(
        data = dat_plot2, aes(x = lon, y = lat, color = environmental_source_mapped), size = 0.5
    ) +
    labs(
        color = "Environmental Type", # Add legend label
        x = "Longitude",
        y = "Latitude"
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))+
    scale_color_manual(values = c("blue", "green", "orange", "red", "purple")) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_text(size = 7)
    )
ncbi_map
