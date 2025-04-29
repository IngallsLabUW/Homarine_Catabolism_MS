# This script makes a map of the lat long coordinates of the NCBI metadata with the most common environmental source for each lat long point

# Load required libraries ----
library(here)
library(tidyverse)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggplot2)
library(ggstar)
fnt_size <- 5
# Define a function to return the mode value from a vector of values ------
get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

# Load the metadata file -----
dat <- read_csv(here("data", "intermediate", "ncbi_metadata_clean.csv"),
                     show_col_types = FALSE) %>%
    rename(assembly = assembly_accession)
dat_isolates <- read_csv(here("figures", "maps", "Isolation_Locations_01142025.csv"),
                     show_col_types = FALSE) %>%
    filter(group != "BH")%>%
    select(Lat, Lon) %>%
    distinct() 

# Summarize data by isolation_source_mapped and environmental_source_mapped and make a pie plot
dat_pie <- dat %>%
    mutate(pie_section = 
               case_when(
            isolation_source_mapped == "non-environmental" ~ "non-environmental",
            isolation_source_mapped == "unknown" ~ "unknown",
            environmental_source_mapped == "marine" ~ "environmental (marine)",
            environmental_source_mapped == "freshwater" ~ "environmental (freshwater)",
            environmental_source_mapped == "sediment" ~ "environmental (sediment)",
            environmental_source_mapped == "soil" ~ "environmental (soil)",
            environmental_source_mapped == "other" ~ "environmental (other)"
            )
    ) %>%
    mutate(pie_section = 
               factor(pie_section, 
                      levels = c("non-environmental", "unknown", 
                                 "environmental (marine)", "environmental (freshwater)", "environmental (sediment)",
                                 "environmental (soil)", "environmental (other)")) )%>%
    group_by(pie_section) %>%
    summarise(n = n()) %>%
    ungroup() 

# Make pie plot
dat_pie %>%
    ggplot(aes(x = "", y = n, fill = pie_section)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("#9467BD", "grey80", "#1F77B4", "#17BECF", "#613613", "#F28E2B", "grey50")) +
    labs(fill = "Isolation Source") 


# Summarize by unique lat and lon, grabbing the most common environmental source for each
dat_plot2 <- dat %>%
    filter(isolation_source_mapped == "environmental") %>% # Only keep environmental sources
    # round lat and lon to 3 decimal places for plotting
    mutate(lat = round(lat, 3), lon = round(lon, 3)) %>%
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
        ))) %>%
    filter(!is.na(lat))

# Plot the lat long on a map!
world <- ne_countries(scale = "medium", returnclass = "sf")

####Scale point size by n
ncbi_map_2 <- ggplot(data = world) +
    geom_sf(fill = "grey85", color = "grey85") +
    coord_sf() + 
    geom_point(
        data = dat_plot2, 
        aes(x = lon, y = lat, fill = environmental_source_mapped, size = n), 
        shape = 21, 
        stroke = 0.15
        # No fixed size here since size is mapped to 'n'
    ) +
    geom_star(
        data = dat_isolates, 
        aes(x = Lon, y = Lat), 
        color = "black",
        starstroke = 0.1,
        fill = "#9467BD", 
        size = 1 # Adjusting star size to be smaller
    ) +
    labs(
        fill = "Environmental Type", 
        size = "# of Genomes",
        x = "Longitude",
        y = "Latitude"
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(size = 3),  # Make the legend key smaller
            title.position = "top", 
            title.hjust = 0.5, 
            nrow = 2
        ),
        size = guide_legend(
            title.position = "top", 
            title.hjust = 0.5
        )
    ) +    
    scale_fill_manual(values = c("#1F77B4", "#17BECF", "#613613", "#F28E2B", "grey40")) +
    scale_size_continuous(transform = "log10", range = c(0.5, 3)) +  # Adjust range to control size scaling
    theme_bw() +
    theme(
        panel.border = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_text(size = fnt_size),
        legend.title = element_text(size = fnt_size),
        legend.text = element_text(size = fnt_size),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key.size = unit(0.1, 'cm'),  # Smaller key size
        legend.margin = margin(-15, 0, 0, 0)
    )
ncbi_map_2

#ggsave(ncbi_map_2, filename = "figures/maps/ncbi_core4_map.png", dpi = 600, height = 7, width = 9, scale = 1.1)







































