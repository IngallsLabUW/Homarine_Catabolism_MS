# Script to make figure 3, which combines th isolate growth heatmap, with the NCBI phylogenetic info, and the map

# Load libraries -----
library(tidyverse)
library(here)
library(patchwork)

# Source figures from other scripts -----
## Isolate growth heatmap figure ----
source(here("figures", "isolates_genesearches", "Isolate_Growth_Heatmap_Figure.R"))
isolate_plot <- combined.plot
# remove everything except isolate_plot
rm(list = ls()[!ls() %in% "isolate_plot"])



## NCBI Map Figure -----
source(here("figures", "maps", "ncbi_location_map.R"))
ncbi_map <- ncbi_map_2 
# remove everything except ncbi_map
rm(list = ls()[!ls() %in% c("ncbi_map", "isolate_plot")])


# Combine figures -----
isolate_plot_legend <- wrap_elements(panel = isolate_plot + theme(legend.position = "right"))
ncbi_map_layout <- wrap_elements(panel = ncbi_map + theme(legend.position = "bottom"))
combined_fig <- isolate_plot_legend + ncbi_map_layout +
    plot_annotation(tag_levels = "A")
combined_fig
