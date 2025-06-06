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


##NCBI Sunburst plot -----
source(here("figures", "isolates_genesearches", "ncbi_genome_taxonomy_summary_sunburst.R"))

## NCBI Map Figure -----
source(here("figures", "maps", "ncbi_location_map.R"))
ncbi_map <- ncbi_map_2 
# remove everything except ncbi_map
rm(list = ls()[!ls() %in% c("ncbi_map", "isolate_plot", "sunburst")])


# Combine figures -----
isolate_plot_legend <- wrap_elements(panel = isolate_plot + theme(legend.position = "bottom",
                                                                  plot.margin = margin(0, 0, 0, 0) ))
ncbi_map_layout <- wrap_elements(panel = ncbi_map + theme(legend.position = "bottom"))
sunburst_layout <- wrap_elements(panel = sunburst + theme(legend.position = "right"))
ncbi_sunburst_layout <- sunburst_layout / ncbi_map_layout +
   plot_layout(ncol = 1, heights = c(1, 3)) 
combined_fig <- isolate_plot_legend + ncbi_sunburst_layout +
  plot_layout(ncol = 2, widths = c(1, 1.3)) +
    plot_annotation(tag_levels = "A")
combined_fig
ggsave(here("figures", "figure_3_combined.pdf"),
       plot = combined_fig,
       width =7, height = 4, dpi = 300,
       device = cairo_pdf)
