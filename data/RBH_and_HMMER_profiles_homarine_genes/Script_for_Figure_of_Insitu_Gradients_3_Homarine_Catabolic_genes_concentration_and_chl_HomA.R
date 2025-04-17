library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)  # for unit()

# Function to process bacterial class data
process_class_data <- function(class_data) {
  class_data %>%
    filter(domain == "Bacteria") %>%
    mutate(filter = ifelse(filter %in% c("2um", "3um"), "0.2_3", filter)) %>%
    mutate(station_numeric = as.numeric(gsub("[^0-9]", "", station))) %>%
    filter(hmm_profile %in% c("OBi1_02315_SPO3188")) %>%
    rowwise() %>%
    mutate(raw_counts_sum = sum(c_across(starts_with("transcripts_L_")), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(hmm_label = case_when(
      hmm_profile == "OBi1_02316_SPO3189" ~ "homB",
      TRUE ~ hmm_profile
    )) %>%
    mutate(class = ifelse(is.na(class) | class == "", "Unknown", class)) %>%  
    group_by(hmm_label, filter, Latitude, class) %>%
    summarise(raw_counts_sum = sum(raw_counts_sum, na.rm = TRUE), .groups = "drop")
}


# ✅ Modified function to accept station data as an argument
create_catabolic_genes_with_chla_plot <- function(class_data, chla_data, station_data, title) {
  stations <- station_data %>%
    distinct(Latitude) %>%
    arrange(Latitude)
  
  chla_summary <- chla_data %>%
    filter(depth == 6) %>%
    group_by(lat) %>%
    summarise(chla_mean = mean(chla_acs, na.rm = TRUE)) %>%
    rename(Latitude = lat)
  
  max_transcripts <- max(class_data$raw_counts_sum, na.rm = TRUE)
  max_chla <- max(chla_summary$chla_mean, na.rm = TRUE)
  chla_scale_factor <- max_transcripts / max_chla
  
  chla_summary <- chla_summary %>%
    mutate(chla_scaled = chla_mean * chla_scale_factor)
  
  ggplot() +
    # geom_rect(aes(xmin = 34, xmax = 36, ymin = -Inf, ymax = Inf),
    #           fill = "gray90", alpha = 0.1) +
    
    geom_line(data = chla_summary,
              aes(x = Latitude, y = chla_scaled),
              color = "#2f9166", size = 1, alpha = 0.25) +
    
    geom_point(data = chla_summary,
               aes(x = Latitude, y = chla_scaled),
               color = "#34B254", size = 2, alpha = 0.25) +
    
    geom_bar(data = class_data,
             aes(x = Latitude, y = raw_counts_sum, fill = class),
             stat = "identity", position = "stack", alpha = 0.8, width = 0.9,
             color = "black", size = 0.1) +
    
    # ✅ Station tick marks
    geom_rug(data = stations, aes(x = Latitude), sides = "b",
             length = unit(0.02, "npc"), color = "black", size = 0.5) +
    
    scale_x_continuous(breaks = seq(24, 40, by = 4), limits = c(24, 40)) +
    scale_y_continuous(
      name = expression("homA transcripts (x10"^6~" per Liter)"),
      labels = function(x) x / 10^6,
      limits = c(0, 4e6),
      expand = c(0.000, 0),
      sec.axis = sec_axis(~ . / chla_scale_factor,
                          name = expression("Chlorophyll-a µg/L"))
    ) +
    scale_fill_manual(values = c("#f57b78", "#be7dff", "#f1a340", "brown")) +
    labs(
      title = paste(" ", title),
      x = "Latitude",
      fill = "Bacterial Class"
    ) +
    theme_classic() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.line.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 24, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(color = "black", size = 26),
      axis.title.y.right = element_text(color = "#2f9166", size = 26),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 24),
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
      legend.position = c(0.2, 0.7),  # x and y coords (fractions of the panel: 0 to 1)
      legend.background = element_rect(fill = alpha("white", 0.7), color = "white"),
      legend.key = element_rect(fill = NA),
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 24),
      axis.ticks.x = element_blank()
    )
}

# Function for homarine concentration dual-axis plot
create_concentration_plot_dual_y <- function(homarine_data, metabs_stations_data, title) {
  km1906_data <- homarine_data %>%
    filter(Cruise == "KM1906", SampID != "MU12") %>% 
    select(Lat, replicate, Part.Conc.nM, Diss.Conc.nM)
  
  max_diss <- max(km1906_data$Diss.Conc.nM, na.rm = TRUE)
  max_part <- max(km1906_data$Part.Conc.nM, na.rm = TRUE)
  scale_factor <- max_diss / max_part
  
  km1906_data <- km1906_data %>%
    mutate(Part.Rescaled = Part.Conc.nM * scale_factor) %>%
    pivot_longer(cols = c(Diss.Conc.nM, Part.Rescaled), 
                 names_to = "Fraction", 
                 values_to = "Concentration") %>%
    mutate(
      Fraction = recode(Fraction,
                        "Diss.Conc.nM" = "Dissolved",
                        "Part.Rescaled" = "Particulate")
    )
  
  ggplot(km1906_data, aes(x = Lat, y = Concentration, color = Fraction, shape = Fraction)) +
    geom_point(size = 5, alpha = 0.6) +
    geom_rug(data = metabs_stations_data %>% distinct(Lat), 
             aes(x = Lat), 
             inherit.aes = FALSE,
             sides = "b", 
             length = unit(0.02, "npc"), 
             color = "black", 
             size = 0.5) +
    scale_color_manual(values = c("Dissolved" = "#f5948e", "Particulate" = "#67a9cf")) +
    scale_shape_manual(values = c("Dissolved" = 16, "Particulate" = 17)) +
    scale_y_continuous(
      name = "Dissolved (nM)",
      sec.axis = sec_axis(~ . / scale_factor, name = "Particulate (nM)")
    ) +
    scale_x_continuous(breaks = seq(24, 40, by = 4), limits = c(24, 40)) +
    labs(
      title = paste(" ", title),
      x = "Latitude",
      color = "Homarine",
      shape = "Homarine"
    ) +
    theme_classic() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = 24, face = "bold"),
      axis.title.x = element_text(size = 26),
      axis.title.y = element_text(color = "#f5948e", size = 26),
      axis.title.y.right = element_text(color = "#67a9cf", size = 26),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
      axis.text.y = element_text(size = 24),
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
      legend.position = c(0.2, 0.7),
      legend.background = element_rect(fill = alpha("white", 0.7), color = "white"),
      legend.key = element_rect(fill = NA),
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 24),
    )
}


# === Load Input Files ===
class_data_files <- c("C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/hmmr_profiles/HMMR_G5/Gradients_3/SPO3188_SPO3189_combined_G3_with_metadata.csv")
homarine_data_path <- "C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/Data/Gradients_cruises/Josh/Homarine_Data_For_Frank.csv"
chlorophyll_path <- "C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/Data/Gradients_cruises/Chlorophyll_a_from_AC_S_Instrument.csv"
stations_path <- "C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/Data/Gradients_cruises/Gradients_3_stations.csv"  # or swap with a different file
metabs_stations_path <- "C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/Data/Gradients_cruises/Gradients_3_stations_metabs.csv"  # or swap with a different file

# === Read Data ===
class_data <- read.csv(class_data_files[1])
homarine_data <- read.csv(homarine_data_path)
chlorophyll_data <- read.csv(chlorophyll_path) %>% filter(depth == 6)
station_data <- read.csv(stations_path)
metabs_stations_data <- read.csv (metabs_stations_path)

# === Process and Plot ===
collapsed_class_data <- process_class_data(class_data)

combined_gene_chla_plot <- create_catabolic_genes_with_chla_plot(
  collapsed_class_data,
  chlorophyll_data,
  station_data,
  " "
)

concentration_plot <- create_concentration_plot_dual_y(homarine_data, metabs_stations_data, "")
spacer <- plot_spacer()

final_plot <- (
  combined_gene_chla_plot +
    spacer +
    concentration_plot
) + plot_layout(ncol = 1, heights = c(1, 0.1, 1))

print(final_plot)

# Optional Save
# ggsave("homarine_combined_plot.png", plot = final_plot, width = 12, height = 16)
