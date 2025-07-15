# === Load Libraries ===
library(ggplot2)
library(dplyr)
library(readr)
library(grid)
library(patchwork)

# === Load Data ===
class_data <- read_csv("C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/hmmr_profiles/HMMR_G5/Gradients_3/SPO3188_SPO3189_combined_G3_revised_annotations_with_metadata_final.csv")
homarine_data <- read_csv("C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/Data/Gradients_cruises/Josh/Homarine_Data_For_Frank.csv")
chlorophyll_data <- read_csv("C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/Data/Gradients_cruises/Chlorophyll_a_from_AC_S_Instrument.csv")

# === Process Data ===
collapsed_class_data <- class_data %>%
  filter(domain == "Bacteria") %>%
  filter(hmm_profile == "OBi1_02316_SPO3189") %>%
  rowwise() %>%
  mutate(raw_counts_sum = sum(c_across(starts_with("transcripts_L_")), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(class = ifelse(is.na(class) | class == "", "Unknown", class)) %>%
  group_by(Latitude, class) %>%
  summarise(raw_counts_sum = sum(raw_counts_sum), .groups = "drop")

homarine_df <- homarine_data %>%
  filter(Cruise == "KM1906", SampID != "MU12") %>%
  mutate(Total = Diss.Conc.nM + Part.Conc.nM) %>%
  rename(Latitude = Lat)

chla_df <- chlorophyll_data %>%
  filter(depth == 6) %>%
  group_by(lat) %>%
  summarise(chla_mean = mean(chla_acs, na.rm = TRUE)) %>%
  rename(Latitude = lat)

# === Define Scaling ===
chla_scale <- 6  # Stretch chlorophyll visually
max_homarine <- max(homarine_df$Total, na.rm = TRUE)
max_transcripts <- max(collapsed_class_data$raw_counts_sum, na.rm = TRUE)
metaT_scale <- (max_homarine / max_transcripts) * 0.6

# === Apply Scaling ===
collapsed_class_data <- collapsed_class_data %>%
  mutate(metaT_scaled = raw_counts_sum * metaT_scale)

chla_df <- chla_df %>%
  mutate(chla_scaled = chla_mean * chla_scale)

# === Main Plot: Homarine, metaT, and scaled Chl-a ===
p_main <- ggplot() +
  # Chlorophyll (scaled, bottom layer)
  geom_line(data = chla_df,
            aes(x = Latitude, y = chla_scaled),
            color = "#2f9166", size = 1.1, alpha = 0.8) +
  geom_point(data = chla_df,
             aes(x = Latitude, y = chla_scaled),
             color = "#34B254", size = 2, alpha = 0.8) +
  
  # metaT
  geom_bar(data = collapsed_class_data,
           aes(x = Latitude, y = metaT_scaled, fill = class),
           stat = "identity", position = "stack", color = "black", size = 0.1,
           alpha = 0.85, width = 0.9) +
  
  # Homarine
  geom_point(data = homarine_df,
             aes(x = Latitude, y = Total, fill = "Homarine (nM)"),
             shape = 21, size = 5, stroke = 0.5, color = "black", alpha = 0.85) +
  
  # Axes
  scale_y_continuous(
    name = "Homarine (nM)",  # Left y-axis
    limits = c(0, max_homarine * 1.1),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . / metaT_scale,
                        name = expression("homB transcripts (x10"^7~" per Liter)"),
                        labels = function(x) x / 1e7)
  ) +
  scale_fill_manual(
    name = "Bacterial Class",
    values = c(
      "Alphaproteobacteria" = "#f57b78",
      "Gammaproteobacteria" = "#be7dff",
      "Unknown" = "#f1a340",
      "Homarine (nM)" = "#67a9cf"
    )
  ) +
  scale_color_manual(
    name = "",  # Legend title (blank or specify e.g. "Measurement")
    values = c("Homarine (nM)" = "#67a9cf")
  ) +
  scale_x_continuous(name = "Latitude", limits = c(24, 42), breaks = seq(24, 40, 4)) +
 theme_classic() +
theme(
  axis.title.y = element_text(color = "#67a9cf", size = 22),
  axis.title.y.right = element_text(color = "black", size = 22),
  axis.title.x = element_text(size = 22),
  axis.ticks.x = element_line(size = 1.2, color = "black"),  # <- RE-ADDED X TICKS
  axis.text = element_text(size = 20),
  legend.position = c(0.15, 0.8),
  legend.title = element_text(size = 20),
  legend.text = element_text(size = 18),
  legend.background = element_rect(fill = alpha("white", 0.7), color = "white"),
  plot.margin = margin(t = 20, r = 80, b = 20, l = 20)
)



# === Fake Plot to Add Chlorophyll Y-Axis ===
p_chla_axis <- ggplot(chla_df, aes(x = Latitude, y = chla_mean)) +
  geom_blank() +
  scale_y_continuous(
    name = "Chlorophyll-a (µg/L)",
    limits = c(0, max(chla_df$chla_mean, na.rm = TRUE) * 1.05),
    expand = c(0, 0),
    position = "right"
  ) +
  theme_void() +
  theme(
    axis.title.y.right = element_text(color = "#2f9166", size = 22, angle = 90, vjust = 1.5),
    axis.text.y.right = element_text(color = "#2f9166", size = 18),
    axis.ticks.y.right = element_line(color = "#2f9166"),
    axis.line.y.right = element_line(color = "#2f9166")
  )

# === Combine ===
final_plot <- p_main + p_chla_axis + plot_layout(ncol = 1)
print(final_plot)

# === Optional Save ===
# ggsave("final_3axis_plot_chla_stretched.png", final_plot, width = 14, height = 8)
