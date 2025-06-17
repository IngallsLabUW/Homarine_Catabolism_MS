library(RaMS)
# library(Rdisop)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
gc()
rm(list = ls())
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))

# Function to plot D2 or D3 chromats
plot_expected_iso <- function(
        mz_iso,
        rt_iso,
        mz_mono,
        sample_key
        ) {
    ppm = 5
  my_y_scale <- scale_y_continuous(
    expand = c(0, 0), limits = c(0, NA), labels = function(x) format(x, scientific = TRUE)
  )
  my_theme <- theme(strip.text = element_blank())

  ## plot unlabeled compound ----
  mz_min <- mz_mono - mz_mono * ppm / 1e6
  mz_max <- mz_mono + mz_mono * ppm / 1e6
  g_ms1_og <- plot_EIC2(
      ms1data = ms1data$MS1 %>%
          filter(
              filename %in% c(
                  basename(g4_c_samples),
                  basename(g4_h_samples)
              )
          ),
      m_z = mz_mono,
      mz_ppm = ppm,
      r_t = rt_iso,
      rt_buffer = 2.5,
      samples_oi = sample_key
  )   +
      annotate("text",
               x = -Inf,
               y = Inf,
               hjust = 0,
               vjust = 2,
               label = paste0(
                   "Monoisotopic (m/z = ", round(mz_min, 5),
                   " - ", round(mz_max, 5), ")"
               )
      ) +
      colScale +
      my_y_scale +
      my_theme
  
  # Plot labeled version
  mz_min <- mz_iso - mz_iso * ppm / 1e6
  mz_max <- mz_iso + mz_iso * ppm / 1e6
  g_ms1_og_label <- plot_EIC2(
      ms1data = ms1data$MS1 %>%
          filter(
              filename %in% c(
                  basename(g4_c_samples),
                  basename(g4_h_samples)
              )
          ),
      m_z = mz_iso,
      mz_ppm = ppm,
      r_t = rt_iso,
      rt_buffer = 2.5,
      samples_oi = sample_key
  )   +
      annotate("text",
               x = -Inf,
               y = Inf,
               hjust = 0,
               vjust = 2,
               label = paste0(
                   "Labeled (m/z = ", round(mz_min, 5),
                   " - ", round(mz_max, 5), ")"
               )
      ) +
      colScale +
      my_y_scale +
      my_theme
  
  # Combine them
  g_out = g_ms1_og / g_ms1_og_label
  
  return(g_out)
  
}

# GET FILES ------
## Set dirs ------
output_dir <- here("figures", "exploratory", "metabolomics", "combined_chromat_figs", "g4_chromats")
g4_meta_data_dir <- here("data", "raw", "metabolomics", "G4")
g4_raw_dat_dir <- here("data", "raw", "metabolomics", "G4", "particulate", "mzml_skyline")
input_dir <- here("data", "intermediate", "metabolomics", "g4")

## Read in data ------
g4_detected_features <- read_csv(
  here(
    input_dir, "grouped_iso_results.csv"
  ),
  show_col_types = FALSE
)

## Read in sample keys and combine -------
g4_sample_key <- read_csv(
  here(
    g4_meta_data_dir, "sample_key.csv"
  ),
  show_col_types = FALSE
) %>%
  rename(filename = replicate_name) %>%
  select(filename, treatment) %>%
  mutate(experiment = "g4") %>%
  mutate(filename = paste0(filename, ".mzML"))

## Set individual example samples ----
g4_h_samples <- c(
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t2_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t2_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t2_C.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_H_t2_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_H_t2_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_H_t2_C.mzML")
)
g4_c_samples <- c(
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C-H_t2_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C-H_t2_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C-H_t2_C.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_C-H_t2_A.mzML")
#  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_C-H_t2_B.mzML")
#  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_C-H_t2_C.mzML")
)

smps_to_plot <- c(g4_h_samples, g4_c_samples)

sample_key <- g4_sample_key %>%
  mutate(treatment_short = case_when(
    filename %in% basename(g4_h_samples) ~ "Plus Homarine",
    filename %in% basename(g4_c_samples) ~ "Control",
    TRUE ~ "NA"
  )) %>%
  # cast treatment_short as factor
  mutate(treatment_short = factor(treatment_short, levels = c("Control", "Plus Homarine"))) %>%
  filter(treatment_short != "NA") %>%
  distinct()

# PLOT POSITIVE MASS FEATURES Read in MS1 data ------
ms1data <- grabMSdata(files = smps_to_plot, grab_what = c("MS1"))

### Prep color scale to share between plots
myColors <- brewer.pal(3, "Set1")
names(myColors) <- levels(sample_key$treatment_short)
colScale <- scale_colour_manual(name = "Treatment", values = myColors, drop = FALSE)

# get unique compounds from g4_detected_features, grabbing name, mz, and rt
compounds_to_plot <- g4_detected_features %>%
  select(id, isotopologue_type, mz_med, scan_time_med, mode, monoisotopic_mass) %>%
  distinct()

# PLOT EACH POTENTIAL HIT -----        

# Note these are all positive so no need to get data for negative mode data
for (i in 1:nrow(compounds_to_plot)) {
  compound <- compounds_to_plot[i, ]
  # check if compound is positive, if so plot
  if (compound$mode == "positive") {
    g <- plot_expected_iso(
        mz_iso = compounds_to_plot$mz_med[i],
        rt_iso = compounds_to_plot$scan_time_med[i],
        mz_mono = compounds_to_plot$monoisotopic_mass[i],
        sample_key = sample_key
    )
    filename <- paste0(
        "for_inspection_D3_D2_ID",
        compounds_to_plot$id[i],
        ".pdf"
    )
    ggsave(
      here(output_dir, filename),
      g,
      width = 12,
      height = 10,
      units = "in",
      dpi = 300
    )
  }
}
