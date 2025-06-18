# This script is used to visualize the chromatograms of detected features in the G4 dataset.
library(RaMS)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
gc()
rm(list = ls())
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))

# GET FILES ------
## Set dirs ------
output_dir <- here("figures", "metabolomics", "exploratory", "g4_chromats_inspect")
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

## LOOP FOR EACH FEATURE
for (i in 1:nrow(g4_detected_features)){
feature_df <- g4_detected_features %>% 
    filter(id == i)
samples <- feature_df$samples %>%
    str_split(",") %>%
    unlist() %>%
    str_trim() %>%
    unique() %>%
    paste0(".mzML")
data_dir <- here(g4_raw_dat_dir, feature_df$mode)
file_paths <- here(data_dir, samples)

ms1data <- grabMSdata(files = file_paths, grab_what = c("MS1"))$MS1

samples_oi <- g4_sample_key %>%
    filter(filename %in% samples)

# Set parameters for plotting ------
r_t <- feature_df$scan_time_med
rt_buffer <- 2.5 # minutes
m_z_iso <- feature_df$mz_med
m_z_mono <- feature_df$monoisotopic_mass
mz_ppm <- 5


rt_dat <- ms1data %>%
    select(rt, filename) %>%
    distinct() %>%
    filter(rt %between% c(r_t - rt_buffer, r_t + rt_buffer))

ms1_data_iso <- ms1data[mz %between% pmppm(m_z_iso, mz_ppm) & rt %between% c(r_t - rt_buffer, r_t + rt_buffer)] %>%
    full_join(rt_dat, by = join_by(rt, filename)) %>%
    left_join(
        samples_oi %>%
            select(filename, experiment),
        by = join_by(filename)
    ) %>%
    mutate(int = ifelse(is.na(int), 0, int)) %>%
    mutate(label = feature_df$isotopologue_type)

ms1_data_mono <- ms1data[mz %between% pmppm(m_z_mono, mz_ppm) & rt %between% c(r_t - rt_buffer, r_t + rt_buffer)] %>%
    full_join(rt_dat, by = join_by(rt, filename)) %>%
    left_join(
        samples_oi %>%
            select(filename, experiment),
        by = join_by(filename)
    ) %>%
    mutate(int = ifelse(is.na(int), 0, int)) %>%
    mutate(label = "monoisotopic")

# Plotting ------
# start pdf
pdf(
    file = here(output_dir, paste0("feature_", feature_df$id, ".pdf")),
    width = 8,
    height = 5
)
for (j in 1:length(file_paths)){
    g <- ggplot(
        data = bind_rows(ms1_data_iso, ms1_data_mono) %>%
            filter(filename == samples[j]),
        aes(x = rt, y = int)
    ) +
        geom_line() +
        geom_vline(
            xintercept = r_t,
            linetype = "dashed",
            color = "red"
        ) +
        facet_wrap(~label, scales = "free_y", ncol = 1) +
        theme_minimal() +
        labs(
            title = paste0("Feature ID: ", feature_df$id, " monoisotopic mass: ", round(m_z_mono, 4), " isotopologue type: ", feature_df$isotopologue_type, " mode: ", feature_df$mode),
            subtitle = paste0("Sample: ", samples[j]),
            x = "Retention Time (min)",
            y = "Intensity"
        ) 
    # add to pdf
    print(g)

    }
# end pdf
dev.off()
}