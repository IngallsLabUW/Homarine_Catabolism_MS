# This will make a table of all the experimental metadata for the metabolomics experiments

# Load required libraries
library(tidyverse)
library(here)

# Get cruise info
cruise_info <- read_csv(here("figures", "maps", "2025_1_14_hom_fate_coords_map.csv"),
                        show_col_types = FALSE)

# Set up data directories --------
## Sample key files are in the raw data directories --------
sub_dirs <- c("rpom", "obi1", "g4", "g5", "rc104")
raw_dir_list <- lapply(sub_dirs, function(x) {
    here("data", "raw", "metabolomics", x)
})
names(raw_dir_list) <- sub_dirs

## Enrichment results are in the intermediate data directories --------
inter_dir_list <- lapply(raw_dir_list, function(x) {
    str_replace(x, "raw", "intermediate")
})

# Grab sample_key.csv from each directory  --------
sample_key_list <- lapply(raw_dir_list, function(x) {
  sample_key_path <- file.path(x, "sample_key.csv")
  read_csv(sample_key_path,
             show_col_types = FALSE)
  })
names(sample_key_list) <- names(raw_dir_list)

# Prepare sample key for g4
cruise_short_id = "g4"

## prepare metadata for each experiment -------
g4_metadata <- read_csv(
    here(inter_dir_list[[cruise_short_id]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    select(replicate_name, sample_fraction, treatment, timepoint, experiment_id) %>%
    distinct() %>%
    left_join(cruise_info, by = join_by(experiment_id)) %>%
    select(cruise, experiment_id, station, lat, lon, date, replicate_name, sample_fraction, treatment, timepoint) %>%
    mutate(treatment = case_when(
        treatment == "Homarine" ~ "Homarine (2H3-labelled)",
        TRUE ~ treatment)) %>%
    filter(!is.na(cruise))

g5_metadata <- read_csv(
    here(inter_dir_list[["g5"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    select(replicate_name, sample_fraction, treatment, timepoint, experiment_id) %>%
    distinct() %>%
    left_join(cruise_info, by = join_by(experiment_id)) %>%
    select(cruise, experiment_id, station, lat, lon, date, replicate_name, sample_fraction, treatment, timepoint) %>%
    # change Homarine treatment to + D3-labelled homarine
    mutate(treatment = case_when(
        treatment == "Homarine" ~ "Homarine (13C7,15N-labelled)",
        TRUE ~ treatment)) %>%
    filter(!is.na(cruise))

rc104_metadata <- read_csv(
    here(inter_dir_list[["rc104"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    select(replicate_name, sample_fraction, treatment, timepoint, experiment_id) %>%
    mutate(experiment_id = "rc104") %>%
    distinct() %>%
    left_join(cruise_info, by = join_by(experiment_id)) %>%
    select(cruise, experiment_id, station, lat, lon, date, replicate_name, sample_fraction, treatment, timepoint) %>%
    mutate(treatment = case_when(
        treatment == "Homarine" ~ "Homarine (13C7,15N-labelled)",
        TRUE ~ treatment)) %>%
    filter(!is.na(cruise))

obi1_metadata <- read_csv(
    here(inter_dir_list[["obi1"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    filter(sample_set == "OBi1_set2") %>%
    select(replicate_name, sample_fraction, treatment) %>%
    mutate(isolate_id = "OBi1") %>%
    mutate(treatment = case_when(
        treatment == "Glucose + NH4" ~ "Glucose",
        treatment == "Glucose + NH4 + Homarine" ~ "Glucose + homarine",
        TRUE ~ treatment)) %>%
    distinct()
    
rpom_metadata <- read_csv(
    here(inter_dir_list[["rpom"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    select(replicate_name, sample_fraction, treatment) %>%
    mutate(isolate_id = "Rpom DSS3") %>%
    filter(treatment %in% c("Glucose", "Homarine", "Glucose and homarine")) %>%
    mutate(treatment = case_when(
        treatment == "Glucose and homarine" ~ "Glucose + homarine",
        TRUE ~ treatment)) %>%
    filter(sample_fraction == "Particulate") %>%
    distinct()

# Combine all metadata into one dataframe
dat_cmb <- obi1_metadata %>%
    bind_rows(rpom_metadata) %>%
    bind_rows(g4_metadata) %>%
    bind_rows(g5_metadata) %>%
    bind_rows(rc104_metadata) %>%
    rename(sample_id = replicate_name) %>%
    select(sample_id, sample_fraction, isolate_id, cruise, experiment_id, treatment, timepoint, station, lat, lon, date)

# Save the combined metadata as an xcel file
write_csv(dat_cmb, here("tables", "metabolomics_experimental_metadata.csv"))

