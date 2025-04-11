# This will make a table of all the adjusted areas and calculated enrichment factors for all except Obi1
# metabolomics data

# Load required libraries
library(tidyverse)
library(here)

# Get metadata file
metadat <- read_csv(here("tables", "metabolomics_experimental_metadata.csv"))

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

# grab mass feature data
mf_info <- read_csv(
    here("data", "intermediate", "metabolomics", "mf_info_curated.csv"),
    show_col_types = FALSE) 

core_metabs <- read_csv(
    here("data", "intermediate", "metabolomics", "coremetabs_working.csv"),
    show_col_types = FALSE) 

## prepare metadata for each experiment -------
g4_mf_info <- read_csv(
    here(inter_dir_list[["g4"]], "targeted/mf_info.csv"),
    show_col_types = FALSE)
g4_area_data <- read_csv(
    here(inter_dir_list[["g4"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    filter(replicate_name %in% metadat$sample_id) %>%
    select(replicate_name, mass_feature, adjusted_area) %>%
    distinct() 

g5_mf_info <- read_csv(
    here(inter_dir_list[["g5"]], "targeted/mf_info.csv"),
    show_col_types = FALSE)
g5_area_data <- read_csv(
    here(inter_dir_list[["g5"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    filter(replicate_name %in% metadat$sample_id) %>%
    select(replicate_name, mass_feature, adjusted_area) %>%
    distinct()

rc104_mf_info <- read_csv(
    here(inter_dir_list[["rc104"]], "targeted/mf_info.csv"),
    show_col_types = FALSE)
rc104_area_data <- read_csv(
    here(inter_dir_list[["rc104"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    filter(replicate_name %in% metadat$sample_id) %>%
    select(replicate_name, mass_feature, adjusted_area) %>%
    distinct()

rpom_mf_info <- read_csv(
    here(inter_dir_list[["rpom"]], "targeted/mf_info.csv"),
    show_col_types = FALSE)
rpom_area_data <- read_csv(
    here(inter_dir_list[["rpom"]], "targeted/combined_tidy_dat_long.csv"),
    show_col_types = FALSE) %>%
    filter(replicate_name %in% metadat$sample_id) %>%
    select(replicate_name, mass_feature, adjusted_area) %>%
    distinct()

#Combine all area data
area_data <- bind_rows(g4_area_data, g5_area_data, rc104_area_data, rpom_area_data) %>%
    pivot_wider(names_from = replicate_name, values_from = adjusted_area)


# Prepare final data frame --------
final_df <- mf_info %>%
    filter(is.na(remove)) %>%
    filter(is.na(internal_standard)) %>%
    filter(stable_isotope_version != "3H2") %>%
    select(mass_feature, molecule_id, stable_isotope_version, retention_time, mz, z, core) %>%
    left_join(area_data, by = "mass_feature") %>%
    arrange(core, molecule_id, stable_isotope_version) %>%
    mutate(retention_time = round(retention_time, 2),
           mz = round(mz, 4))

# Write final data frame to CSV --------
write_csv(
    final_df,
    here("tables", "metabolomics_experiments_targeted_results.csv"),
    na = ""
)

