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
analytical_batch_ids <- c("DSS3", "OBi1", "TN397", "TN412", "RC104")
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
g5_coremetabs <- read_csv(
    here(inter_dir_list[["g5"]], "core_metabs_quality.csv"),
    show_col_types = FALSE) %>%
    pull(core_metab)
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
## Prepare column of core metab info --------
core_info <- mf_info %>%
    filter(core == "core") %>%
    select(mass_feature) %>%
    mutate(studies_used = "")

# For each of the studies, add its study ID to the studies_used column if the mass feature is a core metab in that study
for (i in 1:length(sub_dirs)) {
    study_id <- analytical_batch_ids[i]
    exp_dir <- raw_dir_list[[i]]
    core_metabs_file <- here(inter_dir_list[[i]], "core_metabs_quality.csv")
    study_core_metabs <- read_csv(core_metabs_file, show_col_types = FALSE) %>%
        pull(core_metab)
    core_info <- core_info %>%
            mutate(studies_used = if_else(
                mass_feature %in% study_core_metabs,
                paste0(studies_used, study_id, "; "),
                studies_used
            ))
}
# Remove trailing semicolon and space
core_info <- core_info %>%
    mutate(studies_used = str_remove(studies_used, "; $"))

final_df <- mf_info %>%
    filter(is.na(remove)) %>%
    filter(is.na(internal_standard)) %>%
    filter(stable_isotope_version != "3H2" | is.na(stable_isotope_version)) %>%
    select(mass_feature, molecule_id, stable_isotope_version, retention_time, mz, z, core) %>%
    left_join(area_data, by = "mass_feature") %>%
    arrange(core, molecule_id, stable_isotope_version) %>%
    mutate(retention_time = round(retention_time, 2),
           mz = round(mz, 4))

# Add column that notes if this core metab was used in G5 calculations
final_df <- final_df %>%
    left_join(
        core_info %>%
            select(mass_feature, studies_used),
        by = "mass_feature"
    ) %>%
    select(mass_feature, molecule_id, stable_isotope_version, retention_time, mz, z, core, studies_used, everything()) %>%
    filter(studies_used != "" | core != "core")

# Write final data frame to CSV --------
write_csv(
    final_df,
    here("tables", "metabolomics_experiments_targeted_results.csv"),
    na = ""
)

