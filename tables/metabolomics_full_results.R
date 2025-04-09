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
    here("data", "intermediate", "metabolomics", "intermediates_mf_info2.csv"),
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
mf_data <- area_data %>%
    select(mass_feature) %>%
    left_join( bind_rows(g4_mf_info, g5_mf_info, rc104_mf_info, rpom_mf_info), by = "mass_feature") %>%
    filter(mass_feature %in% area_data$mass_feature) %>%
    group_by(mass_feature) %>%
    summarise(
        retention_time = mean(retention_time, na.rm = TRUE),
        mz = mean(mz, na.rm = TRUE),
        z = mean(z, na.rm = TRUE)
    ) %>%
    ungroup()
## Fill in missing m/z from mf_data
mf_data2 <- mf_data %>%
    rows_patch(
        mf_info %>% 
            rename(mass_feature = name) %>%
            select(mass_feature, mz) %>%
            distinct(),
        by = c("mass_feature"),
        unmatched = "ignore",
    ) %>%
    mutate(core = ifelse(mass_feature %in% core_metabs$core_metabolites_fMsystem, "core", "non-core"))
write_csv(mf_data2, here("data", "intermediate", "metabolomics", "mf_info_curated.csv"))
