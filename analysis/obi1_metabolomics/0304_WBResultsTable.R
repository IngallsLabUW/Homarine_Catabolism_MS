# PURPOSE: provide a summary table of the data for metabolomics workbench submission


# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(here)
library(openxlsx2)


# SET FILE LOCATIONS  -----
meta_data_dir <- here("data", "intermediate", "metabolomics")
output_loc <- here("data", "intermediate", "metabolomics", "obi1")
sample_key <- here("data", "raw", "metabolomics", "obi1", "sample_key.csv")
sample_fraction_oi <- "Particulate"
save_loc <- here("tables")

# IMPORT DATA ----
raw_data <- here(output_loc, "combined_long_dat.csv") %>%
  read_csv(show_col_types = FALSE)
enrichment_data <- here(output_loc, "enrichment_summary.csv") %>%
  read_csv(show_col_types = FALSE)
mf_info <- here(output_loc, "combined_MF_info.csv") %>%
  read_csv(show_col_types = FALSE)
sample_key <- read_csv(sample_key, show_col_types = FALSE)
ingalls_standards <- read_csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv", show_col_types = FALSE) %>%
    clean_names() %>%
    filter(column == "HILIC") %>%
    rename(mass_feature = compound_name) %>%
    select(mass_feature, mz, rt_minute, z, ionization_form, kegg_code, pub_chem_code) %>%
    # Change N-Methyl-L-glutamic acid to n-methyl glutamic acid so they match
    mutate(mass_feature = ifelse(mass_feature == "N-Methyl-L-glutamic acid", "n-methyl glutamic acid", mass_feature)) %>%
    distinct()

supp_t <- raw_data %>%
  filter(sample_set == "OBi1_set2") %>%
  select(mass_feature, replicate_name, adjusted_area) %>%
  pivot_wider(names_from = replicate_name, values_from = adjusted_area) %>%
  # Remove mass features that start with "HILIC" - these are unknowns
  filter(str_starts(mass_feature, "HILIC") == FALSE) %>%
  filter(str_starts(mass_feature, "unknown") == FALSE)

# Add in mass feature info
supp_t2 <- supp_t %>%
  left_join(mf_info %>%
    filter(sample_fraction == sample_fraction_oi), by = "mass_feature") %>%
  select(mass_feature, z, starts_with("2")) %>%
    filter(mass_feature %in% ingalls_standards$mass_feature) %>%
    filter(!is.na(z))

# Prepare mass feature metadata
mf_info <- supp_t2 %>%
  select(mass_feature, z) %>%
  distinct() %>%
left_join(ingalls_standards, by = c("mass_feature", "z"))
  

# Separate results and mass feature metadata into positive and negative mode
t2_pos <- supp_t2 %>%
  filter(z > 0)  %>%
  select(mass_feature,starts_with("2"))
t2_neg <- supp_t2 %>%
  filter(z < 0)  %>%
  select(mass_feature,starts_with("2"))
mf_info_pos <- mf_info %>%
  filter(z > 0) %>%
    select(-z)
mf_info_neg <- mf_info %>%
  filter(z < 0) %>%
    select(-z)

write_tsv(t2_pos, here(save_loc, "WB_OBi1_Positive_MF_Table.tsv"))
write_tsv(t2_neg, here(save_loc, "WB_OBi1_Negative_MF_Table.tsv"))
write_tsv(mf_info_pos, here(save_loc, "WB_OBi1_Positive_MF_Metadata.tsv"))
write_tsv(mf_info_neg, here(save_loc, "WB_OBi1_Negative_MF_Metadata.tsv"))
