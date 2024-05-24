# Purpose:  Combine particulate, dissolved data for targeted data

# History:
# Date	     Remarks
# 20240415   K Heal started script for Rpom, using Obi1 as template

# PACKAGES ----
library(tidyverse)
library(here)

# Read in file  -----
qc_data_loc <- here("data", "intermediate", "metabolomics", "g4", "targeted")
meta_data_loc <- here("data", "raw", "metabolomics", "g4")

dat_part_filename <- here(qc_data_loc, "BMISd_HILIC_Particulate.csv")
sample_key_filename <- here(meta_data_loc, "sample_key.csv")

dat_part <- read_csv(dat_part_filename, show_col_types = FALSE)
sample_key <- read_csv(sample_key_filename, show_col_types = FALSE)

# Combine all data ----
dat_combined <- dat_part %>%
  select(mass_feature, replicate_name, adjusted_area, fraction) %>%
  mutate(sample_fraction = "Particulate") %>%
  mutate(z = case_when(
    str_detect(fraction, "Neg") ~ -1,
    str_detect(fraction, "Pos") ~ 1
  )) %>%
  select(-fraction)

# Add meta data ---
dat_combined2 <- dat_combined %>%
  left_join(
    sample_key %>%
      select(
        replicate_name,
        treatment
      ),
    by = "replicate_name"
  ) %>%
  filter(!is.na(treatment),
         !str_detect(treatment, "blank"),
         !str_detect(treatment, "media"),
         !str_detect(treatment, "in ASW"))

# Add mf data to long data
mf_data <- read_csv(
  here(qc_data_loc, "mf_info.csv"),
  show_col_types = FALSE
)
dat_combined3 <- dat_combined2 %>%
  left_join(mf_data, by = join_by(mass_feature, z)) %>%
  mutate(timepoint = str_extract(replicate_name, "_t\\d+") %>%
             str_replace("_t", "") %>%
             as.numeric()) %>%
    mutate(experiment_id = str_extract(replicate_name, "G4_\\d+")) %>%
    mutate(treatment = case_when(
        str_detect(replicate_name, "Poo") ~ NA_character_,
      str_detect(treatment, "MethylHomFate") ~ "Homarine",
      str_detect(treatment, "Control") ~ "Control"
    ))

glimpse(dat_combined3)
# pivot wider
dat_combined_wide <- dat_combined3 %>%
  pivot_wider(
    id_cols = c("replicate_name", "treatment", "sample_fraction", "timepoint", "experiment_id"),
    names_from = mass_feature,
    values_from = adjusted_area
  )

# Write out long and wide to share
write_csv(dat_combined3, here(qc_data_loc, "combined_tidy_dat_long.csv"))
write_csv(dat_combined_wide, here(qc_data_loc, "combined_tidy_dat_wide.csv"))
