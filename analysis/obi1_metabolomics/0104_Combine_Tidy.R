# Purpose:  Combine particulate, dissolved data for targeted data

# History:
# Date	     Remarks
# 20240312   K Heal transferring script from Obi1 specific repo to general repo
# 20220526   K Heal started script

# PACKAGES ----
library(tidyverse)
library(here)

# Read in file  -----
qc_data_loc <- here("data", "intermediate", "metabolomics", "obi1", "targeted")
meta_data_loc <- here("data", "raw", "metabolomics", "obi1")

dat_diss_filename <- here(qc_data_loc, "BMISd_HILIC_Dissolved_postDissBlk.csv")
dat_part_filename <- here(qc_data_loc, "BMISd_HILIC_Particulate.csv")
sample_key_filename <- here(meta_data_loc, "sample_key.csv")

dat_diss <- read_csv(dat_diss_filename, show_col_types = FALSE)
dat_part <- read_csv(dat_part_filename, show_col_types = FALSE)
sample_key <- read_csv(sample_key_filename, show_col_types = FALSE)

# Combine all data ----
dat_combined <- dat_part %>%
  select(mass_feature, replicate_name, adjusted_area) %>%
  mutate(sample_fraction = "Particulate") %>%
  bind_rows(dat_diss %>%
    select(mass_feature, replicate_name, adjusted_area) %>%
    mutate(sample_fraction = "Dissolved"))

# Add meta data ---
dat_combined2 <- dat_combined %>%
  left_join(
    sample_key %>%
      select(
        replicate_name,
        sample_set,
        treatment
      ),
    by = "replicate_name"
  ) %>%
  filter(!is.na(treatment))

dat_combined_wide <- dat_combined2 %>%
  pivot_wider(
    id_cols = c("replicate_name", "sample_set", "treatment", "sample_fraction"),
    names_from = mass_feature,
    values_from = adjusted_area
  )

# Write out long and wide to share
write_csv(dat_combined2, here(qc_data_loc, "combined_tidy_dat_long.csv"))
write_csv(dat_combined_wide, here(qc_data_loc, "combined_tidy_dat_wide.csv"))
