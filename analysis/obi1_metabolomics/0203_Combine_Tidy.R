# PURPOSE: combine and tidy the untargeted data
# 20240313  KRH moving from Obi1-specific repo

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(here)

# SET FILE LOCATIONS  -----
meta_data_loc <- here("data", "raw", "metabolomics", "obi1")
output_loc <- here("data", "intermediate", "metabolomics", "obi1", "untargeted")

## BMIS files ----
dat_diss_filename <- here(output_loc, "BMISd_HILIC_Dissolved.csv")
dat_part_filename <- here(output_loc, "BMISd_HILIC_Particulate.csv")

## Sample key -----
sample_key_filename <- here(meta_data_loc, "sample_key.csv")

## MF info files -----
dat_diss_MF_filename <- here(output_loc, "HILICdissolved_wide_adducts.csv")
dat_part_MF_filename <- here(output_loc, "HILICparticulate_wide_adducts.csv")

# READ FILES-----
dat_diss <- read_csv(dat_diss_filename, show_col_types = FALSE)
dat_part <- read_csv(dat_part_filename, show_col_types = FALSE)
sample_key <- read_csv(sample_key_filename, show_col_types = FALSE)
dat_diss_MF <- read_csv(dat_diss_MF_filename, show_col_types = FALSE)
dat_part_MF <- read_csv(dat_part_MF_filename, show_col_types = FALSE)


# MUNDGE DATA ----
## Combine data----
dat_combined <- dat_part %>%
  select(mass_feature, replicate_name, adjusted_area) %>%
  mutate(sample_fraction = "Particulate") %>%
  bind_rows(dat_diss %>%
    select(mass_feature, replicate_name, adjusted_area) %>%
    mutate(sample_fraction = "Dissolved"))

## Add meta data ---
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

## Add MF info ----
dat_combined3 <- dat_part_MF %>%
  select(ID:nm_theroetical) %>%
  mutate(sample_fraction = "Particulate") %>%
  bind_rows(dat_diss_MF %>%
    select(ID:nm_theroetical) %>%
    mutate(sample_fraction = "Dissolved")) %>%
  rename(mass_feature = ID) %>%
  semi_join(dat_combined2, by = join_by(mass_feature, sample_fraction)) %>%
    right_join(dat_combined2)

dat_combined_wide <- dat_combined2 %>%
  pivot_wider(
    id_cols = c("replicate_name", "sample_set", "treatment", "sample_fraction"),
    names_from = mass_feature,
    values_from = adjusted_area
  )


# Add MF info -----

# Write out long and wide to share
write_csv(
  dat_combined3,
  here(output_loc, "combined_tidy_dat_long.csv")
)
write_csv(dat_combined_wide, here(output_loc, "combined_tidy_dat_wide.csv"))
