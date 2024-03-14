#
# Copyright: Copyright 2022, Integral Consulting Inc. All rights reserved.
#
# Purpose:
#
# Project Information:
#   Name:
#   Number:
#
# History:
# Date	     Remarks
# YYYY-MM-DD
#

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(integral)
library(janitor)

# SET FILE LOCATIONS  -----
## BMIS files ----
dat_diss_filename <- "02_untargeted/intermediates/BMISd_HILIC_Dissolved.csv"
dat_part_filename <- "02_untargeted/intermediates/BMISd_HILIC_Particulate.csv"

## Sample key -----
sample_key_filename <- "raw_input/sample_key.csv"

## MF info files -----
dat_diss_MF_filename <- "02_untargeted/intermediates/HILICdissolved_wide_adducts.csv"
dat_part_MF_filename <- "02_untargeted/intermediates/HILICparticulate_wide_adducts.csv"

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
  select(ID:add_annotation) %>%
  mutate(sample_fraction = "Particulate") %>%
  bind_rows(dat_diss_MF %>%
    select(ID:add_annotation) %>%
    mutate(sample_fraction = "Dissolved")) %>%
  rename(mass_feature = ID) %>%
  filter(mass_feature %in% dat_combined2$mass_feature) %>%
  right_join(dat_combined2)


dat_combined_wide <- dat_combined2 %>%
  pivot_wider(
    id_cols = c("replicate_name", "sample_set", "treatment", "sample_fraction"),
    names_from = mass_feature,
    values_from = adjusted_area
  )


# Add MF info -----

# Write out long and wide to share
write_csv(dat_combined3, "02_untargeted/intermediates/combined_tidy_dat_long.csv")
write_csv(dat_combined_wide, "02_untargeted/intermediates/combined_tidy_dat_wide.csv")
