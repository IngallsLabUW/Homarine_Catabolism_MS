# _____________________________________________________________________________
#       Copyright 2022, Integral Consulting Inc. All rights reserved.        #
# _____________________________________________________________________________
#
# _____________________________________________________________________________
#  0104_Combine_Tidy.R
# _____________________________________________________________________________
# PURPOSE:
#  Combine particulate, dissolved data
#
#
# PROJECT INFORMATION:
#   Name:
#   Number: C3352
#
# HISTORY:
# 	 Date		        Remarks
# _____________________________________________________________________________
# 	 20220526       Started script
#

# PACKAGES ----
library(tidyverse)
library(here)

# Read in file  -----
dat_diss_filename <- "01_targeted/data_intermediate/BMISd_HILIC_Dissolved_postDissBlk.csv"
dat_part_filename <- "01_targeted/data_intermediate/BMISd_HILIC_Particulate.csv"
sample_key_filename <- "raw_input/sample_key.csv"

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
  left_join(sample_key %>%
              select(replicate_name,
                     sample_set,
                     treatment),
            by = "replicate_name") %>%
  filter(!is.na(treatment))

dat_combined_wide <- dat_combined2 %>%
  pivot_wider(id_cols = c('replicate_name', 'sample_set', 'treatment','sample_fraction'),
              names_from = mass_feature,
              values_from = adjusted_area)

# Write out long and wide to share
write_csv(dat_combined2, "01_targeted/data_intermediate/combined_tidy_dat_long.csv")
write_csv(dat_combined_wide, "01_targeted/data_intermediate/combined_tidy_dat_wide.csv")

