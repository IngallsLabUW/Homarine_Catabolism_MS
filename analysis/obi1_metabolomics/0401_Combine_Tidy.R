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
library(viridis)
library(broom)
library(coin)


# Read in files  -----
dat_target_filename <-
  "01_targeted/data_intermediate/combined_tidy_dat_long.csv"
dat_nontarget_filename <-
  "02_untargeted/intermediates/combined_tidy_dat_long.csv"

dat_target <- read_csv(dat_target_filename, show_col_types = FALSE)
dat_nontarget <-
  read_csv(dat_nontarget_filename, show_col_types = FALSE)


# MUNGE DATA ----
## Combine data -----
dat_cmb <- dat_target %>%
  mutate(dat_type = "targeted") %>%
  bind_rows(dat_nontarget %>%
              mutate(dat_type = "nontargeted")) %>%
  filter(sample_set == "OBi1_set2") %>%
  filter(sample_fraction == "Particulate")

## Calculate fold changes and select median fold change for scaling ------
dat_fc <- dat_cmb %>%
  group_by(mass_feature, treatment, sample_fraction) %>%
  # dump molecules that were not observed in all treatments
  filter(!is.na(adjusted_area)) %>%
  filter(n_distinct(adjusted_area) > 2) %>%
  summarise(mean_area = mean(adjusted_area, rm.na = T)) %>%
  ungroup() %>%
  group_by(mass_feature, sample_fraction) %>%
  summarise(fc_hNH4 = mean_area[treatment == "Glucose + NH4 + Homarine"] /
              mean_area[treatment == "Glucose + NH4"],
            fc_h = mean_area[treatment == "Homarine"] / mean_area[treatment == "Glucose + NH4"]) %>%
  ungroup()

scaler_h <- dat_fc %>%
  group_by(sample_fraction) %>%
  mutate(diff_from_med_h = abs(fc_h-median(fc_h, na.rm = T))) %>%
  filter(diff_from_med_h == min(diff_from_med_h, na.rm = T)) %>%
  ungroup() %>%
  head(1) %>%
  pull(fc_h)

scaler_hNH4 <- dat_fc %>%
  group_by(sample_fraction) %>%
  mutate(diff_from_med_hNH4 = abs(fc_hNH4-median(fc_hNH4, na.rm = T))) %>%
  filter(diff_from_med_hNH4 == min(diff_from_med_hNH4, na.rm = T))%>%
  ungroup() %>%
  head(1) %>%
  pull(fc_hNH4)

## scale data ----
dat_cmb2 <- dat_cmb %>%
  mutate(scaled_area = case_when(
    treatment == "Homarine" ~ adjusted_area/scaler_h,
    treatment == "Glucose + NH4 + Homarine" ~ adjusted_area/scaler_hNH4,
    treatment == "Glucose + NH4"  ~ adjusted_area
  ))


## summarize abundance features ----
high_homarine_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area))%>%
  filter(treatment != "Glucose + NH4 + Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area)))%>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Homarine" & area_rank_sum == 6) %>%
  pull(mass_feature)

high_homarine_NH4_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area))%>%
  filter(treatment != "Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area)))%>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Glucose + NH4 + Homarine" & area_rank_sum == 6) %>%
  pull(mass_feature)

low_homarine_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area))%>%
  filter(treatment != "Glucose + NH4 + Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area)))%>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Homarine" & area_rank_sum == 15) %>%
  pull(mass_feature)

low_homarine_NH4_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area))%>%
  filter(treatment != "Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area)))%>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Glucose + NH4 + Homarine" & area_rank_sum == 15) %>%
  pull(mass_feature)

high_both_compounds <- high_homarine_compounds[high_homarine_compounds %in% high_homarine_NH4_compounds]
low_both_compounds <- low_homarine_compounds[low_homarine_compounds %in% low_homarine_NH4_compounds]

dat_cmb2 <- dat_cmb2 %>%
  mutate(high_homarine_compound = mass_feature %in% high_homarine_compounds,
         low_homarine_compound = mass_feature %in% low_homarine_compounds,
         high_homarine_NH4_compound = mass_feature %in% high_homarine_NH4_compounds,
         low_homarine_NH4_compound = mass_feature %in% low_homarine_NH4_compounds,
         high_both_compound = mass_feature %in% high_both_compounds,
         low_both_compounds = mass_feature %in% low_both_compounds)

## make intermediate of MF info -----
dat_MF_2 <- dat_cmb2 %>%
  select(mass_feature, sample_fraction, sample_set, dat_type, RT, mz, MS2,
         add_annotation, high_homarine_compound:low_both_compounds) %>%
  distinct()

# WRITE DATA -----
write_csv(dat_cmb2, "04_combined/intermediates/combined_long_dat_scaled.csv")
write_csv(dat_MF_2, "04_combined/intermediates/combined_MF_info.csv")
