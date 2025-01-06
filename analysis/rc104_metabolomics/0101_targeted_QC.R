# Purpose: Quality control on all targeted datasets for the RC104 Metabolomics

# History:
# Date	     Remarks
# 20250106   K Heal changing to work with RC104 data

library(tidyverse)
library(janitor)
library(here)
source(here("analysis", "obi1_metabolomics", "functions", "targeted_dataprocessing_fxns.R"))

# Set location of skyline files here
outer_file_loc <- here("data", "raw", "metabolomics", "RC104", "particulate", "mzml")
save_dir <- here("data", "intermediate", "metabolomics", "rc104", "targeted")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# QC ----
## HILIC Neg -----
QE_QC(
  dat_filename = here(
    outer_file_loc,
    "HILIC_negative",
    "HILIC_QE_NEG_Dinimite2_HomarineFate_Report.csv"
  ),
  std_flag <- NA,
  blank_flag <- "_Blk_Blk_",
  sn_min = 3,
  ppm_flex = 6,
  area_min = 40000,
  rt_flex = 2.5,
  blank_ratio_max = 2,
  fileout = here(
    save_dir,
    "QCd_HILIC_Neg.csv"
  )
)

## HILIC Pos -----
QE_QC(
  dat_filename = here(
    outer_file_loc,
    "HILIC_positive",
    "HILIC_QE_POS_Dinimite2_HomarineFate_Report.csv"
  ),
  std_flag <- NA,
  blank_flag <- "_Blk_Blk_",
  sn_min = 3,
  ppm_flex = 6,
  area_min = 40000,
  rt_flex = 2.5,
  blank_ratio_max = 2,
  fileout = here(
    save_dir,
    "QCd_HILIC_Pos.csv"
  )
)

# Get mass feature data for all datasets ----
## Download standards list from https://github.com/IngallsLabUW/Ingalls_Standards/blob/master/Ingalls_Lab_Standards.csv
standards <- read_csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv", show_col_types = FALSE) %>%
  clean_names() %>%
  filter(column == "HILIC") %>%
  select(compound_name, mz, z) %>%
  mutate(compound_name_lower = str_to_lower(compound_name)) %>%
  select(-compound_name) %>%
  distinct()

## HILICNeg -----
dat_filename <- here(
    save_dir,
  "QCd_HILIC_Neg.csv"
)
dat <- read_csv(dat_filename, show_col_types = FALSE, skip = 1) %>%
  clean_names() %>%
  mutate(z = -1)
dat_filename2 <- here(
    save_dir,
  "QCd_HILIC_Pos.csv"
)
dat2 <- read_csv(dat_filename2, show_col_types = FALSE,  skip = 1) %>%
  clean_names() %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = 1)

mf_info <- dat %>%
  bind_rows(dat2) %>%
  filter(area > 0) %>%
  select(compound_name, retention_time, z) %>%
  group_by(compound_name, z) %>%
  summarise(
    retention_time = mean(as.numeric(retention_time), na.rm = TRUE)
  ) %>%
  mutate(compound_name_lower = str_to_lower(compound_name))

mf_info_joined <- mf_info %>%
  left_join(standards, by = c("compound_name_lower", "z")) %>%
  rename(
    mass_feature = compound_name
    ) %>%
  select(mass_feature, retention_time, mz, z)

write_csv(mf_info_joined, here(
  save_dir,
  "mf_info.csv"
))
