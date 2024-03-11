# Purpose: Quality control on all targeted datasets for the Obi1 Metabolomics

# History:
# Date	     Remarks
# 20240311   K Heal transferring script from Obi1 specific repo to general repo

library(tidyverse)
library(janitor)
library(here)
source(here("analysis", "obi1_metabolomics", "functions", "targeted_dataprocessing_fxns.R"))

# Set location of skyline files here
outer_file_loc <- here("data", "raw", "metabolomics", "obi1")

# Fill in blanks with 0s for dissolved data ----
## HILIC Neg -----
dat_filename = here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_Neg_Dissolved_OscarSosaBacteria.csv"
)

dat <- read_csv(dat_filename) %>%
  mutate(Area = ifelse(
    str_detect(`Replicate Name`, "Blk") & Area == 0,
    "Test",
    Area
      )
    )
write_csv(dat, here(
  "raw_input",
  "dissolved",
  "skyline",
  "HILIC_Neg_Dissolved_OscarSosaBacteria_zerod.csv"
))

# Fill in blanks with 0s for dissolved data ----
## HILIC Pos -----
dat_filename = here(
  "raw_input",
  "dissolved",
  "skyline",
  "HILIC_Pos_Dissolved_OscarSosaBacteria.csv"
)
dat <- read_csv(dat_filename) %>%
  mutate(Area = ifelse(
    str_detect(`Replicate Name`, "Blk") & Area == 0,
    "Test",
    Area
  )
  )
write_csv(dat,  here(
  "raw_input",
  "dissolved",
  "skyline",
  "HILIC_Pos_Dissolved_OscarSosaBacteria_zerod.csv"
))


# QC on Dissolved, HILIC Neg -----
QE_QC(
  dat_filename = here(
    "raw_input",
    "dissolved",
    "skyline",
    "HILIC_Neg_Dissolved_OscarSosaBacteria_zerod.csv"
  ),
  std_flag <- "Std_4",
  blank_flag <- "MQBlk",
  sn_min = 3,
  ppm_flex = 6,
  area_min = 40000,
  rt_flex = 2.5,
  blank_ratio_max = 3,
  fileout = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Neg_Dissolved.csv"
  )
)
# QC on Dissolved, HILIC Pos -----
QE_QC(
  dat_filename = here(
    "raw_input",
    "dissolved",
    "skyline",
    "HILIC_Pos_Dissolved_OscarSosaBacteria.csv"
  ),
  std_flag <- "Std_4",
  blank_flag <- "MQBlk",
  sn_min = 3,
  ppm_flex = 6,
  area_min = 40000,
  rt_flex = 2.5,
  blank_ratio_max = 3,
  fileout = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Pos_Dissolved.csv"
  )
)

# QC on Particulate, HILIC Neg -----
QE_QC(
  dat_filename = here(
    "raw_input",
    "particulate",
    "skyline",
    "HILIC_Neg_Particulate_OscarSosaBacteria.csv"
  ),
  std_flag <- "Std_4",
  blank_flag <- "Blk",
  sn_min = 4,
  ppm_flex = 6,
  area_min = 40000,
  rt_flex = 2.5,
  blank_ratio_max = 3,
  fileout = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Neg_Particulate.csv"
  )
)

# QC on Particulate, HILIC Pos -----
QE_QC(
  dat_filename = here(
    "raw_input",
    "particulate",
    "skyline",
    "HILIC_Pos_Particulate_OscarSosaBacteria.csv"
  ),
  std_flag <- "Std_4",
  blank_flag <- "Blk",
  sn_min = 4,
  ppm_flex = 6,
  area_min = 40000,
  rt_flex = 2.5,
  blank_ratio_max = 3,
  fileout = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Pos_Particulate.csv"
  )
)
