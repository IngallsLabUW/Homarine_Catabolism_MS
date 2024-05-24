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
save_dir <- here("data", "intermediate", "metabolomics", "obi1", "targeted")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# Fill in blanks with 0s for dissolved data ----
## HILIC Neg -----
dat_filename <- here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_Neg_Dissolved_OscarSosaBacteria.csv"
)

dat <- read_csv(dat_filename,
  show_col_types = FALSE
) %>%
  mutate(Area = ifelse(
    str_detect(`Replicate Name`, "Blk") & Area == 0,
    "Test",
    Area
  ))
write_csv(dat, here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_Neg_Dissolved_OscarSosaBacteria_zerod.csv"
))

## HILIC Pos -----
dat_filename <- here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_Pos_Dissolved_OscarSosaBacteria.csv"
)
dat <- read_csv(dat_filename) %>%
  mutate(Area = ifelse(
    str_detect(`Replicate Name`, "Blk") & Area == 0,
    "Test",
    Area
  ))
write_csv(dat, here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_Pos_Dissolved_OscarSosaBacteria_zerod.csv"
))


# QC ----
## Dissolved, HILIC Neg -----
### First run regular standards-----
QE_QC(
  dat_filename = here(
    outer_file_loc,
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
    save_dir,
    "QCd_HILIC_Neg_Dissolved.csv"
  )
)

### Second run, intermediates-----
QE_QC(
    dat_filename = here(
        outer_file_loc,
        "dissolved",
        "skyline",
        "HILIC_NEG_DissolvedObi1_Intermediates.csv"
    ),
    std_flag <- NA,
    blank_flag <- "MQBlk",
    sn_min = 3,
    ppm_flex = 6,
    area_min = 40000,
    rt_flex = 2.5,
    blank_ratio_max = 3,
    fileout = here(
        save_dir,
        "QCd_HILIC_Neg_Dissolved_Intermediates.csv"
    )
)


## Dissolved, HILIC Pos -----
### First run regular standards-----
QE_QC(
  dat_filename = here(
    outer_file_loc,
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
    save_dir,
    "QCd_HILIC_Pos_Dissolved.csv"
  )
)

### Second run, intermediates-----
QE_QC(
    dat_filename = here(
        outer_file_loc,
        "dissolved",
        "skyline",
        "HILIC_POS_DissolvedObi1_Intermediates.csv"
    ),
    std_flag <- NA,
    blank_flag <- "MQBlk",
    sn_min = 3,
    ppm_flex = 6,
    area_min = 40000,
    rt_flex = 2.5,
    blank_ratio_max = 3,
    fileout = here(
        save_dir,
        "QCd_HILIC_Pos_Dissolved_Intermediates.csv"
    )
)

## Particulate, HILIC Neg -----
### First run regular standards-----
QE_QC(
  dat_filename = here(
    outer_file_loc,
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
    save_dir,
    "QCd_HILIC_Neg_Particulate.csv"
  )
)

### Second run, intermediates-----
QE_QC(
    dat_filename = here(
        outer_file_loc,
        "particulate",
        "skyline",
        "HILIC_NEG_ParticulateOBi1_Intermediates.csv"
    ),
    std_flag <- NA,
    blank_flag <- "Blk",
    sn_min = 4,
    ppm_flex = 6,
    area_min = 40000,
    rt_flex = 2.5,
    blank_ratio_max = 3,
    fileout = here(
        save_dir,
        "QCd_HILIC_Neg_Particulate_Intermediates.csv"
    )
)


## Particulate, HILIC Pos -----
### First run regular standards-----
QE_QC(
  dat_filename = here(
    outer_file_loc,
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
    save_dir,
    "QCd_HILIC_Pos_Particulate.csv"
  )
)

### Second run, intermediates-----
QE_QC(
    dat_filename = here(
        outer_file_loc,
        "particulate",
        "skyline",
        "HILIC_POS_ParticulateOBi1_Intermediates.csv"
    ),
    std_flag <- NA,
    blank_flag <- "Blk",
    sn_min = 4,
    ppm_flex = 6,
    area_min = 40000,
    rt_flex = 2.5,
    blank_ratio_max = 3,
    fileout = here(
        save_dir,
        "QCd_HILIC_Pos_Particulate_Intermediates.csv"
    )
)

# Get mass feature data for all datasets ----
## Prepare standards lists -----
# from https://github.com/IngallsLabUW/Ingalls_Standards/blob/master/Ingalls_Lab_Standards.csv
standards <- read_csv(
  "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
  show_col_types = FALSE
) %>%
  clean_names() %>%
  filter(column == "HILIC") %>%
  select(compound_name, mz, z) %>%
  mutate(compound_name_lower = str_to_lower(compound_name)) %>%
  select(-compound_name) %>%
  distinct()

intermediate_standards <- read_csv(
  here("data", "intermediate", "metabolomics", "intermediates_mf_info.csv"),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  rename(z = polarity) %>%
  select(name, mz, z) %>%
  mutate(compound_name_lower = str_to_lower(name)) %>%
  select(-name) %>%
  distinct() %>%
  filter(!(compound_name_lower %in% standards$compound_name_lower))

standards <- bind_rows(standards, intermediate_standards)

## Grab data from all datasets -----
## Original run, HILICNeg -----
### Dissolved -----
dat1_filename <- here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_Neg_Dissolved_OscarSosaBacteria_zerod.csv"
)
dat1 <- read_csv(dat1_filename, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Dissolved") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = -1)

### Particulate -----
dat_filename2 <- here(
  outer_file_loc,
  "particulate",
  "skyline",
  "HILIC_Neg_Particulate_OscarSosaBacteria.csv"
)
dat2 <- read_csv(dat_filename2, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Particulate") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = -1)

## Original run, HILICPos -----
### Dissolved -----
dat_filename3 <- here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_Pos_Dissolved_OscarSosaBacteria.csv"
)
dat3 <- read_csv(dat_filename3, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Dissolved") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = 1)

### Particulate -----
dat_filename4 <- here(
  outer_file_loc,
  "particulate",
  "skyline",
  "HILIC_pos_Particulate_OscarSosaBacteria.csv"
)
dat4 <- read_csv(dat_filename4, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Particulate") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = 1)

## Intermediates run, HILICNeg -----
### Dissolved -----
dat5_filename <- here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_NEG_DissolvedObi1_Intermediates.csv"
)
dat5 <- read_csv(dat5_filename, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Dissolved") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = -1)

### Particulate -----
dat6_filename <- here(
  outer_file_loc,
  "particulate",
  "skyline",
  "HILIC_NEG_ParticulateObi1_Intermediates.csv"
)
dat6 <- read_csv(dat6_filename, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Particulate") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = -1)

## Intermediates run, HILICPos -----
### Dissolved -----
dat7_filename <- here(
  outer_file_loc,
  "dissolved",
  "skyline",
  "HILIC_POS_DissolvedObi1_Intermediates.csv"
)
dat7 <- read_csv(dat7_filename, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Dissolved") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = 1)

### Particulate -----
dat8_filename <- here(
  outer_file_loc,
  "particulate",
  "skyline",
  "HILIC_POS_ParticulateObi1_Intermediates.csv"
)
dat8 <- read_csv(dat8_filename, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(fract = "Particulate") %>%
  rename(compound_name = precursor_ion_name) %>%
  mutate(z = 1)

## Combine all data -----
dat_list <- list(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)
## select compound_name, retention_time, area, fract, z
## turn all retention_time, area, fract, z to numeric
dat_list <- lapply(dat_list, function(x) {
  x %>%
    select(compound_name, retention_time, area, fract, z) %>%
    mutate(retention_time = as.numeric(retention_time),
           area = as.numeric(area),
           fract = as.character(fract),
           z = as.numeric(z))
})
mf_info <- bind_rows(dat_list) %>%
  filter(area > 0) %>%
  select(compound_name, retention_time, fract, z) %>%
  group_by(compound_name, fract, z) %>%
  summarise(
    retention_time = mean(as.numeric(retention_time), na.rm = TRUE)
  ) %>%
  mutate(compound_name_lower = str_to_lower(compound_name))

mf_info_joined <- mf_info %>%
  left_join(standards, by = c("compound_name_lower", "z")) %>%
  rename(
    mass_feature = compound_name,
    sample_fraction = fract
  ) %>%
  select(mass_feature, retention_time, mz, z, sample_fraction)

write_csv(mf_info_joined, here(
  save_dir,
  "mf_info.csv"
))

# Combine the standard run and the intermediate run data
### Dissolved, HILIC Neg -----
qc_dat1 <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Neg_Dissolved.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat1_int <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Neg_Dissolved_Intermediates.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat1 <- bind_rows(qc_dat1, qc_dat1_int)
write_csv(qc_dat1, here(
    save_dir,
    "QCd_HILIC_Neg_Dissolved_combined.csv"
))

### Particulate, HILIC Neg -----
qc_dat2 <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Neg_Particulate.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat2_int <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Neg_Particulate_Intermediates.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat2 <- bind_rows(qc_dat2, qc_dat2_int)
write_csv(qc_dat2, here(
    save_dir,
    "QCd_HILIC_Neg_Particulate_combined.csv"
))

### Dissolved, HILIC Pos -----
qc_dat3 <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Pos_Dissolved.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat3_int <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Pos_Dissolved_Intermediates.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat3 <- bind_rows(qc_dat3, qc_dat3_int)
write_csv(qc_dat3, here(
    save_dir,
    "QCd_HILIC_Pos_Dissolved_combined.csv"
))

### Particulate, HILIC Pos -----
qc_dat4 <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Pos_Particulate.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat4_int <- read_csv(
    here(
        save_dir,
        "QCd_HILIC_Pos_Particulate_Intermediates.csv"
    ),
    skip = 1,
    show_col_types = FALSE
)
qc_dat4 <- bind_rows(qc_dat4, qc_dat4_int)
write_csv(qc_dat4, here(
    save_dir,
    "QCd_HILIC_Pos_Particulate_combined.csv"
))

