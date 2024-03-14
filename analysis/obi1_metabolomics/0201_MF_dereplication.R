# Purpose: Dereplicate mass features from MSDial output
# History:
# Date	     Remarks
# 2022-10-25  K Heal starting script

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(Hmisc)
source(here("analysis", "obi1_metabolomics", "functions", "untargeted_dataprocessing_fxns.R"))

# SET FILE LOCS ----
meta_data_loc <- here("data", "raw", "metabolomics", "obi1")
adduct_list_filename <- here("data", "raw", "metabolomics", "adduct_list.csv")
output_loc <- here("data", "intermediate", "metabolomics", "obi1", "untargeted")

# HILIC, particulate ------
## Find adducts, clean up output -----
dat2 <- dereplicate_MFs(
  dat_filename1 = here(
    meta_data_loc,
    "particulate",
    "msdial",
    "HILIC_positive",
    "Area_2_20221018113.txt"
  ),
  dat_type_1 = "HILICPos",
  dat_filename2 = here(
    meta_data_loc,
    "particulate",
    "msdial",
    "HILIC_negative",
    "Area_2_202210181530.txt"
  ),
  dat_type_2 = "HILICNeg",
  comment_tag = "m",
  sample_tag = "2112",
  drop_tags = c("DDA", "_Blk_"),
  RT_tol = 0.3, # minutes to allow for possible adduct/isotope
  cor_tol = 0.85, # minimum correlation coefficient for considering if adduct/iso
  add_ppm_tol = 10, # adduct mass tolerance (in ppm)
  adduct_file = adduct_list_filename,
  date_tag = "211202"
)

## Save output ----
write_csv(dat2, here(output_loc, "HILICparticulate_wide_adducts.csv"))
cat(
  "HILIC particulate:",
  round(sum(!is.na(dat2$add_annotation)) / length(dat2$add_annotation), digits = 2) * 100,
  "% of features marked as adduct or associated with an adduct.\n"
)

# RP, particulate -----
## Find adducts, clean up output -----
dat3 <- dereplicate_MFs(
  dat_filename1 = here(
    meta_data_loc,
    "particulate",
    "msdial",
    "RP_positive",
    "Area_1_202210241613.txt"
  ),
  dat_type_1 = "RPPos",
  dat_filename2 = NA,
  dat_type_2 = "",
  comment_tag = "m",
  sample_tag = "2112",
  drop_tags = c("DDA", "_Blk_"),
  RT_tol = 0.1, # minutes to allow for possible adduct/isotope
  cor_tol = 0.85, # minimum correlation coefficient for considering if adduct/iso
  add_ppm_tol = 15, # adduct mass tolerance (in ppm)
  adduct_file = adduct_list_filename,
  date_tag = "211207"
)

## Save output ----
write_csv(dat3, here(output_loc, "RPparticulate_wide_adducts.csv"))
cat(
  "PR particulate:",
  round(sum(!is.na(dat3$add_annotation)) / length(dat3$add_annotation), digits = 2) * 100,
  "% of features marked as adduct or associated with an adduct.\n"
)

# HILIC, dissolved ----
## Find adducts, clean up output -----
dat4 <- dereplicate_MFs(
  dat_filename1 = here(
    meta_data_loc,
    "dissolved",
    "msdial",
    "HILIC_positive",
    "Area_1_202210181739.txt"
  ),
  dat_type_1 = "HILICPos",
  dat_filename2 = here(
    meta_data_loc,
    "dissolved",
    "msdial",
    "HILIC_negative",
    "Area_1_20221019754.txt"
  ),
  dat_type_2 = "HILICNeg",
  comment_tag = "m",
  sample_tag = "220302",
  drop_tags = c("DDA", "Blk_"),
  RT_tol = 0.3, # minutes to allow for possible adduct/isotope
  cor_tol = 0.85, # minimum correlation coefficient for considering if adduct/iso
  add_ppm_tol = 10, # adduct mass tolerance (in ppm)
  adduct_file = adduct_list_filename,
  date_tag = "220302"
)

## Save output ----
write_csv(dat4, here(output_loc, "HILICdissolved_wide_adducts.csv"))
cat(
  "HILIC dissolved:",
  round(sum(!is.na(dat4$add_annotation)) / length(dat4$add_annotation), digits = 2) * 100,
  "% of features marked as adduct or associated with an adduct.\n"
)
