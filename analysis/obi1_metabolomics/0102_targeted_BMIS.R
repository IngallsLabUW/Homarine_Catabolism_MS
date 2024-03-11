# _____________________________________________________________________________
#       Copyright 2022, Integral Consulting Inc. All rights reserved.        #
# _____________________________________________________________________________
#
# _____________________________________________________________________________
#  0102_targeted_BMIS.R
# _____________________________________________________________________________
# PURPOSE:
#  Perform Best Matched Internal Standard Normalization
#
#
# PROJECT INFORMATION:
#   Name:
#   Number: C3352
#
# HISTORY:
# 	 Date		        Remarks
# _____________________________________________________________________________
# 	 20220420       Using this as a template: https://github.com/kheal/Gradients1_SemiTargeted3/blob/master/SourceCode/BMIS_QEHILIC.R and adapting for these data
#

# PACKAGES ----
library(tidyverse)
library(here)
source(here("01_targeted/targeted_dataprocessing_fxns.R"))

# BMIS on HILIC particulate ------
## Perform BMIS ----
BMIS_results <- BMIS(
  sample_key = "raw_input/sample_key.csv",
  is_names_file = "raw_input/is_names_particulate.csv",
  dat_pos_file = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Pos_Particulate.csv"
  ),
  dat_neg_file = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Neg_Particulate.csv"
  ),
  cut_off1 = 0.2,
  cut_off2 = 0.1,
  smps_to_dump = c(
    "211202_Poo_TruePooBacteria_DDApos20",
    "211202_Poo_TruePooBacteria_DDApos35",
    "211202_Poo_TruePooBacteria_DDApos50",
    "211202_Poo_TruePooBacteria_DDAneg20",
    "211202_Poo_TruePooBacteria_DDAneg35",
    "211202_Poo_TruePooBacteria_DDAneg50"
  ),
  is_to_dump = c("AMP, 15N5", "GMP, 15N5", "Succinic acid, 2H4", "Sulfoacetic acid, 13C2")
)
## Save output ----
### Internal standard reps
ggsave("01_targeted/data_intermediate/BMIS_checks/HILIC_Particulate_ISreps.pdf",
       BMIS_results[[1]],
       width = 15, height = 8, units = "in")

### Quick report
writeLines(BMIS_results[[2]], "01_targeted/data_intermediate/BMIS_checks/HILIC_Particulate_BMISSummary.txt")

### Internal standard reps
ggsave("01_targeted/data_intermediate/BMIS_checks/HILIC_Particulate_ISvIScheck.pdf",
       BMIS_results[[3]],
       width = 10, height = 8, units = "in")

### BMIS output
write_csv(BMIS_results[[4]],
          "01_targeted/data_intermediate/BMISd_HILIC_Particulate.csv")



# BMIS on HILIC dissolved ------
## Perform BMIS ----
BMIS_results <- BMIS(
  sample_key = "raw_input/sample_key.csv",
  is_names_file = "raw_input/is_names_dissolved.csv",
  dat_pos_file = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Pos_Dissolved.csv"
  ),
  dat_neg_file = here(
    "01_targeted",
    "data_intermediate",
    "QCd_HILIC_Neg_Dissolved.csv"
  ),
  cut_off1 = 0.2,
  cut_off2 = 0.1,
  smps_to_dump = c(),
  is_to_dump = c()
)
## Save output ----
### Internal standard reps
ggsave("01_targeted/data_intermediate/BMIS_checks/HILIC_Dissolved_ISreps.pdf",
       BMIS_results[[1]],
       width = 15, height = 8, units = "in")

### Quick report
writeLines(BMIS_results[[2]], "01_targeted/data_intermediate/BMIS_checks/HILIC_Dissolved_BMISSummary.txt")

### Internal standard reps
ggsave("01_targeted/data_intermediate/BMIS_checks/HILIC_Dissolved_ISvIScheck.pdf",
       BMIS_results[[3]],
       width = 10, height = 8, units = "in")

### BMIS output
write_csv(BMIS_results[[4]],
          "01_targeted/data_intermediate/BMISd_HILIC_Dissolved.csv")

