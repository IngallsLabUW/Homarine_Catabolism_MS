# Purpose: BMIS on all targeted datasets for the Obi1 Metabolomics

# History:
# Date	     Remarks
# 20240311   K Heal transferring script from Obi1 specific repo to general repo
# 20220420       Using this as a template: https://github.com/kheal/Gradients1_SemiTargeted3/blob/master/SourceCode/BMIS_QEHILIC.R and adapting for these data


# PACKAGES ----
library(tidyverse)
library(here)
source(here("analysis", "obi1_metabolomics", "functions", "targeted_dataprocessing_fxns.R"))

# SET FILE LOCS ----
meta_data_loc <- here("data", "raw", "metabolomics", "obi1")
qc_data_loc <- here("data", "intermediate", "metabolomics", "obi1", "targeted")
bmis_checks_loc <- here(qc_data_loc, "BMIS_checks")

# MAKE DIR FOR BMIS CHECKS
dir.create(bmis_checks_loc, recursive = TRUE, showWarnings = FALSE)

# BMIS on HILIC particulate ------
## Perform BMIS ----
BMIS_results <- BMIS(
  sample_key = here(
    meta_data_loc,
    "sample_key.csv"
  ),
  is_names_file = here(
    meta_data_loc,
    "is_names_particulate.csv"
  ),
  dat_pos_file = here(
    qc_data_loc,
    "QCd_HILIC_Pos_Particulate_combined.csv"
  ),
  dat_neg_file = here(
    qc_data_loc,
    "QCd_HILIC_Neg_Particulate_combined.csv"
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
  is_to_dump = c("AMP, 15N5", "GMP, 15N5", "Succinic acid, 2H4", "Sulfoacetic acid, 13C2"),
  skip_lines = 0
)
## Save output ----
### Internal standard reps -----
ggsave(
  here(
    bmis_checks_loc, "HILIC_Particulate_ISreps.pdf"
  ),
  BMIS_results[[1]],
  width = 15, height = 8, units = "in"
)

### Quick report----
writeLines(
  BMIS_results[[2]],
  here(
    bmis_checks_loc, "HILIC_Particulate_BMISSummary.txt"
  )
)

### Internal standard reps----
ggsave(
  here(
    bmis_checks_loc, "HILIC_Particulate_ISvIScheck.pdf"
  ),
  BMIS_results[[3]],
  width = 10, height = 8, units = "in"
)

### BMIS output----
write_csv(
  BMIS_results[[4]],
  here(qc_data_loc, "BMISd_HILIC_Particulate.csv")
)
print("BMIS complete on HILIC particulate")


# BMIS on HILIC dissolved ------
## Perform BMIS ----
BMIS_results <- BMIS(
  sample_key = here(
    meta_data_loc,
    "sample_key.csv"
  ),
  is_names_file = here(
    meta_data_loc,
    "is_names_dissolved.csv"
  ),
  dat_pos_file = here(
    qc_data_loc,
    "QCd_HILIC_Pos_Dissolved_combined.csv"
  ),
  dat_neg_file = here(
    qc_data_loc,
    "QCd_HILIC_Neg_Dissolved_combined.csv"
  ),
  cut_off1 = 0.2,
  cut_off2 = 0.1,
  smps_to_dump = c(),
  is_to_dump = c(),
  skip_lines = 0
)
## Save output ----
## Save internal standard reps ----
ggsave(
  here(
    bmis_checks_loc,
    "HILIC_Dissolved_ISreps.pdf"
  ),
  BMIS_results[[1]],
  width = 15, height = 8, units = "in"
)

### Quick report-----
writeLines(
  BMIS_results[[2]],
  here(
    bmis_checks_loc,
    "HILIC_Dissolved_BMISSummary.txt"
  )
)


### Internal standard reps-----
ggsave(
  here(
    bmis_checks_loc,
    "HILIC_Dissolved_ISvIScheck.pdf"
  ),
  BMIS_results[[3]],
  width = 10, height = 8, units = "in"
)

### BMIS output-----
write_csv(
  BMIS_results[[4]],
  here(qc_data_loc, "BMISd_HILIC_Dissolved.csv")
)
print("BMIS complete on HILIC dissolved")
