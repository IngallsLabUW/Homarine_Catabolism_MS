# Purpose: Perform BMIS on untargeted data
# History:
# 20240313      KRH  moving from Obi1 repo

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
source(here("analysis", "obi1_metabolomics", "functions", "untargeted_dataprocessing_fxns.R"))

# SET FILE LOCS ----
meta_data_loc <- here("data", "raw", "metabolomics", "obi1")
adduct_list_filename <- here("data", "raw", "metabolomics", "adduct_list.csv")
output_loc <- here("data", "intermediate", "metabolomics", "obi1", "untargeted")
bmis_checks_loc <- here(output_loc, "BMIS_checks")

# MAKE DIR FOR BMIS CHECKS
dir.create(bmis_checks_loc, recursive = TRUE, showWarnings = FALSE)

# BMIS on HILIC particulate ------
## Perform BMIS ----
BMIS_results_hp <- BMIS_MSdialoutput(
  dat1_filename = here(
    output_loc,
    "HILICparticulate_wide_adducts.csv"
  ),
  sample_key = here(
    meta_data_loc,
    "sample_key.csv"
  ),
  is_names_file = here(
    meta_data_loc,
    "is_names_particulate.csv"
  ),
  is_dat1_filename = here(
    meta_data_loc,
    "particulate",
    "skyline",
    "HILIC_Pos_Particulate_OscarSosaBacteria.csv"
  ),
  is_dat2_filename = here(
    meta_data_loc,
    "particulate",
    "skyline",
    "HILIC_Neg_Particulate_OscarSosaBacteria.csv"
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
  is_to_dump = c("AMP, 15N5", "GMP, 15N5", "Succinic acid, 2H4")
)

## Save output ----
### Internal standard reps ----
ggsave(
  here(
    bmis_checks_loc,
    "HILIC_Particulate_ISreps.pdf"
  ),
  BMIS_results_hp[[1]],
  width = 15, height = 8, units = "in"
)

### Quick report ----
writeLines(BMIS_results_hp[[2]], here(
  bmis_checks_loc,
  "HILIC_Particulate_BMISSummary.txt"
))

### Internal standard reps ----
ggsave(here(bmis_checks_loc, "HILIC_Particulate_ISvIScheck.pdf"),
  BMIS_results_hp[[3]],
  width = 10, height = 8, units = "in"
)

### BMIS output -----
write_csv(
  BMIS_results_hp[[4]],
  here(output_loc, "BMISd_HILIC_Particulate.csv")
)



# BMIS on HILIC dissolved ------
## Perform BMIS ----
BMIS_results_hd <- BMIS_MSdialoutput(
  dat1_filename = here(
    output_loc,
    "HILICdissolved_wide_adducts.csv"
  ),
  sample_key = here(
    meta_data_loc,
    "sample_key.csv"
  ),
  is_names_file = here(
    meta_data_loc,
    "is_names_dissolved.csv"
  ),
  is_dat1_filename = here(
    meta_data_loc,
    "dissolved",
    "skyline",
    "HILIC_Pos_Dissolved_OscarSosaBacteria_zerod.csv"
  ),
  is_dat2_filename = here(
    meta_data_loc,
    "dissolved",
    "skyline",
    "HILIC_Neg_Dissolved_OscarSosaBacteria_zerod.csv"
  ),
  cut_off1 = 0.2,
  cut_off2 = 0.1,
  smps_to_dump = c(),
  is_to_dump = c()
)

## Save output ----
ggsave(
  here(
    bmis_checks_loc,
    "HILIC_Dissolved_ISreps.pdf"
  ),
  BMIS_results_hd[[1]],
  width = 15, height = 8, units = "in"
)

### Quick report ----
writeLines(BMIS_results_hd[[2]], here(
  bmis_checks_loc,
  "HILIC_Dissolved_BMISSummary.txt"
))

### Internal standard reps ----
ggsave(here(bmis_checks_loc, "HILIC_Dissolved_ISvIScheck.pdf"),
  BMIS_results_hd[[3]],
  width = 10, height = 8, units = "in"
)

### BMIS output -----
write_csv(
  BMIS_results_hd[[4]],
  here(output_loc, "BMISd_HILIC_Dissolved.csv")
)


# BMIS on RP particulate ------
## Perform BMIS ----
# Waiting on LTC to give us updated is_names
# BMIS_results_rp <- BMIS_MSdialoutput(
#   dat1_filename = "02_untargeted/intermediates/RPparticulate_wide_adducts.csv",
#   sample_key = "raw_input/sample_key.csv",
#   is_names_filename = "raw_input/is_names_particulate.csv",
#   is_dat1_filename = "raw_input/particulate/skyline/RP_OscarSosa_InternalStandards.csv",
#   is_dat2_filename = NA,
#   cut_off1 = 0.2,
#   cut_off2 = 0.1,
#   smps_to_dump = c(),
#   is_to_dump = c()
# )
#
# ## Save output ----
# ### Internal standard reps ----
# ggsave("02_untargeted/intermediates/BMIS_checks/HILIC_Particulate_ISreps.pdf",
#        BMIS_results_hp[[1]],
#        width = 15, height = 8, units = "in")
#
# ### Quick report ----
# writeLines(BMIS_results_hp[[2]], "02_untargeted/intermediates/BMIS_checks/HILIC_Particulate_BMISSummary.txt")
#
# ### Internal standard reps ----
# ggsave("02_untargeted/intermediates/BMIS_checks/HILIC_Particulate_ISvIScheck.pdf",
#        BMIS_results_hp[[3]],
#        width = 10, height = 8, units = "in")
#
# ### BMIS output -----
# write_csv(BMIS_results_hp[[4]],
#           "02_untargeted/intermediates/BMISd_HILIC_Particulate.csv")
