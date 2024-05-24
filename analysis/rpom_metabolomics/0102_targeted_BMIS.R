# Purpose: BMIS on all targeted datasets for the Rpom Metabolomics

# History:
# Date	     Remarks
# 20240415   K Heal transferring script from Obi1 to Rpom specific repo to general repo


# PACKAGES ----
library(tidyverse)
library(here)
source(here("analysis", "obi1_metabolomics", "functions", "targeted_dataprocessing_fxns.R"))

# SET FILE LOCS ----
meta_data_loc <- here("data", "raw", "metabolomics", "rpom")
qc_data_loc <- here("data", "intermediate", "metabolomics", "rpom", "targeted")
bmis_checks_loc <- here(qc_data_loc, "BMIS_checks")

# MAKE DIR FOR BMIS CHECKS
dir.create(bmis_checks_loc, recursive = TRUE, showWarnings = FALSE)

# GET LISTS OF DISSOLVED/POSITIVE REPLICATES
#TODO KRH: Rerun with updated sample key (need inject volume on all samples)
sample_key <- read_csv(here(meta_data_loc, "sample_key.csv"), show_col_types = FALSE)
samps_dissolved <- sample_key %>%
  filter(str_detect(sample_fraction, "dissolved")) %>%
  pull(replicate_name)
samps_particulate <- sample_key %>%
  filter(str_detect(sample_fraction, "particulate")) %>%
  pull(replicate_name)

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
    "QCd_HILIC_Pos.csv"
  ),
  dat_neg_file = here(
    qc_data_loc,
    "QCd_HILIC_Neg.csv"
  ),
  cut_off1 = 0.2,
  cut_off2 = 0.1,
  smps_to_dump = samps_dissolved,
  is_to_dump = c()
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
    "QCd_HILIC_Pos.csv"
  ),
  dat_neg_file = here(
    qc_data_loc,
    "QCd_HILIC_Neg.csv"
  ),
  cut_off1 = 0.2,
  cut_off2 = 0.1,
  smps_to_dump = samps_particulate,
  is_to_dump = c("L-Methionine, 2H3", "Homarine, 2H3")
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
