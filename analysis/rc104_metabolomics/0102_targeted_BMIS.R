#TODO: Change for RC104 anaylsis
# Purpose: BMIS on all targeted datasets for the rc104 Metabolomics

# History:
# Date	     Remarks
# 20250106   K Heal transferring script from G5 to rc104


# PACKAGES ----
library(tidyverse)
library(here)
source(here("analysis", "obi1_metabolomics", "functions", "targeted_dataprocessing_fxns.R"))

# SET FILE LOCS ----

meta_data_loc <- here("data", "raw", "metabolomics", "RC104")
qc_data_loc <- here("data", "intermediate", "metabolomics", "rc104", "targeted")
bmis_checks_loc <- here(qc_data_loc, "BMIS_checks")

# MAKE DIR FOR BMIS CHECKS
dir.create(bmis_checks_loc, recursive = TRUE, showWarnings = FALSE)

# GET LISTS OF DISSOLVED/POSITIVE REPLICATES
sample_key <- read_csv(here(meta_data_loc, "sample_key.csv"), show_col_types = FALSE)
samps_dissolved <- sample_key %>%
  filter(str_detect(sample_fraction, "dissolved")) %>%
  pull(replicate_name)
samps_particulate <- sample_key %>%
  filter(str_detect(sample_fraction, "particulate")) %>%
  pull(replicate_name)

# BMIS on HILIC Positive particulate ------
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
  cut_off1 = 0.2,
  cut_off2 = 0.1,
  smps_to_dump = samps_dissolved,
  is_to_dump = c()
)


## Save output ----
### Internal standard reps -----
ggsave(
  here(
    bmis_checks_loc, "HILIC_Pos_Particulate_ISreps.pdf"
  ),
  BMIS_results[[1]],
  width = 15, height = 8, units = "in"
)

### Quick report----
writeLines(
  BMIS_results[[2]],
  here(
    bmis_checks_loc, "HILIC_Pos_Particulate_BMISSummary.txt"
  )
)

### Internal standard reps----
ggsave(
  here(
    bmis_checks_loc, "HILIC_Pos_Particulate_ISvIScheck.pdf"
  ),
  BMIS_results[[3]],
  width = 10, height = 8, units = "in"
)

### BMIS output----
write_csv(
  BMIS_results[[4]],
  here(qc_data_loc, "BMISd_Pos_HILIC_Particulate.csv")
)
print("BMIS complete on HILIC particulate positive")


# BMIS on HILIC Negative particulate ------
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
        bmis_checks_loc, "HILIC_Neg_Particulate_ISreps.pdf"
    ),
    BMIS_results[[1]],
    width = 15, height = 8, units = "in"
)

### Quick report----
writeLines(
    BMIS_results[[2]],
    here(
        bmis_checks_loc, "HILIC_Neg_Particulate_BMISSummary.txt"
    )
)

### Internal standard reps----
ggsave(
    here(
        bmis_checks_loc, "HILIC_Neg_Particulate_ISvIScheck.pdf"
    ),
    BMIS_results[[3]],
    width = 10, height = 8, units = "in"
)

### BMIS output----
write_csv(
    BMIS_results[[4]],
    here(qc_data_loc, "BMISd_Neg_HILIC_Particulate.csv")
)
print("BMIS complete on HILIC particulate negative")
