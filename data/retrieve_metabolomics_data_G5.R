# Purpose: Script for populating data folders after project initialization

# History:
# Date	        Remarks
# 20240325      K Heal modifying to pull G5 Homarine experiment metabolomics data

library(tidyverse)
source("R_functions/google_drive_functions.R")


# Functions to help get the folder ID needed for the drive_download_daughters_function
dir <- drive_find(
    pattern = "Labeled_Homarine_Incubations",
    type = "folder",
    shared_drive = "Ingalls Lab"
)
query = paste('"', dir$id, '"',  ' in parents', sep='')
files = drive_ls(q=query, shared_drive = "Ingalls Lab", recursive = TRUE)

# Copy RPom Metabolomics data from Google Drive -----
## mzml data -----
### dissolved HILIC Neg mzml-----
export_folder <- here("data", "raw", "metabolomics", "G5", "particulate","mzml", "HILIC_negative")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "18mIQvL16R-KuWay6-2QogVW03ZCab5TK"
drive_download_daughters(gfolder, export_folder)

### dissolved HILIC Pos mzml-----
export_folder <- here("data", "raw", "metabolomics", "G5", "particulate","mzml", "HILIC_positive")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "12gwnpFB8uqzw7AVs11pPkikv3vkbK0lf"
drive_download_daughters(gfolder, export_folder)

## grab sample key -----
gfile <- drive_get(as_id("15AMo5diyaLdvJIeagqATJQOYt_6WVyY3"), shared_drive = "Ingalls Lab")
drive_download(gfile,
               path = here("data", "raw", "metabolomics", "G5", "sample_key.csv"))