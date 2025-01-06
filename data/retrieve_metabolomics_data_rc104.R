# Purpose: Script for populating data folders after project initialization

# History:
# Date	        Remarks
# 20250106      K Heal modifying to pull RC104 Homarine fate experiment metabolomics data

library(tidyverse)
source("R_functions/google_drive_functions.R")


# Functions to help get the folder ID needed for the drive_download_daughters_function
#dir <- drive_find(
#    pattern = "RC104_HomarineFate",
#    type = "folder",
#    shared_drive = "Ingalls Lab"
#)
# query = paste('"', dir$id, '"',  ' in parents', sep='')
query = "\"1-xt3TXKEtClt5WUEIHc5nMgfhNoompuo\" in parents"
files = drive_ls(q=query, shared_drive = "Ingalls Lab", recursive = TRUE)

# Copy RPom Metabolomics data from Google Drive -----
## mzml data -----
### dissolved HILIC Neg mzml-----
export_folder <- here("data", "raw", "metabolomics", "RC104", "particulate","mzml", "HILIC_negative")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "1JDoNK2W2ps3twLP0VuzaCev_wtlP2cjt"
drive_download_daughters(gfolder, export_folder)

### dissolved HILIC Pos mzml-----
export_folder <- here("data", "raw", "metabolomics", "RC104", "particulate","mzml", "HILIC_positive")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "16oCaT5FE4dnn6FY9qqbJhBQ6TF5hollC"
drive_download_daughters(gfolder, export_folder)

## grab sample key -----
gfile <- drive_get(as_id("10kEQgpEJONqSQqg0m7OCwnSn-SUsxtDR"), shared_drive = "Ingalls Lab")
drive_download(gfile,
               path = here("data", "raw", "metabolomics", "RC104", "sample_key.csv"))

## grab internal standard names file ----
gfile <- drive_get(as_id("1b-DmN87sWI5Lf8xJgNrWSoOInZ-rNGln"), shared_drive = "Ingalls Lab")
drive_download(gfile,
               path = here("data", "raw", "metabolomics", "RC104", "is_names_particulate.csv"))
