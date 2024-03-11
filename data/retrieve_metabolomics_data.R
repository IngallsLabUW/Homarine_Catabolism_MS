# Purpose: Script for populating data folders after project initialization

# History:
# Date	     Remarks
# 20240311   K Heal starting script

library(tidyverse)
source("R_functions/google_drive_functions.R")



# Functions to help get the folder ID needed for the drive_download_daughters_function
# dir = drive_find(
#   pattern='OscarSosa_Bacteria',
#   type='folder',
#   shared_drive = "Ingalls Lab")
# query = paste('"', dir$id, '"',  ' in parents', sep='')
# files = drive_ls(q=query, shared_drive = "Ingalls Lab", recursive = TRUE)


# Copy OBi1 metabolomics data -----
## mzml data -----
### dissolved HILIC Neg mzml-----
export_folder <- here("data","raw", "metabolomics", "obi1", "dissolved", "HILIC_negative_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "17c1AaWndGr4TH_aZZH10_N4ytxnBp9A1"
drive_download_daughters(gfolder, export_folder)

### dissolved HILIC Pos mzml-----
export_folder <-  here("data","raw", "metabolomics", "obi1", "dissolved", "HILIC_positive_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "17btbYSFASS6FuCfTp_3B9c12Jv2SUnka"
drive_download_daughters(gfolder, export_folder)

### particulate HILIC Neg mzml-----
export_folder <-here("data","raw", "metabolomics", "obi1", "particulate", "HILIC_negative_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "11mQ30CtXReGwuTZdyq0NMEb8kX6471G4"
drive_download_daughters(gfolder, export_folder)

### particulate HILIC Pos mzml-----
export_folder <- here("data","raw", "metabolomics", "obi1", "particulate", "HILIC_positive_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "11gy9fnvy9l5gQ0RY9JUlKEXyZA-ml3yc"
drive_download_daughters(gfolder, export_folder)

### particulate RP mzml-----
export_folder <- here("data","raw", "metabolomics", "obi1", "particulate", "RP_positive_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "1-DUbiFRF7rD7XVO2ST6duSO7NuLMGkdR"
drive_download_daughters(gfolder, export_folder, refresh_mzml)


## skyline and associated files ------
### dissolved skyline and associated  files-----
export_folder <- here("data","raw", "metabolomics", "obi1",  "dissolved", "skyline")
dir.create(export_folder, showWarnings = FALSE)
gfolder <- "12iWeaN5D3m7aUt6Q9sRbBxzTELbSQH_Z"
drive_download_daughters(gfolder, export_folder)

### particulate skyline and associated files-----
# HILIC for all samples and compounds; RP for just internal standards
export_folder <- here("data","raw", "metabolomics", "obi1",  "particulate", "skyline")
dir.create(export_folder, showWarnings = FALSE)
gfolder <- "1-EO6dbIP3gb661-HhKfdHxber4IwqwEf"
drive_download_daughters(gfolder, export_folder)