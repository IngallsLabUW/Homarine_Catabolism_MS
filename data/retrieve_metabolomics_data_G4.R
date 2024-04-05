# Purpose: Script for populating data folders after project initialization

# History:
# Date	        Remarks
# 20240325      K Heal modifying to pull G4 Homarine experiment metabolomics data

library(tidyverse)
source("R_functions/google_drive_functions.R")


# Functions to help get the folder ID needed for the drive_download_daughters_function
# dir <- drive_find(
#     pattern = "G4_D3_Homarine_Fate_Incubations",
#     type = "folder",
#     shared_drive = "Ingalls Lab"
# )
# query = paste('"', dir$id, '"',  ' in parents', sep='')
# files = drive_ls(q=query, shared_drive = "Ingalls Lab", recursive = TRUE)

# Copy RPom Metabolomics data from Google Drive -----
## raw data -----
export_folder <- here("data", "raw", "metabolomics", "G4", "D3_Homarine_Fate_Inc","raw")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "1ED4NpEy1DhU875JyJq-tvgrdQeYw8CGy"
drive_download_daughters(gfolder, export_folder)

#19yJ641-QL2j4urP95Ay7KfK3wnpVxqAH is positive, centroid
#19yeehyg64AzgS2yCk8kU6_qs-i8DBEfQ is negative, centroid