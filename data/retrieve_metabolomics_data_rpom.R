# Purpose: Script for populating data folders after project initialization

# History:
# Date	     Remarks
# 20240320   Downloading Rpom metabolomics data

library(tidyverse)
source("R_functions/google_drive_functions.R")


# Functions to help get the folder ID needed for the drive_download_daughters_function
# dir <- drive_find(
#    pattern = "EPIC",
#    type = "folder",
#    shared_drive = "Ingalls Lab"
# )
# query = paste('"', dir$id, '"',  ' in parents', sep='')
# files = drive_ls(q=query, shared_drive = "Ingalls Lab", recursive = TRUE)
# Folder location: https://drive.google.com/drive/u/1/folders/1engPa64gpbKFOhTYUNpTwPTZlknVFbZP

# Copy RPom Metabolomics data from Google Drive -----
## mzml data -----
### dissolved HILIC Neg mzml-----
export_folder <- here("data", "raw", "metabolomics", "rpom", "dissolved", "mzml", "HILIC_negative_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "1-2vr0a4QUVq72FjZCtL17TmPhIcqaX9W" # change this
drive_download_daughters(gfolder, export_folder)

### dissolved HILIC Pos mzml-----
export_folder <- here("data", "raw", "metabolomics", "rpom", "dissolved", "mzml", "HILIC_positive_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "12091AHTwpfcNKBcEm0RG_keKf-gjzcYr" # change this
drive_download_daughters(gfolder, export_folder)

### particulate HILIC Neg mzml-----
export_folder <- here("data", "raw", "metabolomics", "rpom", "particulate", "mzml", "HILIC_negative_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "1g4uLeHkJ5ve7dy7CoElWN6PnpskQywof"
drive_download_daughters(gfolder, export_folder)

### particulate HILIC Pos mzml-----
export_folder <- here("data", "raw", "metabolomics", "rpom", "particulate", "mzml", "HILIC_positive_mzml")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "1dENq4VdLRXw-FNZO0bI1M-xScOTLpKeS"
drive_download_daughters(gfolder, export_folder)


# skyline and associated files ------
## particulate and dissolved skyline and associated  files-----
export_folder <- here("data", "raw", "metabolomics", "rpom", "skyline")
dir.create(export_folder, showWarnings = FALSE)
gfolder <- "1engPa64gpbKFOhTYUNpTwPTZlknVFbZP"
drive_download_daughters(gfolder, export_folder)

## internal standards and sample keys -----
## clean up the internal standards -----
# from copied data and clean names for future use
is_dissolved_dat <- read_csv(here("data", "raw", "metabolomics", "rpom", "skyline", "IS_names_dissolved.csv"),
  show_col_types = FALSE
) %>%
  janitor::clean_names() %>%
  write_csv(here("data", "raw", "metabolomics", "rpom", "is_names_dissolved.csv"))

is_particulate_dat <- read_csv(here("data", "raw", "metabolomics", "rpom", "skyline", "IS_names_particulate.csv"),
  show_col_types = FALSE
) %>%
  janitor::clean_names() %>%
  write_csv(here("data", "raw", "metabolomics", "rpom", "is_names_particulate.csv"))

## grab sample key -----
gfile <- drive_get(as_id("11KjeNysYSckl80N_M5dGxS5VFkhYbGtMmEAvvF69XAQ"), shared_drive = "Ingalls Lab")
drive_download(gfile,
  path = here("data", "raw", "metabolomics", "rpom", "sample_key.csv")
)

