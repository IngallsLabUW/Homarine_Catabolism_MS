# Purpose: Script for populating data folders after project initialization

# History:
# Date	        Remarks
# 20240325      K Heal modifying to pull G4 Homarine experiment metabolomics data

library(tidyverse)
library(msconverteR)
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

# Convert the .raw files to mzML -----
get_pwiz_container()
raw_files <- list.files(export_folder, pattern = ".raw$", full.names = TRUE)

## convert to mzML (positive mode) -----
mzml_folder <- here("data", "raw", "metabolomics", "G4", "D3_Homarine_Fate_Inc","mzML", "positive")
dir.create(mzml_folder, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(raw_files)) {
    raw_file <- raw_files[i]
    output_file <- paste0(mzml_folder, "/", tools::file_path_sans_ext(basename(raw_file)), ".mzML")
    if (!file.exists(output_file)) {
        convert_files(raw_file, outpath =  mzml_folder, msconvert_args = "polarity positive" , docker_args = c())
    }
}

## convert to mzML (negative mode) -----
mzml_folder <- here("data", "raw", "metabolomics", "G4", "D3_Homarine_Fate_Inc","mzML", "negative")
dir.create(mzml_folder, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(raw_files)) {
    raw_file <- raw_files[i]
    output_file <- paste0(mzml_folder, "/", tools::file_path_sans_ext(basename(raw_file)), ".mzML")
    if (!file.exists(output_file)) {
        convert_files(raw_file, outpath =  mzml_folder, msconvert_args = "polarity positive" , docker_args = c())
    }
    convert_files(raw_file, outpath =  mzml_folder, msconvert_args = "polarity negative" , docker_args = c())
}