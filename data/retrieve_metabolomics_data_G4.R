# Purpose: Script for populating data folders after project initialization

# History:
# Date	        Remarks
# 20240325      K Heal modifying to pull G4 Homarine experiment metabolomics data
library(here)
library(tidyverse)
library(msconverteR)
source("R_functions/google_drive_functions.R")
mzml_convert = T


# Functions to help get the folder ID needed for the drive_download_daughters_function
# dir <- drive_find(
#     pattern = "G4_D3_Homarine_Fate_Incubations",
#     type = "folder",
#     shared_drive = "Ingalls Lab"
# )
# query = paste('"', dir$id, '"',  ' in parents', sep='')
# files = drive_ls(q=query, shared_drive = "Ingalls Lab", recursive = TRUE)

# Copy G4 Metabolomics data from Google Drive -----
## raw data -----
export_folder <- here("data", "raw", "metabolomics", "G4", "D3_Homarine_Fate_Inc","raw")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "1ED4NpEy1DhU875JyJq-tvgrdQeYw8CGy"
drive_download_daughters(gfolder, export_folder)

# Convert the .raw files to mzML -----
if (mzml_convert) {
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
}


# Copy .mzml centroided data and skyline files -----

## particulate HILIC Neg mzml and Skyline Files-----
export_folder <- here("data", "raw", "metabolomics", "g4", "particulate", "mzml_skyline", "negative")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "19yeehyg64AzgS2yCk8kU6_qs-i8DBEfQ"
drive_download_daughters(gfolder, export_folder)

## particulate HILIC Pos mzml and Skyline Files-----
export_folder <- here("data", "raw", "metabolomics", "g4", "particulate", "mzml_skyline", "positive")
dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
gfolder <- "19yJ641-QL2j4urP95Ay7KfK3wnpVxqAH"
drive_download_daughters(gfolder, export_folder)


# skyline-associated files ------
## associated  files-----
export_folder <- here("data", "raw", "metabolomics", "g4", "particulate", "mzml_skyline")
dir.create(export_folder, showWarnings = FALSE)
gfolder <- "19wxbPPxBsWX_hRy8sddUqY41qH7d3PoZ"
drive_download_daughters(gfolder, export_folder)

## internal standards and sample keys -----
## clean up the internal standards -----
# from copied data and clean names for future use
is_dissolved_dat <- read_csv(here("data", "raw", "metabolomics", "g4", "particulate", "mzml_skyline", "IS_names_dissolved.csv"),
                             show_col_types = FALSE
) %>%
    janitor::clean_names() %>%
    write_csv(here("data", "raw", "metabolomics", "g4", "is_names_dissolved.csv"))

is_particulate_dat <- read_csv(here("data", "raw", "metabolomics", "g4", "particulate", "mzml_skyline", "IS_names_particulate.csv"),
                               show_col_types = FALSE
) %>%
    janitor::clean_names() %>%
    write_csv(here("data", "raw", "metabolomics", "g4", "is_names_particulate.csv"))

## grab sample key -----
gfile <- drive_get(as_id("1IGbDGnFbXG4DclC9ELjqfftIYaTytZfz"), shared_drive = "Ingalls Lab")
drive_download(gfile,
               path = here("data", "raw", "metabolomics", "g4", "sample_key.csv")
)

