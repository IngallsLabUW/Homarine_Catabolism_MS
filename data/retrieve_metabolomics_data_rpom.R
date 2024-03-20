# Purpose: Script for populating data folders after project initialization

# History:
# Date	     Remarks
# 20240320   Downloading Rpom metabolomics data

library(tidyverse)
source("R_functions/google_drive_functions.R")


# Functions to help get the folder ID needed for the drive_download_daughters_function
dir <- drive_find(
    pattern = "EPIC",
    type = "folder",
    shared_drive = "Ingalls Lab"
)
# query = paste('"', dir$id, '"',  ' in parents', sep='')
# files = drive_ls(q=query, shared_drive = "Ingalls Lab", recursive = TRUE)

#from original script. no dissolved data for rpom yet
# # Copy RPom Metabolomics data from Google Drive -----
# ## mzml data -----
# ### dissolved HILIC Neg mzml-----
# export_folder <- here("data", "raw", "metabolomics", "rpom", "dissolved", "mzml", "HILIC_negative_mzml")
# dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
# gfolder <- "17c1AaWndGr4TH_aZZH10_N4ytxnBp9A1" #change this
# drive_download_daughters(gfolder, export_folder)
#
# ### dissolved HILIC Pos mzml-----
# export_folder <- here("data", "raw", "metabolomics", "rpom", "dissolved", "mzml", "HILIC_positive_mzml")
# dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
# gfolder <- "17btbYSFASS6FuCfTp_3B9c12Jv2SUnka" #change this
# drive_download_daughters(gfolder, export_folder)

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

# commented out below is code from the original script, so it gives obi1 data

# ### particulate RP mzml-----
# export_folder <- here("data", "raw", "metabolomics", "obi1", "particulate", "mzml", "RP_positive_mzml")
# dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
# gfolder <- "1-DUbiFRF7rD7XVO2ST6duSO7NuLMGkdR"
# drive_download_daughters(gfolder, export_folder, refresh_mzml)

# ## skyline and associated files ------
# ### dissolved skyline and associated  files-----
# export_folder <- here("data", "raw", "metabolomics", "obi1", "dissolved", "skyline")
# dir.create(export_folder, showWarnings = FALSE)
# gfolder <- "12iWeaN5D3m7aUt6Q9sRbBxzTELbSQH_Z"
# drive_download_daughters(gfolder, export_folder)
#
# ### particulate skyline and associated files-----
# # HILIC for all samples and compounds; RP for just internal standards
# export_folder <- here("data", "raw", "metabolomics", "obi1", "particulate", "skyline")
# dir.create(export_folder, showWarnings = FALSE)
# gfolder <- "1-EO6dbIP3gb661-HhKfdHxber4IwqwEf"
# drive_download_daughters(gfolder, export_folder)

# ## msdial and associated files ------
# ### dissolved HILIC Neg -----
# export_folder <- here("data", "raw", "metabolomics", "obi1", "dissolved", "msdial", "HILIC_negative")
# dir.create(export_folder, showWarnings = FALSE)
# gfolder <- "1ZNoFGayGlipvys00VtF9uE74bXH4SHyJ"
# drive_download_daughters(gfolder, export_folder)
#
# ### dissolved HILIC Pos -----
# export_folder <- here("data", "raw", "metabolomics", "obi1", "dissolved", "msdial", "HILIC_positive")
# dir.create(export_folder, showWarnings = FALSE, recursive = TRUE)
# gfolder <- "1O4VYXDmKtVpy90qeTktv8wV1DO2WeGdG"
# drive_download_daughters(gfolder, export_folder)
#
# ### particulate HILIC Neg -----
# export_folder <- here("data", "raw", "metabolomics", "obi1", "particulate", "msdial", "HILIC_negative")
# dir.create(export_folder, showWarnings = FALSE, recursive = TRUE)
# gfolder <- "1rg16L01jW2FHSGzd6_PHNB_kzgpFZ58J"
# drive_download_daughters(gfolder, export_folder)
# ### particulate HILIC Pos -----
# export_folder <- here("data", "raw", "metabolomics", "obi1", "particulate", "msdial", "HILIC_positive")
# dir.create(export_folder, showWarnings = FALSE, recursive = TRUE)
# gfolder <- "1rKMrjGk_5KBoip1hDViwcKPJ5vaMbSDt"
# drive_download_daughters(gfolder, export_folder)
# ### particulate RP -----
# export_folder <- here("data", "raw", "metabolomics", "obi1", "particulate", "msdial", "RP_positive")
# dir.create(export_folder, showWarnings = FALSE, recursive = TRUE)
# gfolder <- "1nJzw3bUZ1XNuKvNK89wyKPPn-ezzpRfA"
# drive_download_daughters(gfolder, export_folder)

# # internal standards and sample keys -----
# ## clean up the internal standards -----
# # from copied data and clean names for future use
# is_dissolved_dat <- read_csv(here("data", "raw", "metabolomics", "obi1", "dissolved", "skyline", "IS_names_dissolve.csv"),
#                              show_col_types = FALSE
# ) %>%
#     janitor::clean_names() %>%
#     write_csv(here("data", "raw", "metabolomics", "obi1", "is_names_dissolved.csv"))
#
# is_particulate_dat <- read_csv(here("data", "raw", "metabolomics", "obi1", "particulate", "skyline", "IS_names_particulate.csv"),
#                                show_col_types = FALSE
# ) %>%
#     janitor::clean_names() %>%
#     write_csv(here("data", "raw", "metabolomics", "obi1", "is_names_particulate.csv"))
#
# ## grab sample key -----
# gfile <- drive_get(as_id("1KL7-kJOZrlDtV5wToq_ejkLNtBiGNC5Y"), shared_drive = "Ingalls Lab")
# drive_download(gfile,
#                path = here("data", "raw", "metabolomics", "obi1", "sample_key.csv")
# )
