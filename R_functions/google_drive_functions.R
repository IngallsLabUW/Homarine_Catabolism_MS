# Purpose: Collection of functions related to reading/writing data from google drive

# History:
# Date	     Remarks
# 20230403    K Heal starting script


# PACKAGES & SPECIAL FUNCTIONS ----
library(googledrive)
library(here)
library(dplyr)


#' Download Daughters from a Google Drive Folder
#'
#' @param gfolder_id google drive folder ID
#' @param folder_export folder to export google drive daughters to
#' @param overwrite T/F if would like to overwrite the data, defaults to False
#'
#' @return NULL
#' @export

drive_download_daughters <- function(gfolder_id, folder_export, overwrite = F) {
  folder_id <- drive_get(as_id(gfolder_id), shared_drive = "Ingalls Lab")
  files <- drive_ls(folder_id)
  glimpse(files)
  for (i in seq_along(files$name)) {
    print(paste(i, "of", length(files$name)))
    try(drive_download(files[i, ],
                       path = here(folder_export, files$name[i]),
                       overwrite = overwrite
    ))
  }
}
