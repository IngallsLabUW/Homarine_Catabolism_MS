# This script extracts biosample metadata information from the NCBI datasets
# Input: JSONL files containing metadata information (folder of dataset directories)
# Output: CSV file containing the extracted biosample metadata information

# Load required libraries ----
library(here)
library(tidyverse)
library(jsonlite)
library(future.apply)

# Define function to process biosample data from JSONL files ----
get_biosample_df <- function(file) {
  attributes_oi <- c(
    "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "isolation_source", "lat_lon",
    "geographic location (latitude)", "geographic location (longitude)", "geographic location (region and locality)",
    "env_package"
  )

  if (!file.exists(file)) {
    message("File not found: ", file)
    return(NULL)
  }

  print(paste("Processing:", file))

  dat <- lapply(readLines(file), fromJSON)

  biosamples <- list()
  for (i in seq_along(dat)) {
    if (is.list(dat[[i]]$assemblyInfo$biosample)) {
      if (!is.null(dat[[i]]$assemblyInfo$biosample$attributes)) {
        biosamples[[i]] <- dat[[i]]$assemblyInfo$biosample$attributes
      }else{
          biosamples[[i]] <- list()
      }
        biosamples[[i]]$file <- file
        biosamples[[i]]$taxon_id <- dat[[i]]$taxId
        biosamples[[i]]$assembly_accession <- dat[[i]]$assemblyInfo$assemblyAccession
        biosamples[[i]]$biosample_id <- dat[[i]]$assemblyInfo$biosample$accession
      # Include bioproject information if available
        if (!is.null(dat[[i]]$assemblyInfo$bioprojectLineage$bioprojects)) {
            bioproject_list = dat[[i]]$assemblyInfo$bioprojectLineage$bioprojects
            
            # from each of the dataframes, pull out and join the accesion field
            bioproject_accessions = sapply(bioproject_list, function(x) x$accession)
            
            biosamples[[i]]$bioproject_accession <- 
                paste(bioproject_accessions, collapse = " | ")
            if ("PRJNA391943" %in% biosamples[[i]]$bioproject_accession){
            }
        }
      # Include project name if available
    if (!is.null(dat[[i]]$assemblyInfo$biosample$project_name)) {
        biosamples[[i]]$biosample_project_name <- paste(dat[[i]]$assemblyInfo$biosample$project_name, collapse = " | ")
      }
        
      # Include package and sample IDs if available
      if (!is.null(dat[[i]]$assemblyInfo$biosample$package)) {
        biosamples[[i]]$biosample_package <- dat[[i]]$assemblyInfo$biosample$package
      }
      if (!is.null(dat[[i]]$assemblyInfo$biosample$sampleIds)) {
        # Join all sample IDs into a single string
        biosamples[[i]]$alternative_identifiers <- paste(dat[[i]]$assemblyInfo$biosample$sampleIds$value, collapse = " | ")
      }
    }
  }

  if (length(biosamples) == 0) {
    message("No biosample data found in file: ", file)
    return(NULL)
  }
   # browser()
  biosample_df <- bind_rows(biosamples)

  # Filter only for the attributes of interest
  biosample_filtered <- biosample_df %>%
    filter(name %in% attributes_oi) %>%
    pivot_wider(names_from = name, values_from = value, values_fn = list)

  # Ensure all attributes are of uniform length (convert to single-value per cell)
  cols <- colnames(biosample_filtered)
  cols_to_unnest <- cols[cols %in% attributes_oi]

  biosample_filtered <- biosample_filtered %>%
    mutate(
      across(
        all_of(cols_to_unnest),
        ~ map_chr(
          .,
          ~ ifelse(length(.) == 1, ., paste(., collapse = "; ")),
          .default = NA_character_
        )
      )
    ) %>%
    unnest(cols = all_of(cols_to_unnest))

  return(biosample_filtered)
}

# Locate the dataset directories containing the JSONL files----
dat_loc <- normalizePath("data/raw/genomes/ncbi_metadata/downloaded_taxons/")
dataset_dirs <- list.dirs(dat_loc, recursive = FALSE)
dataset_dirs <- dataset_dirs[grepl("dataset_\\d+$", dataset_dirs)] # Keep only dataset directories

files <- file.path(dataset_dirs, "ncbi_dataset", "data", "assembly_data_report.jsonl")
files <- files[file.exists(files)]
file <- "/Users/heal742/LOCAL/07_UW/Homarine_Catabolism_MS/data/raw/genomes/ncbi_metadata/downloaded_taxons/dataset_1486246/ncbi_dataset/data/assembly_data_report.jsonl"


# Process the JSONL files using parallel processing ----
plan(multisession)
biosample_dfs <- future_lapply(files, get_biosample_df)

# Combine results into a single data frame ----
biosample_data <- bind_rows(biosample_dfs, .id = "dataset")

# Save the output ----
write_csv(biosample_data, here("data", "intermediate", "ncbi_metadata_biosample_info.csv"))

# Show completion message
print(paste("Processing complete. Output saved to:", here("data", "intermediate", "ncbi_metadata_biosample_info.csv")))
