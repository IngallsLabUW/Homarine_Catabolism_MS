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
        "geographic location (latitude)", "geographic location (longitude)", "geographic location (region and locality)"
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
                biosamples[[i]]$taxon_id <- dat[[i]]$taxId
                biosamples[[i]]$assembly_accession <- dat[[i]]$assemblyInfo$assemblyAccession
            }
        }
    }
    
    if (length(biosamples) == 0) {
        message("No biosample data found in file: ", file)
        return(NULL)
    }
    
    biosample_df <- bind_rows(biosamples)
    
    # Filter only for the attributes of interest
    biosample_filtered <- biosample_df %>%
        filter(name %in% attributes_oi) %>%
        pivot_wider(names_from = name, values_from = value, values_fn = list)
    
    # Ensure all attributes are of uniform length (convert to single-value per cell)
    cols <- colnames(biosample_filtered)
    cols_to_unnest <- cols[cols %in% attributes_oi]
    
    biosample_filtered <- biosample_filtered %>%
        mutate(across(all_of(cols_to_unnest), ~ map_chr(., ~ ifelse(length(.) == 1, ., paste(., collapse = "; ")), .default = NA_character_))) %>%
        unnest(cols = all_of(cols_to_unnest))
    
    return(biosample_filtered)
}

# Locate the dataset directories containing the JSONL files----
dat_loc <- normalizePath("C:/Bio/NCBI_datasets/data/taxon_downloads")
dataset_dirs <- list.dirs(dat_loc, recursive = FALSE)
dataset_dirs <- dataset_dirs[grepl("dataset_\\d+$", dataset_dirs)]  # Keep only dataset directories
files <- file.path(dataset_dirs, "ncbi_dataset", "data", "assembly_data_report.jsonl")
files <- files[file.exists(files)]

# Process the JSONL files using parallel processing ----
plan(multisession)
biosample_dfs <- future_lapply(files, get_biosample_df)

# Combine results into a single data frame ----
biosample_data <- bind_rows(biosample_dfs, .id = "dataset")

# Save the output ----
#output_file <- file.path(dat_loc, "ncbi_metadata_biosample_info.csv")
#write_csv(biosample_data, output_file)
write_csv(final_df, here("data", "intermediate", "ncbi_metadata_biosample_info.csv"))

# Show completion message
print(paste("Processing complete. Output saved to:", output_file))
