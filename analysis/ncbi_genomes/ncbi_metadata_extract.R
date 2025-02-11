# Get geo location information from the ncbi_metadata
library(here)
library(tidyverse)
library(jsonlite)
library(future.apply)

# Define function for grabbing the biosample data of interest and associated attributes
get_biosample_df <- function(file) {
  attributes_oi <- c(
    "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "isolation_source", "lat_lon",
    "geographic location (latitude)", "geographic location (longitude)", "geographic location (region and locality)"
  )
  print(file)
  dat <- lapply(readLines(file), fromJSON)

  # from each of the items in dat, grab the biosample info, which is in the assemblyInfo attribute, add the taxon_id as we go
  biosamples <- list()
  for (i in 1:length(dat)) {
    if (!is.null(dat[[i]]$assemblyInfo$biosample$attributes)) {
      biosamples[[i]] <- dat[[i]]$assemblyInfo$biosample$attributes
      biosamples[[i]]$taxon_id <- dat[[i]]$taxId
      biosamples[[i]]$assembly_accession <- dat[[i]]$assemblyInfo$assemblyAccession
    }
  }

  biosample_df <- bind_rows(biosamples)

  # pivot the data
  biosample_filtered <- biosample_df %>%
    filter(name %in% attributes_oi) %>%
    pivot_wider(names_from = name, values_from = value, values_fn = list)

  cols <- colnames(biosample_filtered)
  cols_to_unnest <- cols[cols %in% attributes_oi]

  biosample_filtered <- biosample_filtered %>%
    unnest(cols = all_of(cols_to_unnest))

  return(biosample_filtered)
}

# Load the data
dat_loc <- here("data", "raw", "genomes", "ncbi_metadata", "100_first_taxons")

# Get list of files
dirs <- list.files(dat_loc, full.names = TRUE)
files <- file.path(dirs, "ncbi_dataset", "data", "assembly_data_report.jsonl")

# Get the data using multi-core lapply
plan(multisession)
biosample_dfs <- future_lapply(files, get_biosample_df)
final_df <- bind_rows(biosample_dfs)

# Save the data
write_csv(final_df, here("data", "intermediate", "ncbi_metadata_biosample_info.csv"))
