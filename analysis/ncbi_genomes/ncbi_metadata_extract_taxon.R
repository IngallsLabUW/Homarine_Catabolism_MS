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
  taxon_id <- parse_number(basename(dirname(dirname(dirname(file)))))

  print(paste("Processing:", file))

  dat <- lapply(readLines(file), fromJSON)
  biosamples <- list()
  for (i in seq_along(dat)) {
    if (is.list(dat[[i]]$assemblyInfo$biosample)) {
        # First get attributes 
      if (!is.null(dat[[i]]$assemblyInfo$biosample$attributes)) {
        biosamples[[i]] <- dat[[i]]$assemblyInfo$biosample$attributes
      } else {
        biosamples[[i]] <- data.frame(
          name = NA,
          value = NA
        )
      }
        # Next gather metadata about the assembly
        biosamples[[i]]$file <- file
        biosamples[[i]]$taxon_id <- taxon_id
        #browser()
        if (!is.null(dat[[i]]$accession)) {
            biosamples[[i]]$assembly_accession <- dat[[i]]$accession
        } else if (!is.null(dat[[i]]$assemblyInfo$assemblyAccession)){
            biosamples[[i]]$assembly_accession <- dat[[i]]$assemblyInfo$assemblyAccession
        } else {
            biosamples[[i]]$assembly_accession <- NA
            warning("No assembly accession found for file: ", file)
        }

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
  biosample_df <- bind_rows(biosamples)
  # Filter only for the attributes of interest
    biosample_attributes <- biosample_df %>%
        select(name, value, assembly_accession) %>%
        filter(name %in% attributes_oi)%>%
        pivot_wider(names_from = name, values_from = value, values_fn = list)

  biosample_filtered <- biosample_df %>%
      select(-name, -value) %>%
      distinct() %>%
      left_join(biosample_attributes, by = "assembly_accession")

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
dat_loc <- normalizePath("data/raw/genomes/ncbi_metadata/taxon_downloads_2/")
dataset_dirs <- list.dirs(dat_loc, recursive = FALSE)
dataset_dirs <- dataset_dirs[grepl("dataset_\\d+$", dataset_dirs)] # Keep only dataset directories

#file <- "/Users/heal742/LOCAL/07_UW/Homarine_Catabolism_MS/data/raw/genomes/ncbi_metadata/taxon_downloads_2/ncbi_dataset_303/ncbi_dataset/data/assembly_data_report.jsonl"
#test_df <- get_biosample_df(file)

files <- file.path(dataset_dirs, "ncbi_dataset", "data", "assembly_data_report.jsonl")
files <- files[file.exists(files)]

# Process the JSONL files using parallel processing ----
plan(multisession)
biosample_dfs <- future_lapply(files, get_biosample_df)

# Combine results into a single data frame ----
biosample_data <- bind_rows(biosample_dfs, .id = "dataset")

# Read in the genomes of interest and filter out on accession id
genomes_oi <- read_delim(here("data", "intermediate", "ncbi_genome_counts", "complete_homarine_catabolic_operons_second_round_genome_metadata_and_taxonomy.tsv"), delim = "\t",
                         show_col_types = FALSE) %>%
    rename(assembly_accession = Assembly,
           biosample_id = BioSample)

biosample_data_sub <- biosample_data %>%
    filter(assembly_accession %in% genomes_oi$assembly_accession) 

# Save the output ----
write_csv(biosample_data_sub, here("data", "intermediate", "ncbi_metadata_biosample_info.csv"))

# Show completion message
print(paste("Processing complete. Output saved to:", here("data", "intermediate", "ncbi_metadata_biosample_info.csv")))

