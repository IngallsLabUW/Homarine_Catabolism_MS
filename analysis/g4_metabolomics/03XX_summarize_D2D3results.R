# PURPOSE: calculate enrichment factors for the G4 samples

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(here)

# DATA IMPORT ----
# Import the data
import_dir <- here("data", "intermediate", "metabolomics", "g4", "d3d2_search_results_thresh1")
output_dir <- here("data", "intermediate", "metabolomics", "g4")

# FUNCTION TO SUMMARIZE RESULTS ----
summarise_d2d3_iso_results <- function(file) {
    full_data <- read_csv(file, show_col_types = FALSE)
    samp_name <- str_remove(basename(file), ".csv") %>%
        str_remove("_isos$")
    D2_monos <- full_data %>%
        filter(isotopologue_type == "2H2") %>%
        select(monoisotopic_mf_id) %>%
        pull()
    D3_monos <- full_data %>%
        filter(isotopologue_type == "2H3") %>%
        select(monoisotopic_mf_id) %>%
        pull()
    data_subset <- full_data %>%
        filter(monoisotopic_mf_id %in% c(D2_monos, D3_monos)) %>%
        mutate(samp_name = samp_name)
    return(data_subset)
}

# POSITIVE MODE ----

## Get list of positive files ----
pos_files <- list.files(
    here(import_dir,
         "positive"),
    pattern = ".csv", full.names = TRUE)

## Run the function on all the files ----
positive_summary <- map_df(pos_files, summarise_d2d3_iso_results)

## Save the results ----
write_csv(positive_summary, here(output_dir, "g4_d2d3_summary_positiveMFs.csv"))

# NEGATIVE MODE ----
## Get list of negative files ----
neg_files <- list.files(
    here(import_dir,
         "negative"),
    pattern = ".csv", full.names = TRUE)

## Run the function on all the files ----
negative_summary <- map_df(neg_files, summarise_d2d3_iso_results)

## Save the results ----
write_csv(negative_summary, here(output_dir, "g4_d2d3_summary_negativeMFs.csv"))





