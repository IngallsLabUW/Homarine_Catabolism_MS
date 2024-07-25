# PURPOSE: calculate enrichment factors for the G4 samples

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(here)

# DATA IMPORT ----
# Import the data
import_dir <- here("data", "intermediate", "metabolomics", "g4", "d3d2_search_results")
output_dir <- here("data", "intermediate", "metabolomics", "g4")
meta_data_loc <- here("data", "raw", "metabolomics", "g4")
sample_key_filename <- here(meta_data_loc, "sample_key.csv")
sample_key <- read_csv(sample_key_filename, show_col_types = FALSE) %>%
  mutate(samp_name = replicate_name) %>%
  select(samp_name, treatment) %>%
  distinct()

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
  here(
    import_dir,
    "positive"
  ),
  pattern = ".csv", full.names = TRUE
)

## Run the function on all the files ----
positive_summary <- map_df(pos_files, summarise_d2d3_iso_results)

# NEGATIVE MODE ----
## Get list of negative files ----
neg_files <- list.files(
  here(
    import_dir,
    "negative"
  ),
  pattern = ".csv", full.names = TRUE
)

## Run the function on all the files ----
negative_summary <- map_df(neg_files, summarise_d2d3_iso_results)

# SUMMARIZE RESULTS ----
## Combine the positive and negative results ----
all_summary <- positive_summary %>%
  mutate(mode = "positive") %>%
  bind_rows(negative_summary %>% mutate(mode = "negative")) %>%
  mutate(isotopologue_type = if_else(is.na(isotopologue_type), "C12", isotopologue_type)) %>%
  left_join(sample_key, by = "samp_name")

all_summary2 <- all_summary %>%
  mutate(
    mf_id_samp = paste0(mf_id, "_", samp_name, "_", mode),
    mono_mf_id_samp = paste0(monoisotopic_mf_id, "_", samp_name, "_", mode)
  ) %>%
    mutate(
        homarine = str_detect(treatment, "HomFate")
    )

## Save the results ----
write_csv(all_summary2, here(output_dir, "all_iso_results.csv"))

## Group mass features by scan time, mode, and mz ----
grouped_summary <- all_summary2 %>%
  group_by(
    scan_time_round = round(scan_time, 0),
    mode,
    isotopologue_type,
    mz_round = round(mz, 1),
  ) %>%
  summarise(
    n = n(),
    n_homarine = sum(str_detect(treatment, "HomFate")),
    n_control = sum(str_detect(treatment, "Control_")),
    mz_average = mean(mz),
    scan_time_average = mean(scan_time),
    intensity_average = mean(intensity),
    samples = paste0(samp_name, collapse = ", ")
  ) %>%
  ungroup() %>%
  arrange(desc(n))

## Remove features that were in control samples-----
grouped_summary2 <- grouped_summary %>%
  filter(n_control == 0) %>%
    filter(isotopologue_type == "2H2" | isotopologue_type == "2H3")

## Calculate monoisotopic mass ----
grouped_summary2 <- grouped_summary2 %>%
  mutate(
    monoisotopic_mass = case_when(
      isotopologue_type == "2H2" ~ mz_average - (2.014101777844-1.007825031898)*2,
      isotopologue_type == "2H3" ~ mz_average - (2.014101777844-1.007825031898)*3
    )) %>%
    mutate(id = row_number())

## Save the results ----
write_csv(grouped_summary2, here(output_dir, "grouped_iso_results.csv"))

# MATCH TO INGALLS STANDARDS ----
## Pull in Ingalls standards to try to match ----
ingalls_standards <- read_csv(here(meta_data_loc, "Ingalls_Lab_Standards.csv"), show_col_types = FALSE) %>%
    clean_names() %>%
    filter(priority) %>%
    filter(column == "HILIC")

ingalls_standards2 <- ingalls_standards %>%
    mutate(
        mode = case_when(
            z == 1 ~ "positive",
            z == -1 ~ "negative"
        )
    ) %>%
    select(compound_name, mz, rt_minute, mode) %>%
    rename(mz_expected = mz, rt_expected = rt_minute) %>%
    distinct()


## Match the features to the Ingalls standards ----
grouped_summary_matched <- grouped_summary2 %>%
    left_join(ingalls_standards2, by = c("mode"), relationship = "many-to-many") %>%
    filter(
        abs(mz_expected - monoisotopic_mass) < 0.01,
        abs(rt_expected - scan_time_average) < 1.5
    ) %>%
    mutate(rt_diff = abs(rt_expected - scan_time_average)) %>%
    group_by(id) %>%
    filter(rt_diff == min(rt_diff)) %>%
    ungroup() %>%
    select(id, everything()) %>%
    mutate(
        compound_name = case_when(
            compound_name == "L-2-Aminoadipic acid" ~ "N-methyl glutamic acid",
            TRUE ~ compound_name
        )
    ) %>%
    select(id, n, mode, mz_average, scan_time_average, isotopologue_type, samples, compound_name, rt_expected, mz_expected)

## Save the results ----
write_csv(grouped_summary_matched, here(output_dir, "grouped_iso_results_matched.csv"))

