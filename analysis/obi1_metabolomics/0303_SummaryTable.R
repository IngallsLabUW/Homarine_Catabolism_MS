# PURPOSE: provide a summary table of the data for a supplementary table for the paper


# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(here)
library(openxlsx2)


# SET FILE LOCATIONS  -----
meta_data_dir <- here("data", "intermediate", "metabolomics")
output_loc <- here("data", "intermediate", "metabolomics", "obi1")
sample_key <- here("data", "raw", "metabolomics", "obi1", "sample_key.csv")
sample_fraction_oi <- "Particulate"
save_loc <- here("tables")

# IMPORT DATA ----
raw_data <- here(output_loc, "combined_long_dat.csv") %>%
  read_csv(show_col_types = FALSE)
enrichment_data <- here(output_loc, "enrichment_summary.csv") %>%
  read_csv(show_col_types = FALSE)
mf_info <- here(output_loc, "combined_MF_info.csv") %>%
  read_csv(show_col_types = FALSE)
sample_key <- read_csv(sample_key, show_col_types = FALSE)

i = 5
for (sample_fraction_oi in c("Particulate", "Dissolved")) {
  # SUMMARIZE DATA ----
  # Prepare raw data (pivot to wide)
  supp_t <- raw_data %>%
    filter(sample_set == "OBi1_set2") %>%
    filter(sample_fraction == sample_fraction_oi) %>%
    select(mass_feature, replicate_name, adjusted_area) %>%
    pivot_wider(names_from = replicate_name, values_from = adjusted_area) %>%
    drop_na()

  # Add in mass feature info
  supp_t2 <- supp_t %>%
    left_join(mf_info %>%
      filter(sample_fraction == sample_fraction_oi), by = "mass_feature") %>%
    select(sample_fraction, mass_feature, z, RT, mz, MS2, add_annotation, starts_with("2"))

  # Add in enrichment data
  supp_t3 <- supp_t2 %>%
    left_join(enrichment_data, by = join_by(sample_fraction, mass_feature)) %>%
    relocate(core_metabolite, .after = mass_feature) %>%
    select(-med_pvalue_cmb, -med_qvalue_h, -med_qvalue_hNH4, -med_qvalue_cmb) %>%
    select(sample_fraction:add_annotation, starts_with("med"), starts_with("2")) %>%
    select(-sample_fraction) %>%
    arrange(z, mz)

  # Prepare the table for saving as an excel file
  supp_t4 <- supp_t3 %>%
    rename(
      `Mass Feature` = mass_feature,
      `Charge` = z,
      `Retention Time (min)` = RT,
      `m/z` = mz,
      `MS2` = MS2,
      `Adduct Annotation` = add_annotation,
      `Core Metabolite` = core_metabolite,
      `p (Glucose + NH4 + Homarine)` = med_pvalue_hNH4,
      `p (Homarine)` = med_pvalue_h,
      `fold change (Glucose + NH4 + Homarine)` = med_fc_hNH4,
      `fold change (Homarine)` = med_fc_h
    )


  # Prep sample_key
  sample_key <- sample_key %>%
    filter(replicate_name %in% colnames(supp_t)) %>%
    select(replicate_name, treatment)



  # Write it to an excel workbook
  sheet_name <- paste0("Obi1 ", sample_fraction_oi, " Metab Results")
  row_end <- nrow(supp_t4) + 3
  wb <- openxlsx2::wb_workbook() %>%
    wb_add_worksheet(sheet_name) %>%
    # Add the table
    wb_add_data(supp_t4, start_row = 3, sheet = sheet_name) %>%
    # Add the headers, merge the associated cells
    wb_add_data("Mass Feature Information", start_row = 2, start_col = 1, sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 1:7), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 1:7), sheet = sheet_name) %>%
    wb_add_data("Enrichment Data", start_row = 2, start_col = 8, sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 8:11), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 8:11), sheet = sheet_name) %>%
    wb_add_data("Adjusted areas", start_row = 1, start_col = 12, sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 1:1, cols = 12:20), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:row_end, cols = 12:20), sheet = sheet_name) %>%
    # Add treatments for the enrichment data
    wb_add_data("Glucose + NH4", start_row = 2, start_col = 12, sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 12:14), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:row_end, cols = 12:14), sheet = sheet_name) %>%
    wb_add_data("Homarine", start_row = 2, start_col = 15, sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 15:17), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:row_end, cols = 15:17), sheet = sheet_name) %>%
    wb_add_data("Glucose + NH4 + Homarine", start_row = 2, start_col = 18, sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 18:20), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:row_end, cols = 18:20), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:2, cols = 1:20), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:2, cols = 12:20), sheet = sheet_name) %>%
    # Set the column widths
    wb_set_col_widths(cols = c(1:5, 7:20), widths = "auto", sheet = sheet_name) %>%
    # Make bold and center top 2 rows and first column
    wb_add_font(dims = wb_dims(rows = 1:3, cols = 1:20), bold = "double", sheet = sheet_name) %>%
    wb_add_font(dims = wb_dims(rows = 1:row_end, cols = 1), bold = "single", sheet = sheet_name) %>%
    # Center align the top 2 rows
    wb_add_cell_style(
      wb_dims(rows = 1:2, cols = 1:20),
      horizontal = "center",
      sheet = sheet_name
    )
  # Add a border


  # Save the workbook
  save_path <- here(save_loc, paste0("Table_SX", i, "_Obi1_", sample_fraction_oi, "_metab_results.xlsx"))
  wb$save(save_path)
  i = i + 1
}
