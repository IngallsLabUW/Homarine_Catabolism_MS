# PURPOSE: combine the transcriptomics results for manuscript

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(here)
library(openxlsx2)

# DATA IMPORT ----
# read in xslx file from Table SX. Homarine_vs_Glucose_Obi1.xlsx
homarine_vs_glucose <- read_xlsx(here("data", "intermediate", "transcriptomics", "Table SX. Homarine_vs_Glucose_Obi1.xlsx"), startRow = 2)
homarine_plus_glucose_vs_glucose <- read_xlsx(here("data",
                                                    "intermediate",
                                                    "transcriptomics",
                                                    "Table SX. Homarine_Glucose_vs_Glucose_Obi1.xlsx"),
                                            startRow = 2)

# prepare data for combining in one excel
dat <- homarine_vs_glucose %>%
    select(gene.id,
           strand,
           name,
           COG,
           product,
           KO,
           `KO Definition`,
           baseMean,
           log2FoldChange,
           lfcSE,
           svalue)
dat2 <- homarine_plus_glucose_vs_glucose %>%
    select(gene.id,
           baseMean,
           log2FoldChange,
           lfcSE,
           svalue) %>%
    # rename all except gene.id with _hg tag at end
    rename_with(~str_c(., "_hg"), -gene.id)

# COMBINE DATA ----
# Need to join to make the columns match up
dat <- dat %>%
    full_join(dat2, by = "gene.id")


# WRITE TO EXCEL ----
sheet_name = "Table SX3"
row_end <- nrow(dat) + 3
wb <- openxlsx2::wb_workbook() %>%
    wb_add_worksheet(sheet_name) %>%
    # Add data to the table for Homarine vs Glucose
    wb_add_data(dat %>% select(gene.id:svalue), start_row = 3, sheet = sheet_name) 

# Prepare data for the Homarine+Glucose vs Glucose
dat_sub <- dat %>% select(baseMean_hg:svalue_hg) %>%
    rename_with(~str_remove(., "_hg"))

# Add data to the table for Homarine+Glucose vs Glucose
wb <- wb %>%
    wb_add_data(dat_sub, start_row = 3, start_col = 12, sheet = sheet_name)

# Add headers and formatting
wb <- wb %>%
    # Add header for Gene Annotation and merge cells
    wb_add_data("Gene Annotation", 
                start_row = 2, start_col = 1, 
                sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 1:7), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 1:7), sheet = sheet_name) %>%
    wb_add_data(
        "Homarine vs Glucose statistics", 
        start_row = 2, start_col = 8, 
        sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 8:11), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 8:11), sheet = sheet_name) %>%
    wb_add_data(
        "Homarine + Glucose vs Glucose statistics", 
        start_row = 2, start_col = 12, 
        sheet = sheet_name) %>%
    wb_merge_cells(dims = wb_dims(rows = 2:2, cols = 12:15), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 12:15), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:2, cols = 1:7), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:2, cols = 8:11), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 2:2, cols = 12:15), sheet = sheet_name) %>%
    # Set the column widths
    wb_set_col_widths(cols = c(1:4, 6, 8:15), widths = "auto", sheet = sheet_name) %>%
    # Make bold and center top 2 rows and first column
    wb_add_font(dims = wb_dims(rows = 1:3, cols = 1:15), bold = "double", sheet = sheet_name) %>%
    # Center align the top 2 rows
    wb_add_cell_style(
        wb_dims(rows = 1:2, cols = 1:20),
        horizontal = "center",
        sheet = sheet_name
    )

save_path <- here("tables", "Table_SX2_OBi1_transcript_full_results.xlsx")
wb$save(save_path)
