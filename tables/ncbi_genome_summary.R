# Load required libraries ----
library(here)
library(tidyverse)
library(janitor)
library(openxlsx2)

# Load the data from Oscar of ncbi genomes with full pathway and assocaited taxonomy
dat <- read_delim(here("data", "intermediate", "ncbi_genome_counts", "complete_homarine_catabolic_operons_second_round_genome_metadata_and_taxonomy.tsv"), delim = "\t",
                         show_col_types = FALSE) %>%
    clean_names()
# Load scraped NCBI data
dat_isolation <- read_csv(here("data", "intermediate", "ncbi_metadata_clean.csv"),
                          show_col_types = FALSE) %>%
    rename(assembly = assembly_accession)
# Load homarine operon order
dat_homarine_order <- read_csv(here("tables", "Table_SX8_ncbi_genome_summary_mod.csv"),
                               show_col_types = FALSE) %>%
    rename(assembly = Assembly,
           nucleotide_accession = `Nucleotide Accesion`) %>%
    select(assembly, nucleotide_accession, gene_arrangement) %>%
    distinct() 

# Filter the data to only keep one nucleotide_accession per assembly
dat2 <- dat %>%
    group_by(assembly) %>%
    filter(nucleotide_accession == first(nucleotide_accession)) %>%
    ungroup() 

# Make wide format of assembly, nucleotide_accession, and homA, homB, homC, homD, homE start locations
dat_hom <- dat2 %>%
    filter(gene_name %in% c("homA", "homB", "homC", "homD", "homE")) %>%
    mutate(range = paste0(start, "-", end)) %>%
    select(assembly, nucleotide_accession, gene_name, range) %>%
    pivot_wider(names_from = gene_name, values_from = range) 

# Make long format of assembly, nucleotide_accession, kingdom, phylum, class, order, family, genus, species
dat_taxa <- dat2 %>%
    select(assembly, nucleotide_accession, kingdom, phylum, class, order, family, genus, species) %>%
    distinct() 


# Make long format of assembly and isolation_source informaiton
dat_iso <- dat_isolation %>%
    select(assembly, biosample_id, bioproject_accession, geo_loc_name, isolation_source,  isolation_source_mapped, environmental_source_mapped, lat, lon, notes) %>%
    distinct() %>%
    rename(
        isolation_source_standardized = isolation_source_mapped,
        environmental_source_standardized = environmental_source_mapped
    ) %>%
    # Clean up note - replace "derived from geo_loc_name" with "lat long derived from geographic location name"
    mutate(notes = ifelse(grepl("derived from geo_loc_name", notes), 
                          gsub("derived from geo_loc_name", "lat lon derived from geographic location name", notes), 
                          notes))

# Join the three dataframes together
dat_joined <- dat_hom %>%
    left_join(dat_taxa, by = join_by(assembly, nucleotide_accession)) %>%
    left_join(dat_iso, by = join_by(assembly)) %>%
    left_join(dat_homarine_order, by = join_by(assembly, nucleotide_accession)) %>%
    # move "gene_arrangement" to after "homE"
    relocate(gene_arrangement, .after = homE)

# Print out the % of alphaprotoeobacteria with homE
dat_joined %>%
    group_by(class) %>%
    summarise(homE = sum(!is.na(homE)),
              assembly_n = n_distinct(assembly)) %>%
    ungroup() %>%
    mutate(homE_per = homE / assembly_n*100)


# Rename the columns to be more informative
dat_joined <- dat_joined %>%
    rename(
        `Assembly Accession` = assembly,
        `Nucleotide Accesion` = nucleotide_accession,
        `Kingdom` = kingdom,
        `Phylum` = phylum,
        `Class` = class,
        `Order` = order,
        `Family` = family,
        `Genus` = genus,
        `Species` = species,
        `Biosample ID` = biosample_id,
        `Bioproject Accession` = bioproject_accession,
        `Geographic Location` = geo_loc_name,
        `Isolation Source` = isolation_source,
        `Isolation Source Standardized` = isolation_source_standardized,
        `Environmental Source Standardized` = environmental_source_standardized,
        `Latitude` = lat,
        `Longitude` = lon,
        `Notes` = notes,
        `operon arrangement` = gene_arrangement,
    ) %>%
    # Turn all NAs into empty strings
    mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

# Prepare to export to a nicely formatted excel file
sheet_name <- "NCBI_genome_summary"
row_end <- nrow(dat_joined) + 1 # Add 1 for the header row
col_end <- ncol(dat_joined)

# Construct the workbook
wb <- wb_workbook() %>%
    wb_add_worksheet(sheet_name) %>%
    
    # Add the data
    wb_add_data(dat_joined, sheet = sheet_name, 
                start_row = 2) %>%
    
    # Add header for Gene Location and merge cells
    wb_add_data("Gene Location", sheet = sheet_name, 
                start_col = 3, start_row = 1) %>%
    wb_merge_cells(dims = wb_dims(rows = 1:1, cols = 3:8), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 3:8), sheet = sheet_name) %>%
    
    # Add header for Taxonomy and merge cells
    wb_add_data("NCBI Taxonomy", sheet = sheet_name, 
                start_col = 9, start_row = 1) %>%
    wb_merge_cells(dims = wb_dims(rows = 1:1, cols = 9:15), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 9:15), sheet = sheet_name) %>%
    
    # Add header for Isolation and merge cells
    wb_add_data("NCBI Biosample Information", sheet = sheet_name, 
                start_col = 16, start_row = 1) %>%
    wb_merge_cells(dims = wb_dims(rows = 1:1, cols = 16:19), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 16:19), sheet = sheet_name) %>%
    
    # Add header for Standardized Isolation information and merge cells
    wb_add_data("Standardized Isolation Information", sheet = sheet_name, 
                start_col = 20, start_row = 1) %>%
    wb_merge_cells(dims = wb_dims(rows = 1:1, cols = 20:col_end), sheet = sheet_name) %>%
    wb_add_border(dims = wb_dims(rows = 1:row_end, cols = 20:col_end), sheet = sheet_name) %>%
    
    # Make bold and center top 2 rows and first column
    wb_add_font(dims = wb_dims(rows = 1:2, cols = 1:col_end), bold = "double", sheet = sheet_name) %>%
    
    # Set the column widths
    wb_set_col_widths(cols =  c(1:16, 20:col_end), widths = "auto", sheet = sheet_name) %>%
    wb_set_col_widths(cols = 17:19, widths = 20, sheet = sheet_name) %>%
    
    # Center align the top 2 rows
    wb_add_cell_style(
        wb_dims(rows = 1:2, cols = 1:col_end),
        horizontal = "center",
        sheet = sheet_name
    )

# Save the workbook
wb$save(here("tables", "ncbi_genome_summary.xlsx"))


