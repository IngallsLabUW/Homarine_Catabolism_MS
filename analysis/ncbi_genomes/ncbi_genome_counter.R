library(httr)
library(xml2)
library(dplyr)
library(jsonlite)
library(tidyverse)
library(here)

# Define your functions ------
# Function to get taxon ID from Entrez API
get_taxon_id <- function(name, taxon_level) {
  # convert order_name to url - encoded order name
  name <- URLencode(name)
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=", name, "[", taxon_level, "]&retmode=xml")
  response <- GET(url)

  # Parse XML response
  xml_content <- content(response, "text")
  parsed_xml <- read_xml(xml_content)

  # Extract taxon ID from XML
  taxon_id <- xml_find_first(parsed_xml, ".//Id") %>% xml_text()

  return(taxon_id)
}

# Function to get number of assemblies from NCBI Datasets API
#https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/biosample/-biosample_ids-/dataset_report
get_assemblies_count <- function(taxon_id) {
  url <- paste0("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/", taxon_id)
  response <- GET(url)

  # Parse JSON response
  json_content <- content(response, "text", encoding = "UTF-8")
  parsed_json <- fromJSON(json_content)

  # Extract number of assemblies
  assemblies_count <- parsed_json$taxonomy_nodes$taxonomy$counts %>%
    as.data.frame() %>%
    filter(type == "COUNT_TYPE_ASSEMBLY") %>%
    pull(count)

  return(assemblies_count)
}

# Function to get number of genomes from a name and taxon level,
# using the get_taxon_id and get_assemblies_count functions
get_genomes_count <- function(name, taxon_level) {
  taxon_id <- get_taxon_id(name, taxon_level)
  assemblies_count <- get_assemblies_count(taxon_id)

  return(assemblies_count)
}

# Get data --------
dat <- read_delim(
  here(
    "data",
    "intermediate",
    "ncbi_genome_counts",
    "ncbi_homarine_operons_loci_with_taxonomy.txt"
  ),
  show_col_types = FALSE
)

# Find the unique genomes per phylum -------
phyla <- dat %>%
  filter(!is.na(phylum)) %>%
  select(phylum) %>%
  distinct() %>%
  mutate(genomes = map_dbl(phylum, get_genomes_count, taxon_level = "Phylum"))
write_csv(phyla, here("data", "intermediate", "ncbi_genome_counts", "ncbi_genomes_per_phylum.csv"))

# Find the unique genomes per class -------
classes <- dat %>%
  filter(!is.na(class)) %>%
  select(class) %>%
  distinct() %>%
  mutate(genomes = map_dbl(class, get_genomes_count, taxon_level = "Class"))
write_csv(classes, here("data", "intermediate", "ncbi_genome_counts", "ncbi_genomes_per_class.csv"))

# Find the unique genomes per order -------
orders <- dat %>%
  filter(!is.na(order)) %>%
  select(order) %>%
  distinct() %>%
  mutate(genomes = map_dbl(order, get_genomes_count, taxon_level = "Order"))
write_csv(orders, here("data", "intermediate", "ncbi_genome_counts", "ncbi_genomes_per_order.csv"))

# Find the unique genomes per family -------
families <- dat %>%
  filter(!is.na(family)) %>%
  select(family) %>%
  distinct() %>%
  mutate(genomes = map_dbl(family, get_genomes_count, taxon_level = "Family"))
write_csv(families, here("data", "intermediate", "ncbi_genome_counts", "ncbi_genomes_per_family.csv"))
