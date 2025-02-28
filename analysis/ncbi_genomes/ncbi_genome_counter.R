library(httr)
library(xml2)
library(dplyr)
library(jsonlite)

# List of Order names to query
orders <- c("Cellvibrionales", "Chromatiales", "Pseudomonadales")  # Replace with your list of Order names

# Function to get taxon ID from Entrez API
get_taxon_id <- function(order_name) {
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=", order_name, "[Order]&retmode=xml")
    response <- GET(url)
    
    # Parse XML response
    xml_content <- content(response, "text")
    parsed_xml <- read_xml(xml_content)
    
    # Extract taxon ID from XML
    taxon_id <- xml_find_first(parsed_xml, ".//Id") %>% xml_text()
    
    return(taxon_id)
}

# Function to get number of assemblies from NCBI Datasets API
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

# Loop through each order, get taxon ID and number of assemblies
results <- tibble(Order = character(), Taxon_ID = character(), Number_of_Assemblies = integer())

for (order in orders) {
    cat("Processing:", order, "\n")
    
    # Get taxon ID from Entrez API
    taxon_id <- get_taxon_id(order)
    if (taxon_id == "") {
        cat("No Taxon ID found for", order, "\n")
        next
    }
    
    # Get number of assemblies from Datasets API
    assemblies_count <- get_assemblies_count(taxon_id)
    
    # Store the result
    results <- results %>%
        add_row(Order = order, Taxon_ID = taxon_id, Number_of_Assemblies = assemblies_count)
}

# View the results
print(results)
