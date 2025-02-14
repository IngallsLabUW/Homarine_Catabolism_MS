# This script reads in the metadata from NCBI and cleans it up for plotting and analysis purposes.
# Input: CSV file containing the extracted biosample metadata information (output of ncbi_metadata_extract.R)
# Output: Cleaned and processed metadata information for downstream analysis - data frame with

# Libraries -----
library(tidyverse)
library(here)

# Special Functions-----
## function to parse lat_lon into lat and lon for plotting
parse_lat_lon <- function(lat_lon) {
  # split lat_lon into lat, dir, lon, dir2
  lat_lon_split <- str_split(lat_lon, " ")
  lat <- as.numeric(lat_lon_split[[1]][1])
  dir <- lat_lon_split[[1]][2]
  lon <- as.numeric(lat_lon_split[[1]][3])
  dir2 <- lat_lon_split[[1]][4]
  # if dir is S, make lat negative
  if (dir == "S") {
    lat <- -lat
  }
  # if dir2 is W, make lon negative
  if (dir2 == "W") {
    lon <- -lon
  }
  return(c(lat, lon))
}

# Get Data -----
dat <- read_csv(here("data", "intermediate", "ncbi_metadata_biosample_info.csv"),
  show_col_types = FALSE
) %>%
    mutate(isolation_source = tolower(isolation_source))

# get unique isolation_source with counts
unique_isolation_source <- dat %>%
  count(isolation_source) %>%
  mutate(perc = n / sum(n) * 100) %>%
  arrange(desc(n))

# NOTE: isolation_source is a free text field and contains a lot of variation in 
# the data, but is the most completely filled field in the dataset and is mostly 
# human-interpretable. We will use this field for our mapping of the biosamples 
# to different types

# Map data to categories ----
## Generate regex lookups -----
# These will be used as regex patterns to match the isolation_source field into 
# "environment", "non-environment", and "unknown" categories
unknown_isolation_source <- c(
  "missing", "not applicable", "unknown", "not collected", "not available",
  "not provided"
)
environment_isolation_source <- c(
  "soil", "root nodule", "seawater", "dock", "lake water", "sediment", 
  "groundwater", "marine", "estuary", "river", "environment", "water", 
  "ocean", "sea", "aquatic", "hydrothermal", "sludge", "root", "ocean", 
  "cold seep", "sand", "rhizhosphere",
  "rhizosphere", "coastal",
  "nodules", "grassland", "rock", "terrestrial", "lake", "aquifer", "compost",
  "phycosphere", "wetland", "woodland", "plankton", "hardwood", "hot spring",
  "tidal flat", "fjord", "lagoon", "glacier", "chesapeake bay", "microbial mat",
  "salt marsh", "stromatolite", "cave", "intertidal", "mud", "wga-x",
  "beach", "forest"
)
non_environment_isolation_source <- c(
  "patient", "feces", "stool", "food",
  "plastic particle", "clinical", "urban", "blood", "sputum", "urine",
  "wound", "intestines", "porcine", "human", "oral", "skin", "anal", "clinic",
  "gut", "metal", "cheese", "saliva", "intestinal", "intestine", "space station",
  "fracture fluid", "lpde", "bioreactor", "cleanroom floor", "rectal",
  "epithelial", "mucus", "fecal", "stomach", "chicken", "hospital",
  "digester", "groin", "laboratory", "milk", "ear", "polyp", "sewage",
  "mouth", "tongue", "bile", "occiput", "forearm", "lung", "rumen",
  "tooth", "clean room floor", "liquor", "urinary tract",
  "oyster", "fish", "shrimp", "clams", "tissue",
  "eye"
)

# These will be used as regex patterns to further categorize the environmental 
# samples into marine, freshwater, soil, and sediment
marine_islation_source <- c(
  "seawater",
  "marine",
  "estuary",
  "ocean",
  "coastal",
  "wga-x",
  "sea"
  )
freshwater_isolation_source <- c(
  "lake water",
  "groundwater",
  "river",
  "lake",
  "aquifer"
)
soil_isolation_source <- c(
  "soil",
  "root nodule",
  "sand",
  "rhizhosphere",
  "rhizosphere",
  "terrestrial",
  "compost",
  "root"
)
sediment_isolation_source <- c(
  "sediment"
)

## Apply lookups -----
## Generate lookup table for mapping isolation_source to categories
iso_source_mapped <- unique_isolation_source %>%
  mutate(isolation_source = tolower(isolation_source)) %>%
  mutate(isolation_source_mapped = case_when(
    is.na(isolation_source) ~ "unknown",
    str_detect(isolation_source, paste(non_environment_isolation_source, collapse = "|")) ~ "non-environmental",
    str_detect(isolation_source, paste(environment_isolation_source, collapse = "|")) ~ "environmental",
    str_detect(isolation_source, paste(unknown_isolation_source, collapse = "|")) ~ "unknown",
    TRUE ~ "unknown"
  )) %>%
  mutate(environmental_source_mapped = case_when(
      !str_detect(isolation_source_mapped, "environmental") ~ NA_character_,
    str_detect(isolation_source, paste(marine_islation_source, collapse = "|")) ~ "marine",
    str_detect(isolation_source, paste(freshwater_isolation_source, collapse = "|")) ~ "freshwater",
    str_detect(isolation_source, paste(soil_isolation_source, collapse = "|")) ~ "soil",
    str_detect(isolation_source, paste(sediment_isolation_source, collapse = "|")) ~ "sediment",
    TRUE ~ "other"
  ))

# Add resulting mappings
dat2 <- dat %>%
    left_join(
        iso_source_mapped %>%
            select(isolation_source, isolation_source_mapped, environmental_source_mapped),
        by = "isolation_source"
    )

# Add lat_lon for environmental samples -----
## Generate lat lon lookup table -----
lat_lon_lookup <- dat2 %>%
  select(lat_lon) %>%
  filter(!is.na(lat_lon)) %>%
  filter(str_detect(lat_lon, " N| S| W| E")) %>%
  distinct() %>%
  mutate(lat = NA, lon = NA)

for (i in 1:nrow(lat_lon_lookup)) {
  lat_lon <- lat_lon_lookup$lat_lon[i]
  lat_lon_parsed <- parse_lat_lon(lat_lon)
  lat_lon_lookup$lat[i] <- lat_lon_parsed[1]
  lat_lon_lookup$lon[i] <- lat_lon_parsed[2]
}

## Apply lat lon lookup  -----
dat3 <- dat2 %>%
  left_join(lat_lon_lookup, by = "lat_lon") %>%
    select(
        taxon_id, assembly_accession, isolation_source, lat_lon, 
        isolation_source_mapped, environmental_source_mapped, lat, lon
        )

# Save the data -----
write_csv(dat3, here("data", "intermediate", "ncbi_metadata_clean.csv"))
