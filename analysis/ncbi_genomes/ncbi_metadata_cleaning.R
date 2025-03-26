# This script reads in the metadata from NCBI and cleans it up for plotting and analysis purposes.
# Input: CSV file containing the extracted biosample metadata information (output of ncbi_metadata_extract.R)
# Output: Cleaned and processed metadata information for downstream analysis - data frame with

# Libraries -----
library(tidyverse)
library(here)
library(httr)
library(jsonlite)

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

## function to convert lat or lon from DMS to decimal
convert_lat_lon <- function(lat) {
  # split lat_lon into lat, dir, lon, dir2
  lat_split <- str_split(lat, " ")
  dir <- lat_split[[1]][1]
  lat_deg <- parse_number(lat_split[[1]][2])
  lat_min <- parse_number(lat_split[[1]][3])
  lat_sec <- parse_number(lat_split[[1]][4])
  lat_dec <- lat_deg + lat_min / 60 + lat_sec / 3600
  if (dir == "S") {
    lat_dec <- -lat_dec
  }
  if (dir == "W") {
    lat_dec <- -lat_dec
  }
  return(lat_dec)
}
 
# Read in and parse lookup -------
geo_loc_lookup <- read_csv(here("data", "intermediate", "ncbi_metadata_geo_loc_name_lookup.csv"),
  show_col_types = FALSE
) %>%
    filter(!is.na(lat) & !is.na(lon)) %>%
    select(geo_loc_name, lat, lon) %>%
    mutate(lat_new = NA_real_, lon_new = NA_real_)

# parse lat_lon in lookup using parse_lat_lon2
for (i in 1:nrow(geo_loc_lookup)) {
  lat_new <- convert_lat_lon(geo_loc_lookup$lat[i])
  lon_new <- convert_lat_lon(geo_loc_lookup$lon[i])
  geo_loc_lookup$lat_new[i] <- lat_new
  geo_loc_lookup$lon_new[i] <- lon_new
}
geo_loc_lookup <- geo_loc_lookup %>%
  select(-lat, -lon) %>%
  rename(lat = lat_new, lon = lon_new) %>%
    mutate(location_qual = "derived from geo_loc_name")

# Get Data -----
dat <- read_csv(here("data", "intermediate", "ncbi_metadata_biosample_info.csv"),
  show_col_types = FALSE
) %>%
    mutate(isolation_source = tolower(isolation_source))

genomes_oi <- read_delim(here("data", "intermediate", "ncbi_genome_counts", "ncbi_complete_homarine_operons_loci_with_taxonomy.txt"), delim = "\t",
  show_col_types = FALSE
)

## filter data to only those in genomes_oi
dat <- dat %>%
  filter(assembly_accession %in% genomes_oi$Assembly)

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
    str_detect(isolation_source, paste(environment_isolation_source, collapse = "|")) ~ "environmental",
    str_detect(isolation_source, paste(non_environment_isolation_source, collapse = "|")) ~ "non-environmental",
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

# Map bioproject to categories ----
## Read in bioproject lookup table -----
bioproject_lookup <- read_csv(here("data", 
                                   "intermediate",
                                   "ncbi_genome_counts",
                                   "ncbi_bioproject_lookup.csv"),
                              show_col_types = FALSE)

# update records with bioproject-level annotations of interest
dat2 <- dat2 %>%
    rows_update(bioproject_lookup, by = "bioproject_accession")
  
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
  left_join(lat_lon_lookup, by = "lat_lon")

## Use and apply `geographic location (latitude)`, `geographic location (longitude)`
# if lat_lon is missing
# Go from ~30% coverage to 37% coverage
dat3 <- dat3 %>%
  mutate(lat = case_when(
    !is.na(lat) ~ lat,
    !is.na(`geographic location (latitude)`) ~ as.numeric(`geographic location (latitude)`),
    TRUE ~ NA_real_
  )) %>%
  mutate(lon = case_when(
    !is.na(lon) ~ lon,
    !is.na(`geographic location (longitude)`) ~ as.numeric(`geographic location (longitude)`),
    TRUE ~ NA_real_
  )) %>%
    # if lat is out of range, replace with NA
    mutate(lat = case_when(
      lat < -90 | lat > 90 ~ NA_real_,
      lon < -180 | lon > 180 ~ NA_real_,
      TRUE ~ lat
    )) %>%
    # if lon is out of range, replace with NA
    mutate(lon = case_when(
      lon < -180 | lon > 180 ~ NA_real_,
      lat < -90 | lat > 90 ~ NA_real_,
      TRUE ~ lon
    ))

## Apply lat lon lookup from geo_loc_name -----
### Pull out records with missing lat_lon and attempt to fill in with geo_loc_name lookup
dat3_no_lat_lon <- dat3 %>%
  filter(is.na(lat) | is.na(lon)) %>%
    select(-lat, -lon) %>%
    left_join(geo_loc_lookup, by = "geo_loc_name")

# nearly 80% coverage for lat lon!
dat4 <- dat3 %>%
    filter(!is.na(lat) | !is.na(lon)) %>%
    bind_rows(dat3_no_lat_lon)

## Summarize data -----
dat_summary <- dat4 %>%
  group_by(isolation_source_mapped, environmental_source_mapped) %>%
  summarize(n = n(),
            n_with_lat_lon = sum(!is.na(lat) & !is.na(lon)),
            n_with_lat_lon_pct = n_with_lat_lon / n(),
            ) %>%
  arrange(desc(n))

# Get a df of just the environmental samples -----
#only about 30% of these have lat_lon data :(
dat5 <- dat4 %>%
    filter(isolation_source_mapped == "environmental")

# Make a plot of dat4, colored by environmental_source_mapped
dat5 %>%
    ggplot(aes(x = lon, y = lat, color = environmental_source_mapped)) +
    geom_point() +
    theme_minimal()
dat4 %>%
    ggplot(aes(x = lon, y = lat, color = environmental_source_mapped)) +
    geom_point() +
    theme_minimal()


# PROJECT NUMBERS OF INTEREST: 
#PRJNA596592, n = 642 "Genomic analysis of coral associated bacteria in Hong Kong waters" Assign lat long from Hong Kong: Wong Wan Chau, with note "Hong Kong waters, lat long derived from biopackage information"
#PRJNA391943 (tara, metagenome assembled genomes); PRJNA417962 (metagenome-assembled genomed (marine atlas project))



#%>%
#    select(
#        taxon_id, assembly_accession, isolation_source, lat_lon, 
#        isolation_source_mapped, environmental_source_mapped, lat, lon
#        )

# Save the data -----
write_csv(dat3, here("data", "intermediate", "ncbi_metadata_clean.csv"))
