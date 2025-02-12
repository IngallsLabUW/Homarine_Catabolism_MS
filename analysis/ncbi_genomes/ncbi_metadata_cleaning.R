# libraries
library(tidyverse)
library(here)

# get data
dat <- read_csv(here("data", "intermediate", "ncbi_metadata_biosample_info.csv"),
                show_col_types = FALSE)


# get unique lat_lon
unique_lat_lon <- dat %>%
  distinct(lat_lon) %>%
  filter(!is.na(lat_lon))

# parse lat_lon (from ##.### N ##.### W to lat and lon)
unique_lat_lon2 <- unique_lat_lon %>%
  separate(lat_lon, into = c("lat", "dir", "lon", "dir2"), sep = " ", remove = FALSE) %>%
    filter(dir == "N" | dir == "S")

# parse envo_broad_scale
unique_envo_broad_scale <- dat %>%
  distinct(env_broad_scale) %>%
  filter(!is.na(env_broad_scale)) %>%
    mutate(envo_term = str_extract(env_broad_scale, "ENVO.{1,2}\\d+"))

