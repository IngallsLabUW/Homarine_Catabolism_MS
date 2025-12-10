library(dplyr)
library(tidyr)
library(readr)

setwd("C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project")  

# Load your file
df <- read_csv("ncbi_genomes_summary.csv")

# Extract start and end for homA and homB
df_lengths <- df %>%
  separate(homA, into = c("homA_start", "homA_end"), sep = "-", convert = TRUE) %>%
  separate(homB, into = c("homB_start", "homB_end"), sep = "-", convert = TRUE) %>%
  mutate(
    homA_length = homA_end - homA_start + 1,
    homB_length = homB_end - homB_start + 1
  )

# Summaries
summary_lengths <- df_lengths %>%
  summarise(
    homA_mean = mean(homA_length, na.rm = TRUE),
    homA_sd   = sd(homA_length, na.rm = TRUE),
    homB_mean = mean(homB_length, na.rm = TRUE),
    homB_sd   = sd(homB_length, na.rm = TRUE)
  )

summary_lengths
# # A tibble: 1 × 4
# homA_mean homA_sd homB_mean homB_sd
# <dbl>   <dbl>     <dbl>   <dbl>
#   1      704.    63.3     1721.    89.9