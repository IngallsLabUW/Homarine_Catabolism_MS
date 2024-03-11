# _____________________________________________________________________________
#       Copyright 2022, Integral Consulting Inc. All rights reserved.        #
# _____________________________________________________________________________
#
# _____________________________________________________________________________
#  0103_Dissolved_QC.R
# _____________________________________________________________________________
# PURPOSE:
#  Perform post-BMIS QC on dissolved sampls
#
#
# PROJECT INFORMATION:
#   Name:
#   Number: C3352
#
# HISTORY:
# 	 Date		        Remarks
# _____________________________________________________________________________
# 	 20220505       Using this as a guide: https://github.com/jssacks/CX_SPE_Method/blob/main/R_Code/EX_CX_HILIC_Blk_QC.R per J. Sacks' suggestion
#

# PACKAGES ----
library(tidyverse)
library(here)

# Read in file  -----
dat_filename <- "01_targeted/data_intermediate/BMISd_HILIC_Dissolved.csv"

dat <- read_csv(dat_filename, show_col_types = FALSE)

# Calculate batch blank limit of detection -----
blk_ave_dat <- dat %>%
  filter(replicate_name %>%
           str_detect("MQBlk")) %>%
 # filter(!is.na(adjusted_area)) %>%
  group_by(mass_feature) %>%
  summarise(
    blk_n = n(),
    blk_ave = mean(area,na.rm = TRUE),
            blk_sd = sd(area,na.rm = TRUE),
            blk_max = max(area,na.rm = TRUE)) %>%
  mutate(blk_lod = blk_ave + (3.182 * (blk_sd/sqrt(9))))

# Overwrite areas when blk LOD < adjusted_area LOD ------
dat_wdissblanks <- dat %>%
  left_join(blk_ave_dat %>%
              select(mass_feature,
                     blk_lod)) %>%
  mutate(adjusted_area = ifelse(adjusted_area<blk_lod, NA, adjusted_area),
         QC_area = ifelse(QC_area<blk_lod, NA, QC_area))

# Write out results -----
write_csv(dat_wdissblanks, "01_targeted/data_intermediate/BMISd_HILIC_Dissolved_postDissBlk.csv")


