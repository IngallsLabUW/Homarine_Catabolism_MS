# PURPOSE: combine and tidy the untargeted data with the targeted data
# 20240313  KRH moving from Obi1-specific repo

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(here)

# SET FILE LOCATIONS  -----
output_loc <- here("data", "intermediate", "metabolomics", "obi1")

# Read in files  -----
dat_target_filename <- here(
  output_loc,
  "targeted",
  "combined_tidy_dat_long.csv"
)
dat_nontarget_filename <- here(
  output_loc,
  "untargeted",
  "combined_tidy_dat_long.csv"
)

dat_target <- read_csv(dat_target_filename, show_col_types = FALSE) %>%
  rename(RT = retention_time)
dat_nontarget <-
  read_csv(dat_nontarget_filename, show_col_types = FALSE) %>%
  mutate(z = case_when(
    grepl("Pos", mass_feature) ~ 1,
    grepl("Neg", mass_feature) ~ -1
  ))

# MUNGE DATA ----
## Combine data -----
dat_cmb <- dat_target %>%
  mutate(dat_type = "targeted") %>%
  bind_rows(dat_nontarget %>%
    mutate(dat_type = "nontargeted")) %>%
  filter(sample_set == "OBi1_set2")

## make intermediate of MF info -----
dat_MF_2 <- dat_cmb %>%
  select(
    mass_feature, z, sample_fraction, sample_set, dat_type, RT, mz, MS2,
    add_annotation
  ) %>%
  distinct() %>%
  mutate(row_id = row_number())
dat_cmb2 <- dat_cmb %>%
  left_join(dat_MF_2, by = join_by(mass_feature, sample_fraction, z, sample_set, RT, mz, dat_type, MS2, add_annotation))


# DEREPLICATE UNTARGETED FROM TARGETED------
# remove from targeted and rename the mass feature in untargeted to the targeted mass feature's name
## first do particulate
target_MFs <- dat_MF_2 %>%
  filter(dat_type == "targeted")
untarget_MFs <- dat_MF_2 %>%
  filter(dat_type == "nontargeted")
names_to_change <- list()
to_drop <- c()
for (i in 1:nrow(target_MFs)) {
  # for each targeted MF
  target_oi <- target_MFs[i, ]
  # filter untarget_MF for mz within 0.01 of target_oi mz and RT within 0.1 of target_oi RT
  untarget_oi <- untarget_MFs %>%
    filter(
      abs(mz - target_oi$mz) < 0.005,
      abs(RT - target_oi$RT) < 0.25,
      z == target_oi$z,
      sample_fraction == target_oi$sample_fraction
    )
  # if there is a match, drop the untargeted MF
  if (nrow(untarget_oi) > 0) {
    untarget_oi <- untarget_oi %>%
      mutate(
        mass_feature = target_oi$mass_feature,
        dat_type = "targeted"
      )
    names_to_change[[i]] <- untarget_oi
    to_drop <- c(to_drop, target_oi$row_id)
  }
}

# drop untargeted MFs that are dereplicated
new_names <- bind_rows(names_to_change) %>%
  group_by(row_id) %>%
  filter(row_number() == 1)
dat_MF_update <- dat_MF_2 %>%
  filter(!(row_id %in% to_drop)) %>%
  filter(!(row_id %in% new_names$row_id)) %>%
  bind_rows(new_names)
dat_cmb3 <- dat_cmb2 %>%
  filter(!(row_id %in% to_drop)) %>%
  left_join(new_names %>%
    select(mass_feature, sample_fraction, row_id, dat_type) %>%
    rename(
      new_mass_feature = mass_feature,
      new_dat_type = dat_type
    )) %>%
  mutate(
    mass_feature = ifelse(!is.na(new_mass_feature), new_mass_feature, mass_feature),
    dat_type = ifelse(!is.na(new_dat_type), new_dat_type, dat_type)
  ) %>%
  select(-new_mass_feature, -new_dat_type)

# WRITE DATA -----
write_csv(dat_cmb3, here(output_loc, "combined_long_dat.csv"))
write_csv(dat_MF_update, here(output_loc, "combined_MF_info.csv"))
