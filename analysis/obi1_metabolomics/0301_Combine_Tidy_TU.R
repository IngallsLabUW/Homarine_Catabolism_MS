# PURPOSE: combine and tidy the untargeted data with the targeted data
# 20240313  KRH moving from Obi1-specific repo

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(viridis)
library(broom)
library(coin)
library(here)

# SET FILE LOCATIONS  -----
output_loc <- here("data", "intermediate", "metabolomics", "obi1")

# Read in files  -----
dat_target_filename <- here(
    output_loc,
    "targeted",
    "combined_tidy_dat_long.csv"
)
dat_nontarget_filename <-here(
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
  filter(sample_set == "OBi1_set2") %>%
  filter(sample_fraction == "Particulate")

## Calculate fold changes and select median fold change for scaling ------
dat_fc <- dat_cmb %>%
  group_by(mass_feature, treatment, sample_fraction) %>%
  summarise(mean_area = mean(adjusted_area, rm.na = T)) %>%
  ungroup()%>%
  group_by(mass_feature, sample_fraction) %>%
  summarise(
    fc_hNH4 = mean_area[treatment == "Glucose + NH4 + Homarine"] /
      mean_area[treatment == "Glucose + NH4"],
    fc_h = mean_area[treatment == "Homarine"] / mean_area[treatment == "Glucose + NH4"]
  ) %>%
  ungroup()

scaler_h <- dat_fc %>%
  group_by(sample_fraction) %>%
  mutate(diff_from_med_h = abs(fc_h - median(fc_h, na.rm = T))) %>%
  filter(diff_from_med_h == min(diff_from_med_h, na.rm = T)) %>%
  ungroup() %>%
  head(1) %>%
  pull(fc_h)

scaler_hNH4 <- dat_fc %>%
  group_by(sample_fraction) %>%
  mutate(diff_from_med_hNH4 = abs(fc_hNH4 - median(fc_hNH4, na.rm = T))) %>%
  filter(diff_from_med_hNH4 == min(diff_from_med_hNH4, na.rm = T)) %>%
  ungroup() %>%
  head(1) %>%
  pull(fc_hNH4)

## scale data ----
dat_cmb2 <- dat_cmb %>%
  mutate(scaled_area = case_when(
    treatment == "Homarine" ~ adjusted_area / scaler_h,
    treatment == "Glucose + NH4 + Homarine" ~ adjusted_area / scaler_hNH4,
    treatment == "Glucose + NH4" ~ adjusted_area
  ))


## summarize abundance features ----
high_homarine_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area)) %>%
  filter(treatment != "Glucose + NH4 + Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area))) %>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Homarine" & area_rank_sum == 6) %>%
  pull(mass_feature)

high_homarine_NH4_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area)) %>%
  filter(treatment != "Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area))) %>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Glucose + NH4 + Homarine" & area_rank_sum == 6) %>%
  pull(mass_feature)

low_homarine_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area)) %>%
  filter(treatment != "Glucose + NH4 + Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area))) %>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Homarine" & area_rank_sum == 15) %>%
  pull(mass_feature)

low_homarine_NH4_compounds <- dat_cmb2 %>%
  filter(!is.na(scaled_area)) %>%
  filter(treatment != "Homarine") %>%
  group_by(mass_feature) %>%
  mutate(area_rank = rank(desc(scaled_area))) %>%
  ungroup() %>%
  group_by(mass_feature, treatment) %>%
  summarise(area_rank_sum = sum(area_rank)) %>%
  ungroup() %>%
  filter(treatment == "Glucose + NH4 + Homarine" & area_rank_sum == 15) %>%
  pull(mass_feature)

high_both_compounds <- high_homarine_compounds[high_homarine_compounds %in% high_homarine_NH4_compounds]
low_both_compounds <- low_homarine_compounds[low_homarine_compounds %in% low_homarine_NH4_compounds]

dat_cmb2 <- dat_cmb2 %>%
  mutate(
    high_homarine_compound = mass_feature %in% high_homarine_compounds,
    low_homarine_compound = mass_feature %in% low_homarine_compounds,
    high_homarine_NH4_compound = mass_feature %in% high_homarine_NH4_compounds,
    low_homarine_NH4_compound = mass_feature %in% low_homarine_NH4_compounds,
    high_both_compound = mass_feature %in% high_both_compounds,
    low_both_compounds = mass_feature %in% low_both_compounds
  )

## make intermediate of MF info -----
dat_MF_2 <- dat_cmb2 %>%
  select(
    mass_feature, z, sample_fraction, sample_set, dat_type, RT, mz, MS2,
    add_annotation, high_homarine_compound:low_both_compounds
  ) %>%
  distinct() %>%
    mutate(row_id = row_number())
dat_cmb2 <- dat_cmb2 %>%
  left_join(dat_MF_2)


# DEREPLICATE UNTARGETED FROM TARGETED------
# remove from targeted and rename the mass feature in untargeted to the targeted mass feature's name
target_MFs <- dat_MF_2 %>%
  filter(dat_type == "targeted")
untarget_MFs <- dat_MF_2 %>%
  filter(dat_type == "nontargeted")
names_to_change = list()
to_drop <- c()
for (i in 1:nrow(target_MFs)) {
  # for each targeted MF
    target_oi <- target_MFs[i,]
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
            mutate(mass_feature = target_oi$mass_feature,
                   dat_type = "targeted")
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
    filter(!(row_id %in% to_drop))%>%
    left_join(new_names %>%
                  select(mass_feature, row_id, dat_type) %>%
                  rename(new_mass_feature = mass_feature,
                         new_dat_type = dat_type)) %>%
    mutate(mass_feature = ifelse(!is.na(new_mass_feature), new_mass_feature, mass_feature),
           dat_type = ifelse(!is.na(new_dat_type), new_dat_type, dat_type)) %>%
    select(-new_mass_feature, -new_dat_type)

# WRITE DATA -----
write_csv(dat_cmb3, here(output_loc, "combined_long_dat_scaled.csv"))
write_csv(dat_MF_update, here(output_loc, "combined_MF_info.csv"))
