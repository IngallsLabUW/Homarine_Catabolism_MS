library(RaMS)
library(Rdisop)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))


# GET FILES ------
## Set dirs ------
obi1_meta_data_dir <- here("data", "raw", "metabolomics", "obi1")
obi1_raw_dat_dir <- here("data", "raw", "metabolomics", "obi1", "particulate", "mzml")
rpom_metat_data_dir <- here("data", "raw", "metabolomics", "rpom")
rpom_raw_dat_dir <- here("data", "raw", "metabolomics", "rpom", "particulate", "mzml")
g4_meta_data_dir <- here("data", "raw", "metabolomics", "G4")
g4_raw_dat_dir <- here("data", "raw", "metabolomics", "G4", "particulate", "mzml_skyline")
g5_meta_data_dir <- here("data", "raw", "metabolomics", "G5")
g5_raw_dat_dir <- here("data", "raw", "metabolomics", "G5", "particulate", "mzml")

## Read in sample keys and combine -------
obi1_sample_key <- read_csv(
  here(
    obi1_meta_data_dir, "sample_key.csv"
  ),
  show_col_types = FALSE
) %>%
  rename(filename = replicate_name) %>%
  select(filename, treatment) %>%
  mutate(experiment = "obi1") %>%
  mutate(filename = paste0(filename, ".mzml"))
rpom_sample_key <- read_csv(
  here(
    rpom_metat_data_dir, "sample_key.csv"
  ),
  show_col_types = FALSE
) %>%
  rename(filename = replicate_name) %>%
  select(filename, treatment) %>%
  mutate(experiment = "rpom") %>%
  mutate(filename = paste0(filename, ".mzML"))
g4_sample_key <- read_csv(
  here(
    g4_meta_data_dir, "sample_key.csv"
  ),
  show_col_types = FALSE
) %>%
  rename(filename = replicate_name) %>%
  select(filename, treatment) %>%
  mutate(experiment = "g4") %>%
  mutate(filename = paste0(filename, ".mzML"))
g5_sample_key <- read_csv(
  here(
    g5_meta_data_dir, "sample_key.csv"
  ),
  show_col_types = FALSE
) %>%
  rename(filename = replicate_name) %>%
  select(filename, treatment) %>%
  mutate(experiment = "g5") %>%
  mutate(filename = paste0(filename, ".mzML"))

## Set individual example samples -----
obi1_h_samples <- c(
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Smp_OBi1_4.mzml"),
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Smp_OBi1_5.mzml"),
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Smp_OBi1_6.mzml")
)
obi1_c_samples <- c(
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Smp_OBi1_1.mzml"),
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Smp_OBi1_2.mzml"),
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Smp_OBi1_3.mzml")
)
obi1_b_samples <- c(
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Blk_BlankTube_1.mzml"),
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Blk_BlankTube_2.mzml"),
  here(obi1_raw_dat_dir, "HILIC_positive_mzml", "211202_Blk_BlankTube_3.mzml")
)

rpom_c_samples <- c(
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_WT_Gluc_A.mzML"),
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_WT_Gluc_B.mzML"),
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_WT_Gluc_C.mzML")
)
rpom_h_samples <- c(
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_WT_Hom_A.mzML"),
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_WT_Hom_B.mzML"),
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_WT_Hom_C.mzML")
)
rpom_b_samples <- c(
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_Ctrl_B3.mzML"),
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_Ctrl_A1.mzML"),
  here(rpom_raw_dat_dir, "HILIC_positive_mzml", "240304_Smp_Ctrl_B3.mzML")
)

g4_h_samples <- c(
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t48_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t48_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t48_C.mzML")
)
g4_c_samples <- c(
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C_t0_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C_t0_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C_t0_C.mzML")
)


blank_smps <- c(obi1_b_samples, rpom_b_samples)
plus_h_smps <- c(obi1_h_samples, rpom_h_samples, g4_h_samples)
c_smps <- c(obi1_c_samples, rpom_c_samples, g4_c_samples)
smps_to_plot <- c(plus_h_smps, c_smps, blank_smps)

sample_key <- obi1_sample_key %>%
  bind_rows(rpom_sample_key) %>%
  bind_rows(g4_sample_key) %>%
  mutate(treatment_short = case_when(
    filename %in% basename(plus_h_smps) ~ "Plus Homarine",
    filename %in% basename(c_smps) ~ "Control",
    filename %in% basename(blank_smps) ~ "Blank",
    TRUE ~ "NA"
  )) %>%
    distinct()

## POSITIVE MASS FEATURES Read in MS1 data ------
ms1data <- grabMSdata(files = smps_to_plot, grab_what = c("MS1"))

### Homarine figure ----
mz_oi <- 138.0549
rt_oi <- 6.2
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 100,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Homarine (m/z 138.0549)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "homarine_eic.pdf"), g_ms1, width = 8, height = 3)

## Lysine figure ----
mz_oi <- 147.113353
rt_oi <- 18.5
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Lysine (m/z 147.1128)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "lysine_eic.pdf"), g_ms1, width = 8, height = 3)

## Serine
mz_oi <- 106.050419
rt_oi <- 11.9
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Proline (m/z 116.0711)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "serine_eic.pdf"), g_ms1, width = 8, height = 3)


### N-methyl glutamic acid figure ----
mz_oi <- 162.0759
rt_oi <- 10
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "N-methyl glutamic acid (m/z 162.0759)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "n_methyl_glutamic_acid_eic.pdf"), g_ms1, width = 8, height = 3)


### N-methyl glutamine figure ----
mz_oi <- 161.0921
rt_oi <- 9.5
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "N-methyl glutamine (m/z 161.0921)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "n_methyl_glutamine_eic.pdf"), g_ms1, width = 8, height = 3)

### unknown_C6H9NO4_pos ----
mz_oi <- 160.061
rt_oi <- 12
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Unknown C6H9NO4 (m/z 160.061)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "unknown_C6H9NO4_pos_eic.pdf"), g_ms1, width = 8, height = 3)

rt_oi <- 3.6
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Unknown C6H9NO4 (m/z 160.061)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "unknown_C6H9NO4_pos_early_eic2.pdf"), g_ms1, width = 8, height = 3)


### unknown_C7H9NO5_pos ----
mz_oi <- 188.055899
rt_oi <- 12
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Unknown C7H9NO5 (m/z 188.055899)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "unknown_C7H9NO5_pos_eic.pdf"), g_ms1, width = 8, height = 3)

### unknown_C7H11NO3_pos ----
mz_oi <- 158.0817
rt_oi <- 8.4
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Unknown C7H11NO3 (m/z 158.0817)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "unknown_C7H11NO3_pos_eic.pdf"), g_ms1, width = 8, height = 3)



## NEGATIVE MASS FEATURES Read in MS1 data ------
smps_to_plot_neg <- smps_to_plot %>%
  str_replace("positive", "negative")
ms1data <- grabMSdata(files = smps_to_plot_neg, grab_what = c("MS1"))

### unknown_C6H9NO4_neg ----
mz_oi <- 158.0453
rt_oi <- 3.5
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Unknown C6H9NO4 (m/z 160.061)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "unknown_C6H9NO4_neg_eic.pdf"), g_ms1, width = 8, height = 3)

### unknown_C6H9NO3_1_neg ----
mz_oi <- 142.0504
rt_oi <- 11
g_ms1 <- plot_EIC2(
  ms1data = ms1data$MS1,
  m_z = mz_oi,
  mz_ppm = 5,
  r_t = rt_oi,
  rt_buffer = 2.5,
  samples_oi = sample_key
) +
  labs(title = "Unknown C6H9NO3 (m/z 142.0504)", x = "Retention Time (min)", y = "Intensity")
ggsave(here("figures", "exploratory", "metabolomics", "unknown_C6H9NO3_all_neg_eic.pdf"), g_ms1, width = 8, height = 3)
