library(RaMS)
library(Rdisop)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))

# Function to plot super EICs
plot_combined_EICs <- function(mz_oi, rt_oi, name_oi, label_expected, mol_formula) {
  # Plot OBi1 and Rpom (no isotopes)
  g_ms1 <- plot_EIC2(
    ms1data = ms1data$MS1 %>%
      filter(!(filename %in% c(
        basename(g4_c_samples), basename(g4_h_samples),
        basename(g5_c_samples), basename(g5_h_samples)
      ))),
    m_z = mz_oi,
    mz_ppm = 10,
    r_t = rt_oi,
    rt_buffer = 2.5,
    samples_oi = sample_key
  ) +
    colScale +
    labs(title = paste0(name_oi, " (m/z ", mz_oi, ")"))

  # if labeled_expected is TRUE, plot the expected m/z of the D3 labeled compound and annotate
  if (label_expected) {
    ## plot D3 labeled compound
    g_ms1_D3 <- plot_EIC2(
      ms1data = ms1data$MS1 %>%
        filter(filename %in% c(basename(g4_c_samples), basename(g4_h_samples))),
      m_z = mz_oi + 1.006276746 * 3, # 3 Ds heavier
      mz_ppm = 10,
      r_t = rt_oi,
      rt_buffer = 2.5,
      samples_oi = sample_key
    ) +
        geom_text(
            data = . %>% filter(int == max(int)),
            parse = TRUE,
            size = 3,
            aes(
                x = rt_oi-1.5*.7,
                y = max(int)*.9,
                label = paste0("  italic(m/z)==", round((mz_oi + 1.006276746 * 3), 4))))+
        colScale +
        theme(legend.position = "none")
    # get the mz difference between the unlabeled and labeled compound
    Cs <- str_extract(mol_formula, "C\\d+") %>%
      str_extract("\\d+") %>%
      as.numeric()
    if (Cs > 7) {
      Cs <- 7
    }
    Ns <- ifelse(str_detect(mol_formula, "N"), 1, 0)
    mz_del <- 1.003355 * Cs + 0.997035 * Ns
    g_ms1_CNlab <- plot_EIC2(
      ms1data = ms1data$MS1 %>%
        filter(filename %in% c(basename(g5_c_samples), basename(g5_h_samples))),
      m_z = mz_oi + mz_del,
      mz_ppm = 10,
      r_t = rt_oi,
      rt_buffer = 2.5,
      samples_oi = sample_key
    ) +
        geom_text(
            data = . %>% filter(int == max(int)),
            parse = TRUE,
            size = 3,
            aes(
                x = rt_oi-1.5*.7,
                y = max(int)*.9,
                label = paste0("  italic(m/z)==", round((mz_oi + mz_del), 4))))+
        colScale +
      theme(legend.position = "none")
  } else {
    g_ms1_D3 <- plot_EIC2(
      ms1data = ms1data$MS1 %>%
        filter(filename %in% c(basename(g4_c_samples), basename(g4_h_samples))),
      m_z = mz_oi, # 3 Ds heavier
      mz_ppm = 10,
      r_t = rt_oi,
      rt_buffer = 2.5,
      samples_oi = sample_key
    ) + colScale +
      theme(legend.position = "none")
    g_ms1_CNlab <- plot_EIC2(
      ms1data = ms1data$MS1 %>%
        filter(filename %in% c(basename(g5_c_samples), basename(g5_h_samples))),
      m_z = mz_oi,
      mz_ppm = 10,
      r_t = rt_oi,
      rt_buffer = 2.5,
      samples_oi = sample_key
    ) + colScale +
      theme(legend.position = "none")
  }

  g_cmb <- g_ms1 + g_ms1_D3 + g_ms1_CNlab +
    plot_layout(
      guides = "collect",
      axis_title = "collect",
      ncol = 3,
      widths = c(2, 1, 1)
    )
  return(g_cmb)
}

# GET FILES ------
## Set dirs ------
output_dir <- here("figures", "exploratory", "metabolomics", "combined_chromat_figs")
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
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_C-H_t48_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_C-H_t48_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_1_C-H_t48_C.mzML")
)

g5_h_samples <- c(
  here(g5_raw_dat_dir, "HILIC_positive", "230717_Smp_G5_Hom2_T48_A.mzML"),
  here(g5_raw_dat_dir, "HILIC_positive", "230717_Smp_G5_Hom2_T48_B.mzML"),
  here(g5_raw_dat_dir, "HILIC_positive", "230717_Smp_G5_Hom2_T48_C.mzML")
)
g5_c_samples <- c(
  here(g5_raw_dat_dir, "HILIC_positive", "230717_Smp_G5_Hom2_C_A.mzML"),
  here(g5_raw_dat_dir, "HILIC_positive", "230717_Smp_G5_Hom2_C_B.mzML"),
  here(g5_raw_dat_dir, "HILIC_positive", "230717_Smp_G5_Hom2_C_C.mzML")
)


blank_smps <- c(obi1_b_samples, rpom_b_samples)
plus_h_smps <- c(obi1_h_samples, rpom_h_samples, g4_h_samples, g5_h_samples)
c_smps <- c(obi1_c_samples, rpom_c_samples, g4_c_samples, g5_c_samples)
smps_to_plot <- c(plus_h_smps, c_smps, blank_smps)

sample_key <- obi1_sample_key %>%
  bind_rows(rpom_sample_key) %>%
  bind_rows(g4_sample_key) %>%
  bind_rows(g5_sample_key) %>%
  mutate(treatment_short = case_when(
    filename %in% basename(plus_h_smps) ~ "Plus Homarine",
    filename %in% basename(c_smps) ~ "Control",
    filename %in% basename(blank_smps) ~ "Blank",
    TRUE ~ "NA"
  )) %>%
  # cast treatment_short as factor
  mutate(treatment_short = factor(treatment_short, levels = c("Control", "Plus Homarine", "Blank"))) %>%
  filter(treatment_short != "NA") %>%
  distinct()

# PLOT POSITIVE MASS FEATURES Read in MS1 data ------
ms1data <- grabMSdata(files = smps_to_plot, grab_what = c("MS1"))

### Prep color scale to share between plots
myColors <- brewer.pal(3, "Set1")
names(myColors) <- levels(sample_key$treatment_short)
colScale <- scale_colour_manual(name = "Treatment", values = myColors, drop = FALSE)

### Homarine figure ----
g_cmb <- plot_combined_EICs(
    mz_oi = 138.0549,
    rt_oi = 6.2,
    name_oi = "Homarine",
    label_expected = T,
    mol_formula = "C7H7NO2")
ggsave(here(output_dir, "homarine_eic.pdf"), g_cmb, width = 14, height = 4)

### Lysine figure ----
g_cmb <- plot_combined_EICs(
    mz_oi = 147.113353,
    rt_oi = 19,
    name_oi = "Lysine",
    label_expected = F,
    mol_formula = "C6H14N2O2")
ggsave(here(output_dir, "lysine_eic.pdf"), g_cmb, width = 14, height = 4)

## Serine figure ----
g_cmb <- plot_combined_EICs(
    mz_oi = 106.050419,
    rt_oi = 11.9,
    name_oi = "Serine",
    label_expected = F)
ggsave(here(output_dir, "serine_eic.pdf"), g_cmb, width = 14, height = 4)

### N-methyl glutamic acid figure ----
g_cmb <- plot_combined_EICs(
    mz_oi = 162.0759,
    rt_oi = 10.5,
    name_oi = "N-methyl glutamic acid",
    label_expected = T,
    mol_formula = "C6H11NO4")
ggsave(here(output_dir, "n_methyl_glutamic_acid_eic.pdf"), g_cmb, width = 14, height = 4)

### N-methyl glutamine figure ----
g_cmb <- plot_combined_EICs(
    mz_oi = 161.0921,
    rt_oi = 10,
    name_oi = "N-methyl glutamine",
    label_expected = T,
    mol_formula = "C6H12N2O3")
ggsave(here(output_dir, "n_methyl_glutamine_eic.pdf"), g_cmb, width = 14, height = 4)

### unknown_C6H9NO4_pos ----
g_cmb <- plot_combined_EICs(
    mz_oi = 160.061,
    rt_oi = 12,
    name_oi = "Unknown C6H9NO4 at 12 minutes",
    label_expected = T,
    mol_formula = "C6H9NO4")
ggsave(here(output_dir, "unknown_C6H9NO4_pos_12min_eic.pdf"), g_cmb, width = 14, height = 4)

### unknown_C6H9NO4_pos early----
g_cmb <- plot_combined_EICs(
    mz_oi = 160.061,
    rt_oi = 3.6,
    name_oi = "Unknown C6H9NO4 at 3.6 minutes",
    label_expected = T,
    mol_formula = "C6H9NO4")
ggsave(here(output_dir, "unknown_C6H9NO4_pos_3min_eic.pdf"), g_cmb, width = 14, height = 4)

### unknown_C7H9NO5_pos ----
g_cmb <- plot_combined_EICs(
    mz_oi = 188.055899,
    rt_oi = 12.5,
    name_oi = "Unknown C7H9NO5 at 12 minutes",
    label_expected = T,
    mol_formula = "C7H9NO5")
ggsave(here(output_dir, "unknown_C7H9NO5_pos_12min_eic.pdf"), g_cmb, width = 14, height = 4)

### unknown_C7H11NO3_pos ----
g_cmb <- plot_combined_EICs(
    mz_oi = 158.0817,
    rt_oi = 8.4,
    name_oi = "Unknown C7H11NO3 at 8 minutes",
    label_expected = T,
    mol_formula = "C7H11NO3")
ggsave(here(output_dir, "unknown_C7H11NO3_pos_8min_eic.pdf"), g_cmb, width = 14, height = 4)

# PLOT NEGATIVE MASS FEATURES Read in MS1 data ------
smps_to_plot_neg <- smps_to_plot %>%
  str_replace("positive", "negative")
ms1data <- grabMSdata(files = smps_to_plot_neg, grab_what = c("MS1"))

### unknown_C6H9NO4_neg ----
g_cmb <- plot_combined_EICs(
    mz_oi = 158.0453,
    rt_oi = 3.5,
    name_oi = "Unknown C6H9NO4 at 3.5 minutes",
    label_expected = T,
    mol_formula = "C6H9NO4")
ggsave(here(output_dir, "unknown_C6H9NO4_neg_3min_eic.pdf"), g_cmb, width = 14, height = 4)


### unknown_C6H9NO3_1_neg ----
g_cmb <- plot_combined_EICs(
    mz_oi = 142.0504,
    rt_oi = 11.5,
    name_oi = "Unknown C6H9NO3s around 11 minutes",
    label_expected = T,
    mol_formula = "C6H9NO3")
ggsave(here(output_dir, "unknown_C6H9NO3_11min_eic.pdf"), g_cmb, width = 14, height = 4)

