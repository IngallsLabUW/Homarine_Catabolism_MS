library(RaMS)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))

# Prepare function and themes for plotting -----
my_chromat_theme <- theme_bw() +
    theme(
        plot.title = element_text(
            size = 7,
            margin=margin(b=-30),
            hjust = 0.5,
            face = "bold"
            ),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        legend.margin=margin(0,0,0,0),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        panel.border = element_blank()
            ) +
    # Make all margins 0
    theme(plot.margin = margin(0, 0, -50, 0))

## Function to plot super EICs ----
plot_combined_EICs <- function(mz_oi, rt_oi, name_oi, label_expected, mol_formula) {
  # Plot OBi1 and Rpom (no isotopes)
  g_ms1 <- plot_EIC2(
    ms1data = ms1data$MS1,
    m_z = mz_oi,
    mz_ppm = 10,
    r_t = rt_oi,
    rt_buffer = 2.3,
    samples_oi = sample_key
  ) +
    colScale +
    labs(title = paste0(name_oi)) +
      
  my_chromat_theme

  return(g_ms1)
}

# Get Data ------
## Set dirs ------
output_dir <- here("figures", "exploratory", "metabolomics", "combined_chromat_figs")
obi1_meta_data_dir <- here("data", "raw", "metabolomics", "obi1")
obi1_raw_dat_dir <- here("data", "raw", "metabolomics", "obi1", "particulate", "mzml")
rpom_metat_data_dir <- here("data", "raw", "metabolomics", "rpom")
rpom_raw_dat_dir <- here("data", "raw", "metabolomics", "rpom", "particulate", "mzml")

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


plus_h_smps <- c(obi1_h_samples, rpom_h_samples)
c_smps <- c(obi1_c_samples, rpom_c_samples)
smps_to_plot <- c(plus_h_smps, c_smps)

sample_key <- obi1_sample_key %>%
  bind_rows(rpom_sample_key) %>%
  mutate(treatment_short = case_when(
    filename %in% basename(plus_h_smps) ~ "homarine",
    filename %in% basename(c_smps) ~ "glucose",
    TRUE ~ "NA"
  )) %>%
  # cast treatment_short as factor
  mutate(treatment_short = factor(treatment_short, levels = c("glucose", "homarine"))) %>%
  filter(treatment_short != "NA") %>%
  mutate(experiment = factor(experiment, levels = c("obi1", "rpom"), labels = c("OBi1", "DSS-3"))) %>%
  distinct()

# Make plots ------
ms1data <- grabMSdata(files = smps_to_plot, grab_what = c("MS1"))
myColors <- c("#9467BD", "black")
names(myColors) <- levels(sample_key$treatment_short)
colScale <- scale_colour_manual(name = "Treatment", values = myColors, drop = FALSE)
### N-methyl glutamic acid figure ----
eic_acid <- plot_combined_EICs(
    mz_oi = 162.0759,
    rt_oi = 10,
    name_oi = "N-methlyglutamic acid",
    label_expected = T,
    mol_formula = "C6H11NO4") +
    theme(legend.position = "none",
          axis.title.x = element_blank())

### N-methyl glutamine figure ----
eic_amine <- plot_combined_EICs(
    mz_oi = 161.0921,
    rt_oi = 9.5,
    name_oi = "N-methylglutamine",
    label_expected = T,
    mol_formula = "C6H12N2O3") +
    theme(legend.position = "bottom")


### Combine plots ----
#g_chromat <- eic_acid / eic_amine +
#    # collect legends and axes titles
#    plot_layout(axis = "collect") +
#    plot_annotation(tag_levels = list(c("B", "C")))

