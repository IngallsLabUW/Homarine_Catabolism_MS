library(RaMS)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))


# Gather data
data_loc <- here("data", "raw", "metabolomics", "obi1", "nmga_standard", "positive")
obi1_file <- here(data_loc, "230406_Smp_OBi1_4.mzML") # This file shows that the retention time of the NMGA standard is 10.2 minutes, matching the standard
nmga_standard <- here(data_loc, "230406_Std_5uM_Methylaminopentanedioic.mzML")
obi1_file_with_ms2 <- here(
  "data",
  "raw", "metabolomics", "obi1", "particulate", "mzml", "HILIC_positive_mzml", "211202_Poo_TruePooBacteria_DDApos50.mzml"
)

# Get data
msdata <- grabMSdata(files = c(obi1_file, nmga_standard, obi1_file_with_ms2), grab_what = c("MS1", "MS2"))

# Make a data frame of filename, sample_name for renaming
sample_names <- data.frame(
  filename = c(
      "230406_Smp_OBi1_4.mzML", 
      "230406_Std_5uM_Methylaminopentanedioic.mzML", 
      "211202_Poo_TruePooBacteria_DDApos50.mzml"
      ),
  sample_name = c(
      "Cobetia sp. OBi1", 
      "N-methylglutamic Acid Standard",
      "Cobetia sp. OBi1")
)

# n-methyl-glutamic acid information
r_t <- 10.2
filename <- "230406_Smp_OBi1_4.mzML"
rt_buffer <- 5
m_z <- 162.0759
mz_ppm <- 5

# Plot EIC - matching run for NMGA to show time is same as standard
rt_dat <- msdata$MS1 %>%
  select(rt, filename) %>%
  distinct() %>%
  filter(rt %between% c(r_t - rt_buffer, r_t + rt_buffer))

ms1_data <- msdata$MS1[mz %between% pmppm(m_z, mz_ppm) & rt %between% c(r_t - rt_buffer, r_t + rt_buffer)] %>%
  full_join(rt_dat, by = join_by(rt, filename)) %>%
  mutate(int = ifelse(is.na(int), 0, int)) %>%
    left_join(sample_names, by = "filename")

my_theme <- theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 8)) 
my_scale <- scale_y_continuous(expand = c(0, NA))

g_chromat_nmethyl_glutamic_acid <- 
    ggplot(ms1_data %>%
               # drop original run
               filter(filename != "211202_Poo_TruePooBacteria_DDApos50.mzml")) +
  geom_line(aes(x = rt, y = int, group = filename)) +
  facet_wrap(~sample_name, ncol = 1, scales = "free_y") +
    my_theme +
    my_scale+
  labs(
    y = "Intensity",
    x = "Retention Time (min)"
  )

# Plot MS2
ms2_data <- msdata$MS2[premz %between% pmppm(m_z, mz_ppm) & rt %between% c(r_t - 1, r_t + 1)]%>%
    left_join(sample_names, by = "filename")

g_ms2_nmethyl_glutamic_acid <- ggplot(ms2_data) +
  geom_segment(aes(x = fragmz, xend = fragmz, yend = 0, y = int, group = rt)) +
  facet_wrap(~sample_name, ncol = 1, scales = "free_y") +
    my_theme +
    
    my_scale+
    labs(
    y = "Intensity",
    x = "m/z"
  )

# Combine
g_cmb <- g_chromat_nmethyl_glutamic_acid + g_ms2_nmethyl_glutamic_acid +
    plot_layout(guides = "collect", axis = "collect")
g_cmb

# Save
ggsave(
    here(
        "figures", "metabolomics", "Figure_SX2_nmga_identification.pdf"), g_cmb, width = 6, height = 4, dpi = 300)
