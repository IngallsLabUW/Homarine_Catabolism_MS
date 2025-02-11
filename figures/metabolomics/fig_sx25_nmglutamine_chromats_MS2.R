library(RaMS)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))

# Gather data
data_loc <- here("data", "raw", "metabolomics", "obi1", "nmga_standard", "positive")
obi1_file <- here(data_loc, "230406_Smp_OBi1_4.mzML")
glutamind_std_file <- here("data", "raw", "metabolomics", "RC104", "particulate", "mzml", "HILIC_positive", "240926_Std_4uMStdsMix1InH2O-1_posDDA.mzml")

# Get data
msdata <- grabMSdata(files = c(obi1_file), grab_what = c("MS2"))
msdata_gluatmine <- grabMSdata(files = c(glutamind_std_file), grab_what = c("MS2"))

# MS2 data for glutamine from https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-BGC_Munich-RP001301
# ms2_data_glutamine <- data.frame(
#  fragmz = c(81.0449, 83.0606, 84.0442, 85.0269, 85.9987, 101.0712, 102.0547, 116.0205, 124.0215, 130.05, 132.0537, 147.076),
#  int = c(48, 52, 21494, 168, 42, 3440, 338, 50, 40, 13228, 62, 10940)
# )

# Make a data frame of filename, sample_name for renaming
sample_names <- data.frame(
  filename = c(
    "230406_Smp_OBi1_4.mzML",
    "240926_Std_4uMStdsMix1InH2O-1_posDDA.mzml"
  ),
  sample_name = c(
    "Cobetia sp. OBi1",
    "Glutamine Standard"
  )
)

# n-methyl-glutamine information
r_t <- 9.5
filename <- "230406_Smp_OBi1_4.mzML"
m_z <- 161.0921

# glutamine information
m_z_glutamine <- 147.0764
r_t_glutamine <- 11.7

mz_ppm <- 5
rt_buffer <- 5

# Plot EIC - matching run for NMGA to show time is same as standard

# Plot MS2
ms2_data <- msdata$MS2[
  premz %between% pmppm(m_z, mz_ppm) &
    rt %between% c(r_t - 0.5, r_t + 0.5)
] 
ms2_data_glutamine <- msdata_gluatmine$MS2[
  premz %between% pmppm(m_z_glutamine, mz_ppm) &
    rt %between% c(r_t_glutamine - 0.5, r_t_glutamine + 0.5)
] %>%
    filter(rt == unique(.$rt[1]))


# Normalize the intensity to max of 100
ms2_data <- ms2_data %>%
  group_by(filename) %>%
  mutate(int = int / max(int) * 100) %>%
  ungroup()
ms2_data_glutamine <- ms2_data_glutamine %>%
    group_by(filename) %>%
    mutate(int = int / max(int) * 100) %>%
    ungroup()

# Prepare single data frame for plotting
ms2_data <- ms2_data %>%
  mutate(treatment = "Observed") %>%
  bind_rows(ms2_data_glutamine %>%
    filter(fragmz < 147) %>%
    mutate(int = -int) %>%
    mutate(treatment = "Observed")) %>%
  bind_rows(ms2_data_glutamine %>%
    filter(fragmz < 147) %>%
    mutate(
      fragmz = fragmz + (m_z - 147.0764),
      int = -int
    ) %>%
    mutate(
      treatment = "MS2 fragments + 14"
    )) %>%
    left_join(sample_names, by = "filename")

# Prepare plot
my_theme <- theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 8)
  )
my_yscale <- scale_y_continuous(expand = c(0, NA))
my_xscale <- scale_x_continuous(limits = c(50, 165))

g_ms2_both <- ggplot() +
  geom_segment(
    data = ms2_data %>%
      filter(abs(int) > 1.5),
    aes(
      x = fragmz,
      xend = fragmz,
      yend = 0,
      y = int,
      color = treatment
    )
  ) +
  geom_segment(
    aes(x = m_z, xend = m_z, y = 0, yend = 100),
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = 147.076, xend = 147.0764, y = 0, yend = -100),
    linetype = "dashed"
  ) +
  # Add arrows to match peaks between pre- and post- modified masses
  geom_segment(
    data = ms2_data %>%
      filter(sample_name == "Glutamine Standard") %>%
      filter(treatment == "Observed") %>%
      filter(int < -1.5, fragmz < 147),
    aes(
      x = fragmz,
      xend = fragmz + (m_z - 147.0764),
      yend = int,
      y = int
    ),
    arrow = arrow(), size = 0.5, color = "grey"
  ) +
  geom_hline(yintercept = 0, color = "grey40") +
  annotate(
    "text",
    x = 50,
    y = -80,
    label = "glutamine standard",
    hjust = 0,
    size = 3
  ) +
    annotate(
        "text",
        x = 50,
        y = 80,
        label = "putative n-methyl glutamine \nin Cobetia sp. OBi1",
        hjust = 0,
        size = 3
    ) +
  my_theme +
  scale_color_manual(values = c("Observed" = "black", "MS2 fragments + 14" = "red")) +
  labs(y = "Normalized intensity", x = "m/z") +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.7, 0.9),
    legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size = 8))
  
g_ms2_both

# Save
ggsave(
  here(
    "figures", "metabolomics", "Figure_SX25_nmethylglutamine_identification.pdf"
  ), g_ms2_both,
  width = 6, height = 4, dpi = 300
)
