library(RaMS)
library(Rdisop)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
gc()
rm(list = ls())
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))

# Function to plot super EICs
plot_2H2_13C2 <- function(
        mz_oi,
        rt_oi,
        name_oi,
        ion_formula,
        ppm = 5
        ) {
    my_y_scale <- scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    my_theme <- theme(strip.text = element_blank())

    ## plot unlabeled compound ----
    mz_min = mz_oi - mz_oi * ppm / 1e6
    mz_max = mz_oi + mz_oi * ppm / 1e6
    g_ms1_og <- plot_EIC2(
        ms1data = ms1data$MS1 %>%
            filter(
                filename %in% c(
                    basename(g4_c_samples),
                    basename(g4_h_samples))
            ),
        m_z = mz_oi,
        mz_ppm = ppm,
        r_t = rt_oi,
        rt_buffer = 2.5,
        samples_oi = sample_key
    ) +
        annotate("text",
                 x=-Inf,
                 y=Inf,
                 hjust=0,
                 vjust=2,
                 label=paste0("Monoisotopic (m/z = " , round(mz_min, 5),
                              " - ", round(mz_max, 5), ")"))+
        colScale +
        my_y_scale +
        my_theme +
        labs(title = name_oi)
    g_ms1_og <- g_ms1_og + plot_spacer() + plot_spacer()

    ## +2 plots -----
    ### D2 labeled compound ----
    mz_2D_min = mz_oi + 1.006276746 * 2 - (mz_oi + 1.006276746 * 2) * ppm / 1e6
    mz_2D_max = mz_oi + 1.006276746 * 2 + (mz_oi + 1.006276746 * 2) * ppm / 1e6
    g_ms1_2H2 <- plot_EIC2(
      ms1data = ms1data$MS1 %>%
        filter(
            filename %in% c(basename(g4_c_samples), basename(g4_h_samples))
            ),
      m_z = mz_oi + 1.006276746 * 2, # 2 Ds heavier
      mz_ppm = ppm,
      r_t = rt_oi,
      rt_buffer = 2.5,
      samples_oi = sample_key
    ) +
      annotate("text",
               x=-Inf,
               y=Inf,
               hjust=0,
               vjust=2,
               label=paste0("D2 labeled (m/z = " , round(mz_2D_min, 5),
                            " - ", round(mz_2D_max, 5), ")"))+

        colScale +
        my_y_scale +
        my_theme

    ### 13C2 labeled compound----
    mz_13C2_min = mz_oi + 1.003355 * 2 - (mz_oi + 1.003355 * 2) * ppm / 1e6
    mz_13C2_max = mz_oi + 1.003355 * 2 + (mz_oi + 1.003355 * 2) * ppm / 1e6
    g_ms1_13C2 <- plot_EIC2(
        ms1data = ms1data$MS1 %>%
            filter(
                filename %in% c(basename(g4_c_samples), basename(g4_h_samples))
            ),
        m_z = mz_oi + 1.003355 * 2, # 2 Ds heavier
        mz_ppm = ppm,
        r_t = rt_oi,
        rt_buffer = 2.5,
        samples_oi = sample_key
    ) +
        annotate("text",
                 x=-Inf,
                 y=Inf,
                 hjust=0,
                 vjust=2,
                 label=paste0("13C2 labeled (m/z = " , round(mz_13C2_min, 5),
                              " - ", round(mz_13C2_max, 5), ")"))+
        colScale +
        my_y_scale +
        my_theme

    ### 13C2 simulated compound----
    # Get number of Cs from mol formula
    iso_dat <- t(getMolecule(ion_formula)$isotopes[[1]])
    mz_13C2 <- mz_oi + 1.003355 * 2
    # get the row of iso_dat that is closest to mz_13C2
    scaler = iso_dat[which.min(abs(iso_dat[,1] - mz_13C2)),2]/iso_dat[1,2]

    g_ms1_13C2_sim <- plot_EIC2(
        ms1data = ms1data$MS1 %>%
            filter(
                filename %in% c(basename(g4_c_samples), basename(g4_h_samples))
            ),
        m_z = mz_oi, # 2 Ds heavier
        mz_ppm = ppm,
        r_t = rt_oi,
        rt_buffer = 2.5,
        samples_oi = sample_key,
        scaler = scaler
    ) +
        annotate("text",
                 x=-Inf,
                 y=Inf,
                 hjust=0,
                 vjust=2,
                 label=paste0("13C2 simulated (m/z = " , round(mz_13C2_min, 5),
                              " - ", round(mz_13C2_max, 5), ")"))+
        colScale+
        my_y_scale +
        my_theme

  # Combine plots and fix the y-axis on the +2 plots
    g_plus2 <- g_ms1_2H2 + g_ms1_13C2 + g_ms1_13C2_sim +
        plot_layout(
            guides = "collect",
            axis_title = "collect",
            axes = "collect"
        )
    p_ranges_y <- c(ggplot_build(g_plus2[[1]])$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(g_plus2[[2]])$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(g_plus2[[3]])$layout$panel_scales_y[[1]]$range$range)
    g_plus2 <- g_plus2 &
        ylim(0, max(p_ranges_y)) &
        theme(legend.position = "none")

    ## +3 plots ----
    ### D3 labeled compound ----
    mz_3D_min = mz_oi + 1.006276746 * 3 - (mz_oi + 1.006276746 * 3) * ppm / 1e6
    mz_3D_max = mz_oi + 1.006276746 * 3 + (mz_oi + 1.006276746 * 3) * ppm / 1e6
    g_ms1_3D <- plot_EIC2(
        ms1data = ms1data$MS1 %>%
            filter(
                filename %in% c(basename(g4_c_samples), basename(g4_h_samples))
            ),
        m_z = mz_oi + 1.006276746 * 3, # 3 Ds heavier
        mz_ppm = ppm,
        r_t = rt_oi,
        rt_buffer = 2.5,
        samples_oi = sample_key
    ) +
        annotate("text",
                 x=-Inf,
                 y=Inf,
                 hjust=0,
                 vjust=2,
                 label=paste0("D3 labeled (m/z = " , round(mz_3D_min, 5),
                              " - ", round(mz_3D_max, 5), ")"))+
        colScale +
        my_y_scale +
        my_theme

    ### 13C3 labeled compound ----
    mz_13C3_min = mz_oi + 1.003355 * 3 - (mz_oi + 1.003355 * 3) * ppm / 1e6
    mz_13C3_max = mz_oi + 1.003355 * 3 + (mz_oi + 1.003355 * 3) * ppm / 1e6
    g_ms1_13C3 <- plot_EIC2(
        ms1data = ms1data$MS1 %>%
            filter(
                filename %in% c(basename(g4_c_samples), basename(g4_h_samples))
            ),
        m_z = mz_oi + 1.003355 * 3, # 3 Ds heavier
        mz_ppm = ppm,
        r_t = rt_oi,
        rt_buffer = 2.5,
        samples_oi = sample_key
    ) +
        annotate("text",
                 x=-Inf,
                 y=Inf,
                 hjust=0,
                 vjust=2,
                 label=paste0("13C3 labeled (m/z = " , round(mz_13C3_min, 5),
                              " - ", round(mz_13C3_max, 5), ")"))+
        colScale +
        my_y_scale +
        my_theme

    ### 13C3 simulated compound----
    # Get number of Cs from mol formula
    mz_13C3 <- mz_oi + 1.003355 * 3
    # get the row of iso_dat that is closest to mz_13C3
    scaler = iso_dat[which.min(abs(iso_dat[,1] - mz_13C3)),2]/iso_dat[1,2]
    g_ms1_13C3_sim <- plot_EIC2(
        ms1data = ms1data$MS1 %>%
            filter(
                filename %in% c(basename(g4_c_samples), basename(g4_h_samples))
            ),
        m_z = mz_oi, # 3 Ds heavier
        mz_ppm = ppm,
        r_t = rt_oi,
        rt_buffer = 2.5,
        samples_oi = sample_key,
        scaler = scaler
    ) +
        annotate("text",
                 x=-Inf,
                 y=Inf,
                 hjust=0,
                 vjust=2,
                 label=paste0("13C3 simulated (m/z = " , round(mz_13C3_min, 5),
                              " - ", round(mz_13C3_max, 5), ")"))+
        colScale+
        my_y_scale +
        my_theme

    # Combine plots and fix the y-axis on the +3 plots
    g_plus3 <- g_ms1_3D + g_ms1_13C3 + g_ms1_13C3_sim +
        plot_layout(
            guides = "collect",
            axis_title = "collect",
            axes = "collect"
        )
    p_ranges_y_3da <- c(ggplot_build(g_plus3[[1]])$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(g_plus3[[2]])$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(g_plus3[[3]])$layout$panel_scales_y[[1]]$range$range)
    g_plus3 <- g_plus3 &
        ylim(0, max(p_ranges_y_3da)) &
        theme(legend.position = "none")


    # Combine all plots----
  g_cmb <-
      g_ms1_og /
      g_plus2 /
      g_plus3 +
    plot_layout(
      axis_title = "collect",
      ncol = 1
    )
  #TODO KRH: Add title here
  return(g_cmb)
}

# GET FILES ------
## Set dirs ------
output_dir <- here("figures", "exploratory", "metabolomics", "combined_chromat_figs", "g4_chromats")
g4_meta_data_dir <- here("data", "raw", "metabolomics", "G4")
g4_raw_dat_dir <- here("data", "raw", "metabolomics", "G4", "particulate", "mzml_skyline")

## Read in sample keys and combine -------
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

## Set individual example samples ----
g4_h_samples <- c(
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t48_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t48_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_H_t48_C.mzML")
)
g4_c_samples <- c(
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C-H_t48_A.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C-H_t48_B.mzML"),
  here(g4_raw_dat_dir, "positive", "230724_Smp_G4_2_C-H_t48_C.mzML")
)

smps_to_plot <- c(g4_h_samples, g4_c_samples)

sample_key <- g4_sample_key %>%
  mutate(treatment_short = case_when(
    filename %in% basename(g4_h_samples) ~ "Plus Homarine",
    filename %in% basename(g4_c_samples) ~ "Control",
    TRUE ~ "NA"
  )) %>%
  # cast treatment_short as factor
  mutate(treatment_short = factor(treatment_short, levels = c("Control", "Plus Homarine"))) %>%
  filter(treatment_short != "NA") %>%
  distinct()

# PLOT POSITIVE MASS FEATURES Read in MS1 data ------
ms1data <- grabMSdata(files = smps_to_plot, grab_what = c("MS1"))

### Prep color scale to share between plots
myColors <- brewer.pal(3, "Set1")
names(myColors) <- levels(sample_key$treatment_short)
colScale <- scale_colour_manual(name = "Treatment", values = myColors, drop = FALSE)


### Homarine ------
g_h <- plot_2H2_13C2(
    mz_oi = 138.055503,
    rt_oi = 6.19,
    name_oi = "Homarine",
    ion_formula = "C7H8NO2")
ggsave(
    here(output_dir, "homarine.pdf"),
    g_h,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300)

### Glycine betaine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 118.08636474609375,
    rt_oi = 8,
    name_oi = "Glycine Betaine",
    ion_formula = "C5H12NO2")
ggsave(
    here(output_dir, "glycine_betaine.pdf"),
    g_cmb,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300)

### DMSP  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 135.047977,
    rt_oi = 9.325,
    name_oi = "DMSP",
    ion_formula = "C5H11O2S")

### Adenine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 136.06232,
    rt_oi = 4,
    name_oi = "Adenine",
    ion_formula = "C5H6N5")
ggsave(
    here(output_dir, "adenine.pdf"),
    g_cmb,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300)

### Adenosine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 268.10458,
    rt_oi = 5,
    name_oi = "Adenosine",
    ion_formula = "C10H14N5O4")
ggsave(
    here(output_dir, "adenosine.pdf"),
    g_cmb,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300)

### Hypoxanthine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 137.046336,
    rt_oi = 6.2,
    name_oi = "Hypoxanthine",
    ion_formula = "C5H5N4O")

### N6-Methyladenine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 150.07797,
    rt_oi = 4,
    name_oi = "N6-Methyladenine",
    ion_formula = "C6H8N5")

## Inosine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 269.088596,
    rt_oi = 7.8,
    name_oi = "Inosine",
    ion_formula = "C10H13N4O5")
ggsave(
    here(output_dir, "inosine.pdf"),
    g_cmb,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300)

### Xanthine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 153.041251,
    rt_oi = 7.5,
    name_oi = "Xanthine",
    ion_formula = "C5H5N4O2")

### Caffeine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 195.088201,
    rt_oi = 3.9,
    name_oi = "Caffeine",
    ion_formula = "C8H11N4O2")

### guanine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 152.057235,
    rt_oi = 8.1,
    name_oi = "Guanine",
    ion_formula = "C5H6N5O")

### guanosine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 284.099495,
    rt_oi = 9,
    name_oi = "Guanosine",
    ion_formula = "C10H14N5O5")


### Glycerophosphocholine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 258.110652,
    rt_oi = 10.86,
    name_oi = "Glycerophosphocholine",
    ion_formula = "C8H21NO6P")


### N-methyl-glutamic acid  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 162.076634,
    rt_oi = 10.5,
    name_oi = "N-methyl-glutamic acid",
    ion_formula = "C6H12NO4")
ggsave(
    here(output_dir, "n_methyl_glutamic_acid.pdf"),
    g_cmb,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300)

### Creatine  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 132.077302,
    rt_oi = 11.34,
    name_oi = "Creatine",
    ion_formula = "C4H10N3O2")

### Trimethylamine N-oxide  ------
g_cmb <- plot_2H2_13C2(
    mz_oi = 76.076239,
    rt_oi = 6.9,
    name_oi = "Trimethylamine N-oxide",
    ion_formula = "C3H9NO")



