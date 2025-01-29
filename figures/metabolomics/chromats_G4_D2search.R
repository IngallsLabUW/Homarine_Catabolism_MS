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
        empirical_formula,
        ionization_form,
        ppm = 10
        ) {
    # format y scale to be scientific notation
    my_y_scale <- scale_y_continuous(expand = c(0, 0), limits = c(0, NA), labels = function(x) format(x, scientific = TRUE))
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


  # Combine plots and fix the y-axis on the +2 plots
    g_plus2 <- g_ms1_2H2 + g_ms1_13C2 + plot_spacer()+
        plot_layout(
            guides = "collect",
            axis_title = "collect",
            axes = "collect"
        )
    p_ranges_y <- c(ggplot_build(g_plus2[[1]])$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(g_plus2[[2]])$layout$panel_scales_y[[1]]$range$range)
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



    ### D2_13C1 labeled compound ----
    mz_2D13C1_min = (mz_oi + 1.003355 + 1.006276746 * 2) - (mz_oi + 1.003355 * 1.006276746 * 2) * ppm / 1e6
    mz_2D13C1_max = (mz_oi + 1.003355 + 1.006276746 * 2) + (mz_oi + 1.003355 * 1.006276746 * 2) * ppm / 1e6

    g_ms1_2D13C1 <- plot_EIC2(
        ms1data = ms1data$MS1 %>%
            filter(
                filename %in% c(basename(g4_c_samples), basename(g4_h_samples))
            ),
        m_z = mz_oi + 1.003355 + 1.006276746 * 2, # 2 Ds and 1 13C heavier
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
                 label=paste0("D2_13C1 labeled (m/z = " , round(mz_2D13C1_min, 5),
                              " - ", round(mz_2D13C1_max, 5), ")"))+
        colScale +
        my_y_scale +
        my_theme

    # Combine plots and fix the y-axis on the +3 plots
    g_plus3 <- g_ms1_3D + g_ms1_13C3 + g_ms1_2D13C1 +
        plot_layout(
            guides = "collect",
            axis_title = "collect",
            axes = "collect"
        )
    p_ranges_y_3da <- c(ggplot_build(g_plus3[[1]])$layout$panel_scales_y[[1]]$range$range,
                        ggplot_build(g_plus3[[2]])$layout$panel_scales_y[[1]]$range$range,
                        ggplot_build(g_plus3[[2]])$layout$panel_scales_y[[1]]$range$range)
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
input_dir <- here("data", "intermediate", "metabolomics", "g4")

## Read in data ------
g4_detected_features <- read_csv(
  here(
    input_dir, "grouped_iso_results_matched.csv"
  ),
  show_col_types = FALSE
)

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

# get unique compounds from g4_detected_features, grabbing name, mz, and rt
compounds_to_plot <- g4_detected_features %>%
  select(compound_name, mz_expected, rt_expected, mode, ionization_form, empirical_formula) %>%
  distinct()

# drop Trigonelline
compounds_to_plot <- compounds_to_plot %>%
  filter(compound_name != "Trigonelline")

# Note these are all positive so no need to get data for negative mode data
for (i in 1:nrow(compounds_to_plot)) {
  compound <- compounds_to_plot[i, ]
  # check if compound is positive, if so plot
  if (compound$mode == "positive") {
      g <- plot_2H2_13C2(
        mz_oi = compound$mz_expected,
        rt_oi = compound$rt_expected,
        name_oi = compound$compound_name,
        empirical_formula = compound$empirical_formula,
        ionization_form = compound$ionization_form)
      filename <- paste0(tolower(compound$compound_name), ".pdf") %>%
        str_replace_all(" ", "_")
      ggsave(
        here(output_dir, filename),
        g,
        width = 12,
        height = 10,
        units = "in",
        dpi = 300)
  }
}