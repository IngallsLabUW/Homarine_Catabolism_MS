#PURPOSE: To generate basic plots for the combined Obi1 Metabolomics data

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(RColorBrewer)
library(htmlwidgets)
library(plotly)
library(here)

# Set file locs ------
data_dir <- here("data", "intermediate", "metabolomics", "obi1")

# Read in files  -----
dat_long_filename <- here(data_dir, "combined_long_dat_scaled.csv")
dat_MF_filename <- here(data_dir, "combined_MF_info.csv")

dat_long <- read_csv(dat_long_filename, show_col_types = FALSE)
dat_MF <-
  read_csv(dat_MF_filename, show_col_types = FALSE)

# SUMMARY TABLES------
## Table of MF / mz / MS2 / FC h / FC hNH4 / ave h / ave hNH4 / ave control / pvalue h / pvalue hNH4
dat_tab <- dat_long %>%
  mutate(scaled_area = ifelse(scaled_area == 0, NA, scaled_area)) %>%
  group_by(mass_feature, treatment, sample_fraction) %>%
  # dump molecules that were not observed in all treatments
  summarise(
    mean_area = mean(scaled_area[!is.na(scaled_area)]),
    n_detected = n_distinct(scaled_area[!is.na(scaled_area)])
  ) %>%
  ungroup() %>%
  group_by(mass_feature, sample_fraction) %>%
  mutate(
    fc_hNH4 = mean_area[treatment == "Glucose + NH4 + Homarine"] /
      mean_area[treatment == "Glucose + NH4"],
    fc_h = mean_area[treatment == "Homarine"] / mean_area[treatment == "Glucose + NH4"],
    mean_overall = mean(mean_area, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(pivot_labels = case_when(
    treatment == "Glucose + NH4 + Homarine" ~ "hNH4",
    treatment == "Homarine" ~ "h",
    treatment == "Glucose + NH4" ~ "control",
  ))

dat_tab2 <- dat_tab %>%
  pivot_wider(
    id_cols = c(mass_feature, fc_hNH4, fc_h, mean_overall),
    names_from = pivot_labels,
    values_from = mean_area,
    names_glue = "ave_{pivot_labels}"
  ) %>%
  left_join(dat_tab %>%
    pivot_wider(
      id_cols = c(mass_feature),
      names_from = pivot_labels,
      values_from = n_detected,
      names_glue = "n_{pivot_labels}"
    )) %>%
  left_join(dat_MF %>%
    select(mass_feature, RT:add_annotation, high_homarine_compound:low_both_compounds))

# SUMMARY FIGS ----
dat_plot <- dat_tab2 %>%
  mutate(plot_label = ifelse(is.na(RT), mass_feature,
    paste0(mass_feature, ", RT = ", RT, ", mz = ", mz)
  )) %>%
  mutate(color_label = case_when(
    high_both_compound ~ "High in both homarine treatments",
    low_both_compounds ~ "Low in both homarine treatments",
    (high_homarine_compound | high_homarine_NH4_compound) & !(high_both_compound) ~ "High in single homarine treatment",
    (low_homarine_compound | low_homarine_NH4_compound) & !(low_both_compounds) ~ "Low in single homarine treatment",
    TRUE ~ "No change"
  ) %>%
    factor(
      levels = c(
        "High in both homarine treatments",
        "High in single homarine treatment",
        "No change",
        "Low in single homarine treatment",
        "Low in both homarine treatments"
      )
    )) %>%
  mutate(shape_label = case_when(
    str_detect(mass_feature, "HILIC") ~ "untargeted",
    TRUE ~ "targeted"
  ))



## Two way plots -----
g_h <- ggplot(
  data = dat_plot,
  aes(
    x = log10(mean_overall),
    y = log2(fc_h),
    fill = color_label,
    label = plot_label,
    shape = shape_label
  ),
  color = "black"
) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_brewer(type = "div", palette = "RdBu") +
  labs(
    title = "Particulate Samples only",
    x = "Log10(Mean peak area), across all treatments",
    y = "Log2(Fold change), of +homarine treatment"
  ) +
    theme_bw()

g_h_interactive <- ggplotly(g_h)

g_h2 <- ggplot(
    data = dat_plot,
    aes(
        x = RT,
        y = mz,
        fill = color_label,
        label = plot_label,
        shape = shape_label
    ),
    color = "black"
) +
    geom_point(size = 2) +
    scale_shape_manual(values = c(22, 21)) +
    scale_fill_brewer(type = "div", palette = "RdBu") +
    labs(
        title = "Particulate Samples only",
        x = "RT",
        y = "mz"
    ) +
    theme_bw()

g_h2_interactive <- ggplotly(g_h2)





# SAVE OUTPUTS -----
ggsave("04_combined/intermediates/FCvsMeanArea.pdf", plot = g_h, height = 6, width = 8, units = "in")
saveWidget(g_h_interactive, "04_combined/intermediates/FCvsMeanArea.html")
ic_write_csv_wlog(dat_tab2, "04_combined/intermediates/summary_MF_table.csv")
