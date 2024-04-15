# PURPOSE: To generate basic plots for the combined Obi1 Metabolomics data

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(RColorBrewer)
library(htmlwidgets)
library(plotly)
library(here)
library(patchwork)

# Set file locs ------
data_dir <- here("data", "intermediate", "metabolomics", "obi1")
figure_dir <- here("figures", "exploratory", "metabolomics", "obi1")

# Read in files  -----
dat_MF_filename <- here(data_dir, "combined_MF_info.csv") #PARTICULATE ONLY
dat_MF_enrichment_filename <- here(data_dir, "enrichment_results.csv")

dat_MF <-
  read_csv(dat_MF_filename, show_col_types = FALSE)
dat_enrichment_results <-
  read_csv(dat_MF_enrichment_filename, show_col_types = FALSE)

# GET LIST OF PSEUDOS OR TARGETED ------
pseudos <- dat_MF %>%
  filter(!str_detect(mass_feature, "HILIC") | !is.na(add_annotation)) %>%
  pull(mass_feature)


# SUMMARY TABLES------
dat_enrich_summary <- dat_enrichment_results %>%
    filter(mass_feature %in% pseudos) %>%
    filter(sample_fraction == "Particulate") %>%
  group_by(mass_feature, sample_fraction) %>%
    summarise(
      pvalue_h = median(ttest_p_h),
      pvalue_hNH4 = median(ttest_p_hNH4),
      ave_h = median(enrich_h_fc),
      ave_hNH4 = median(enrich_hNH4_fc)
    ) %>%
    left_join(dat_MF, by = c("mass_feature", "sample_fraction")) %>%
    mutate(plot_label = ifelse(!str_detect(mass_feature, "HILIC"), mass_feature,
                               paste0(mass_feature, "\n RT = ", RT, ", mz = ", mz)))

library(gghighlight)
library(ggrepel)
g4 <- ggplot(
    data = dat_enrich_summary,
    aes(
        x = log2(ave_h),
        y = log2(ave_hNH4),
        text = plot_label
    )
) +
    geom_point(color = "red") +
    gghighlight(pvalue_h < 0.05 & log2(ave_h) > 0,
                use_direct_label = FALSE) +
    geom_text_repel(aes(label = plot_label), colour = "black", size = 2.5)+
    theme_bw() +
    labs(
        x = "Median enrichment factor in homarine treatment (log2)",
        y = "Median enrichment factor in homarine + NH4 treatment (log2)",
    )

g4



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

# Subset to only mass features that we can get pseudomolecular formula for
dat_plot_psuedos_only <- dat_plot %>%
  filter(!str_detect(mass_feature, "HILIC") | !is.na(add_annotation))

## Two way plots -----
g_h <- ggplot(
  data = dat_plot_psuedos_only,
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


# SAVE OUTPUTS -----
ggsave(here(figure_dir, "FCvsMeanArea_particulate_pseudos.pdf"), plot = g_h, height = 6, width = 8, units = "in")
saveWidget(g_h_interactive, here(figure_dir, "FCvsMeanArea_particulate_pseudos.html"))
write_csv(dat_tab2, here(data_dir, "summary_MF_table.csv"))
