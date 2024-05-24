# PURPOSE: To generate basic plots for the combined Obi1 Metabolomics data

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(here)
library(patchwork)
library(gghighlight)
library(ggrepel)

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
  select(mass_feature, sample_fraction)


# SUMMARY TABLES------
dat_enrich_summary <- dat_enrichment_results %>%
    inner_join(pseudos, by = c("mass_feature", "sample_fraction")) %>%
  group_by(mass_feature, sample_fraction) %>%
    summarise(
      med_pvalue_h = median(ttest_p_h),
      med_pvalue_hNH4 = median(ttest_p_hNH4),
      med_fc_h = median(enrich_h_fc),
      med_fc_hNH4 = median(enrich_hNH4_fc),
      med_pvalue_cmb = med_pvalue_h*med_pvalue_hNH4
    ) %>%
    ungroup() %>%
    mutate(
        med_qvalue_h = p.adjust(med_pvalue_h, method = "fdr"),
        med_qvalue_hNH4 = p.adjust(med_pvalue_hNH4, method = "fdr"),
        med_qvalue_cmb = p.adjust(med_pvalue_cmb, method = "fdr")
    ) %>%
    left_join(dat_MF, by = c("mass_feature", "sample_fraction")) %>%
    mutate(plot_label =
               case_when(
                   !str_detect(mass_feature, "HILIC") ~ mass_feature,
                   str_detect(mass_feature, "HILICPos") ~ paste0("\nRT = ", round(RT, digits = 1), "\n mz = ", round(mz, digits = 3), "+"),
                   str_detect(mass_feature, "HILICNeg") ~ paste0("\nRT = ", round(RT, digits = 1), "\n mz = ", round(mz, digits = 3), "-")
               ))


set.seed(123)
g4 <- ggplot(
    data = dat_enrich_summary,
    aes(
        x = log2(med_fc_h),
        y = log2(med_fc_hNH4),
        color = sample_fraction,
        shape = sample_fraction
    )
) +
    geom_point() +
    gghighlight(med_pvalue_h < 0.05 & log2(med_fc_h) > 0,
                use_direct_label = FALSE) +
    geom_text_repel(
        aes(label = plot_label),
        colour = "black",
        size = 2,
        )+
    theme_bw() +
    labs(
        x = "Median enrichment factor in homarine treatment (log2)",
        y = "Median enrichment factor in homarine + NH4 treatment (log2)",
    ) +
    scale_color_manual(values = c("aquamarine3", "darkorchid3", "black")) +
    scale_shape_manual(values = c(17, 16)) +
    theme(legend.position = c(0.1, 0.8),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "white", colour = "black"),
          text = element_text(size = 7))

plot(g4)

# SAVE OUTPUTS -----
ggsave(here(figure_dir, "EnrichmentFactors_scatter_psuedos.pdf"), plot = g4, height = 5, width = 6, units = "in")
