# PURPOSE: combine and tidy the untargeted data with the targeted data
# 20240313  KRH moving from Obi1-specific repo

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(viridis)
library(broom)
library(coin)
library(here)

# SET FILE LOCATIONS  -----
meta_data_dir <- here("data", "intermediate", "metabolomics")
output_loc <- here("data", "intermediate", "metabolomics", "obi1")

# Read in files  -----
dat_target_filename <- here(
  output_loc,
  "targeted",
  "combined_tidy_dat_long.csv"
)
dat_nontarget_filename <- here(
  output_loc,
  "untargeted",
  "combined_tidy_dat_long.csv"
)

dat_target <- read_csv(dat_target_filename, show_col_types = FALSE) %>%
  rename(RT = retention_time)
dat_nontarget <-
  read_csv(dat_nontarget_filename, show_col_types = FALSE) %>%
  mutate(z = case_when(
    grepl("Pos", mass_feature) ~ 1,
    grepl("Neg", mass_feature) ~ -1
  ))

# MUNGE DATA ----
## Combine data -----
dat_cmb <- dat_target %>%
  mutate(dat_type = "targeted") %>%
  bind_rows(dat_nontarget %>%
    mutate(dat_type = "nontargeted")) %>%
  filter(sample_set == "OBi1_set2")

# If adjusted_area is NA, replace it with 1/50 of the minimum adjusted_area for that mass_feature
dat_cmb <- dat_cmb %>%
    group_by(mass_feature, sample_fraction) %>%
    # If all NA, drop from dataframe
    filter(!all(is.na(adjusted_area))) %>%
    mutate(adjusted_area = ifelse(is.na(adjusted_area), 0.02*min(adjusted_area, na.rm = TRUE), adjusted_area)) %>%
    ungroup()


## Pull out core metab ------
core_metabs <- read_csv(
    here(
    meta_data_dir,
    "coremetabs_working.csv"
), show_col_types = FALSE) %>%
    filter(core_metabolites_fMsystem %in% dat_cmb$mass_feature) %>%
    # Drop citrulline
#    filter(core_metabolites_fMsystem != "Citrulline") %>%
    pull(core_metabolites_fMsystem)


dat_enrich_list <- list()
dat_core_list <- list()
for (i in 1:length(core_metabs)) {
    dat_core_list[[i]] <- dat_cmb %>%
        filter(mass_feature == core_metabs[i]) %>%
        mutate(core_metab = core_metabs[i]) %>%
        rename(core_adjusted_area = adjusted_area) %>%
        select(replicate_name, sample_fraction,core_metab, core_adjusted_area)

    dat_enrich_list[[i]] <- dat_cmb %>%
        filter(mass_feature != core_metabs[i]) %>%
        left_join(dat_core_list[[i]], by = c("replicate_name", "sample_fraction")) %>%
        mutate(scaled_area = adjusted_area/core_adjusted_area)
}

dat_enrich <- bind_rows(dat_enrich_list)
dat_core <- bind_rows(dat_core_list)


## Calculate fold changes and select median fold change for scaling ------
dat_enrich2 <- dat_enrich %>%
    group_by(mass_feature, treatment, sample_fraction, core_metab) %>%
    summarise(mean_scaled_area = mean(scaled_area, rm.na = T)) %>%
    ungroup()

dat_enrich3 <- dat_enrich2 %>%
    group_by(mass_feature, sample_fraction, core_metab) %>%
    summarise(enrich_hNH4_fc = mean_scaled_area[treatment == "Glucose + NH4 + Homarine"] /  mean_scaled_area[treatment == "Glucose + NH4"],
              enrich_h = mean_scaled_area[treatment == "Homarine"] / mean_scaled_area[treatment == "Glucose + NH4"]) %>%
    pivot_longer(cols = c("enrich_hNH4_fc", "enrich_h"), names_to = "treatment", values_to = "enrichment")

# get order of mass_features by median enrichment
dat_enrich3 %>%
    filter(sample_fraction == "Particulate") %>%
    group_by(mass_feature) %>%
    summarise(median_enrichment = median(enrichment, na.rm = T)) %>%
    arrange(median_enrichment) %>%
    pull(mass_feature) %>%
    unique() -> mass_feature_order
# set order of mass_features
dat_enrich3$mass_feature <- factor(dat_enrich3$mass_feature, levels = mass_feature_order)

dat_enrich_summary <- dat_enrich3 %>%
    group_by(mass_feature, sample_fraction, treatment) %>%
    summarise(median_enrichment = median(enrichment, na.rm = T),
              mean_enrichment = mean(enrichment, na.rm = T),
              sd_enrichment = sd(enrichment, na.rm = T),
              n = sum(!is.na(enrichment))) %>%
    ungroup()

# Plot boxplots of enrichment ------
g <- ggplot(
    dat_enrich_summary %>%
        filter(sample_fraction == "Particulate") %>%
        filter(!(mass_feature %>% str_detect("HILIC"))),
    aes(
        x = mass_feature,
        y = log2(median_enrichment),
        color = treatment,
        )) +
    geom_segment(
        aes( y = log2(median_enrichment- sd_enrichment),
            yend = log2(median_enrichment + sd_enrichment))
    ) +
    geom_point(size = 1) +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Median enrichment compared to control (log2)", x = "Mass feature")
g

c2 <geom_point()c2 <- ggplot(data = dat3cul, aes(factor(Identification), intracell_conc_umolCL/1000)) +
    geom_boxplot(outlier.size = 0.3, lwd = 0.3) +
    scale_y_log10(breaks = c(0.00001,  0.001,  0.1, 10, 1000), labels = function(x) sprintf("%g", x))+
    coord_flip() +
    theme(axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text (size = 5),
          axis.title.x= element_text (size = 6, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          panel.grid.major.y = element_line(colour="lightgrey", size = rel(0.2)),
          plot.margin = unit(c(0,0,0,0), "mm"))+
    labs(y = "mmol C / L \nintracellcular concentrations")
