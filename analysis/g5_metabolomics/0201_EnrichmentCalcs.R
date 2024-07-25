# PURPOSE: calculate enrichment factors for the G5 samples

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(viridis)
library(broom)
library(coin)
library(here)

# SET FILE LOCATIONS  -----
meta_data_dir <- here("data", "intermediate", "metabolomics")
output_loc <- here("data", "intermediate", "metabolomics", "g5")

# Read in files  -----
dat_long_filename <- here(
  output_loc,
  "targeted",
  "combined_tidy_dat_long.csv"
)

# MUNGE DATA ----
## Combine data -----
dat_cmb <- read_csv(dat_long_filename, show_col_types = FALSE)

# If adjusted_area is NA, replace it with 1/50 of the minimum adjusted_area for that mass_feature
dat_cmb <- dat_cmb %>%
    group_by(mass_feature, sample_fraction) %>%
    # If all NA, drop from dataframe
    filter(!all(is.na(adjusted_area))) %>%
    mutate(adjusted_area = ifelse(is.na(adjusted_area), 0.02*min(adjusted_area, na.rm = TRUE), adjusted_area)) %>%
    ungroup() %>%
    mutate(replicate_name_short = str_remove(replicate_name, "^\\d\\d\\d\\d\\d\\d_")) %>%
    # drop problematic mass features
    filter(!(mass_feature %in% c(
        "L-Cysteic acid",
        "UDP-N-acetylglucosamine", "L-Homoserine")
        )) %>%
    # drop all _D2 mass features
    filter(!grepl("_D2", mass_feature))


## Pull out core metab ------
core_metabs <- read_csv(
    here(
    meta_data_dir,
    "coremetabs_working.csv"
), show_col_types = FALSE) %>%
    filter(core_metabolites_fMsystem %in% dat_cmb$mass_feature) %>%
    pull(core_metabolites_fMsystem)


## Calculate metabolite / core area for each compound x sample-----
# Keep cores and enrichments separate so we can test for trends in core metabolites
dat_enrich_list <- list()
dat_core_list <- list()
for (i in 1:length(core_metabs)) {
    dat_core_list[[i]] <- dat_cmb %>%
        filter(mass_feature == core_metabs[i]) %>%
        filter(sample_fraction == "Particulate") %>%
        mutate(core_metab = core_metabs[i]) %>%
        rename(core_adjusted_area = adjusted_area) %>%
        select(replicate_name_short, core_metab, core_adjusted_area)

    dat_enrich_list[[i]] <- dat_cmb %>%
        filter(mass_feature != core_metabs[i]) %>%
        left_join(dat_core_list[[i]], by = c("replicate_name_short")) %>%
        mutate(scaled_area = adjusted_area/core_adjusted_area)
}
dat_enrich <- bind_rows(dat_enrich_list) %>%
    filter(!is.na(core_metab))

# Pull out the Control samples
dat_enrich_control <- dat_enrich %>%
    filter(treatment == "Control") %>%
    select(-timepoint) %>%
    full_join(dat_enrich %>% select(timepoint, experiment_id) %>% distinct() %>%
                  filter(!is.na(timepoint)), by = "experiment_id",
              relationship = "many-to-many")
dat_enrich_control_fix <- dat_enrich %>%
    filter(treatment != "Control") %>%
    bind_rows(dat_enrich_control)


## Calculate per-sample ttest and mann-whitney test for each compound x core_metab combo ------
dat_enrich_h_tests <- dat_enrich_control_fix %>%
    select(mass_feature, core_metab, scaled_area, treatment, experiment_id, timepoint, sample_fraction) %>%
    filter(treatment %in% c("Control", "Homarine")) %>%
    group_by(mass_feature, core_metab, sample_fraction, experiment_id, timepoint) %>%
    summarise(ttest_p_h = t.test(scaled_area ~ treatment, paired = FALSE)$p.value,
              p_value_mannwhitney_h = wilcox.test(scaled_area ~ treatment)$p.value)

## Calculate fold changes for each compound x core_metab x timepoint x experiment combo  ------
dat_enrich3 <- dat_enrich_control_fix %>%
    filter(!is.na(treatment)) %>%
    group_by(mass_feature, treatment, experiment_id, timepoint, sample_fraction, core_metab) %>%
    summarise(mean_scaled_area = mean(scaled_area, rm.na = T),
              mean_adjusted_area = mean(adjusted_area, rm.na = T)) %>%
    ungroup() %>%
    group_by(mass_feature, core_metab, sample_fraction, experiment_id, timepoint) %>%
    reframe(enrich_h_fc = mean_scaled_area[treatment == "Homarine"] / mean_scaled_area[treatment == "Control"]) %>%
   # ungroup() %>%
    left_join(dat_enrich_h_tests, by = c("mass_feature", "core_metab", "sample_fraction", "experiment_id", "timepoint"))

## Test for outliers among core metabolites ------
core_outlier_tests <- list()
for (i in 1:length(unique(dat_enrich3$mass_feature))){
    # test if any of the core_metab values are outliers for the fold change
    core_outlier_tests[[i]] <- dat_enrich3 %>%
        filter(mass_feature == unique(dat_enrich3$mass_feature)[i]) %>%
        group_by(experiment_id, timepoint) %>%
        # Calculate the grubb's test for outliers
        mutate(outlier_h =
                   ifelse(
                       abs(enrich_h_fc - mean(enrich_h_fc)) > 2*sd(enrich_h_fc),
                       TRUE, FALSE)
        )

}
core_outlier_tests <- bind_rows(core_outlier_tests)
core_outliers <- core_outlier_tests %>%
    group_by(core_metab, experiment_id, timepoint) %>%
    summarise(outlier_f_h = sum(outlier_h,  na.rm = T)/n()) %>%
    ungroup() %>%
    group_by(core_metab) %>%
    summarise(outlier_f_h = mean(outlier_f_h, na.rm = T))

# pull core metabs in which less than 10% of the values are outliers
core_metabs_quality <- core_outliers %>%
    filter( outlier_f_h < 0.1) %>%
    filter(!is.na(core_metab))
write_csv(core_metabs_quality, here(output_loc, "core_metabs_quality.csv"))

## Remove core metabolites with high outlier rates ------
dat_enrich4 <- dat_enrich3 %>%
    filter(core_metab %in% core_metabs_quality$core_metab)
write_csv(dat_enrich4, here(output_loc, "enrichment_results.csv"))

## Summarize enrichment results ------
dat_enrich5 <- dat_enrich4 %>%
    group_by(mass_feature, sample_fraction, experiment_id, timepoint) %>%
    summarise(
        med_pvalue_h = median(ttest_p_h),
        med_fc_h = median(enrich_h_fc),
    ) %>%
    ungroup() %>%
    mutate(
        med_qvalue_h = p.adjust(med_pvalue_h, method = "fdr"),
    ) %>%
    mutate(core_metabolite = ifelse(mass_feature %in% core_metabs_quality$core_metab, "core", "non-core"))
write_csv(dat_enrich5, here(output_loc, "enrichment_summary.csv"))

