# PURPOSE: make tile plots of the enrichment of core metabolites and conserved intermediates
# 20240528  KRH Making a tile plot of the enrichment of core metabolites and conserved intermediates

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(here)
library(patchwork)
library(scico)

# SET FILE LOCATIONS  -----
output_dir <- here("figures", "exploratory", "metabolomics")
rpom_data_dir <- here("data", "intermediate", "metabolomics", "rpom")
obi1_dat_dir <- here("data", "intermediate", "metabolomics", "obi1")
g4_dat_dir <- here("data", "intermediate", "metabolomics", "g4")
g5_dat_dir <- here("data", "intermediate", "metabolomics", "g5")
rc104_dat_dir <- here("data", "intermediate", "metabolomics", "rc104")

# Read in files  -----
rpom_dat <- read_csv(
    here(
        rpom_data_dir,
        "enrichment_summary.csv"
    ), show_col_types = FALSE
) %>%
    mutate(
        data_origin = "rpom"
    )
obi1_dat <- read_csv(
    here(
        obi1_dat_dir,
        "enrichment_summary.csv"
    ), show_col_types = FALSE
) %>%
    mutate(
        data_origin = "obi1"
)
g4_dat <- read_csv(
    here(
        g4_dat_dir,
        "enrichment_summary.csv"
    ), show_col_types = FALSE
) %>%
    mutate(
        data_origin = "g4",
        tag = "_D3"
)
g5_dat <- read_csv(
    here(
        g5_dat_dir,
        "enrichment_summary.csv"
    ), show_col_types = FALSE
) %>%
    mutate(
        data_origin = "g5",
        tag = "_13C_15N"
)
rc104_dat <- read_csv(
    here(
        rc104_dat_dir,
        "enrichment_summary.csv"
    ), show_col_types = FALSE
) %>%
    mutate(
        data_origin = "rc104",
        tag = "_13C_15N"
)


# Combine data -----
dat_cmb <- bind_rows(rpom_dat, obi1_dat) %>%
    bind_rows(g4_dat) %>%
    bind_rows(g5_dat) %>%
    bind_rows(rc104_dat)

# Add new column that is mass feature after removing the suffix tags
dat_cmb <- dat_cmb %>%
    mutate(mass_feature = ifelse(mass_feature == "N-Methyl-L-glutamic acid", "n-methyl glutamic acid", mass_feature),
           mass_feature = ifelse(mass_feature == "Homarine", "homarine", mass_feature)) %>%
    mutate(
        mass_feature_short = case_when(
            str_detect(mass_feature, tag) ~ str_remove(mass_feature, tag),
            TRUE ~ NA
        )
    ) %>%
    mutate(
        mass_feature_oi = case_when(
            !is.na(mass_feature_short) ~ mass_feature_short,
            tolower(mass_feature) %in% tolower(.$mass_feature_short) & is.na(tag) ~ mass_feature,
            TRUE ~ NA_character_
        )
    )

# Pull out core metabolites and mass features of interest
dat_sub <- dat_cmb %>%
    filter((core_metabolite == "core" &
                sample_fraction == "Particulate") | !is.na(mass_feature_oi))%>%
    filter(sample_fraction == "Particulate")

# pivot longer, use _h or _hNH4 as the suffix
dat_long <- pivot_longer(
    dat_sub %>%
        select(mass_feature, core_metabolite, data_origin, sample_fraction, med_fc_h, med_fc_hNH4, experiment_id, timepoint),
    cols = contains("med_fc"),
    names_to = "treatment",
    values_to = "enrichment",
    names_pattern = "med_fc_(.*)"
)

# Summarize the core metabolites' max and min values
dat_long_core <- dat_long %>%
    filter(core_metabolite == "core") %>%
    group_by(data_origin, sample_fraction, treatment, experiment_id, timepoint) %>%
    filter(!is.na(enrichment)) %>%
        summarize(
        core_median = median(enrichment, na.rm = TRUE),
        core_min = min(enrichment, na.rm = TRUE),
        core_max = max(enrichment, na.rm = TRUE)
    ) %>%
    pivot_longer(
        cols = c(core_median, core_min, core_max),
        names_to = "mass_feature",
        values_to = "enrichment"
    )

# Summarize pvalues
dat_long2 <- pivot_longer(
    dat_sub %>%
        select(mass_feature, mass_feature_oi, core_metabolite, data_origin, sample_fraction, med_pvalue_h, med_pvalue_hNH4, experiment_id, timepoint),
    cols = contains("med_pvalue"),
    names_to = "treatment",
    values_to = "q_value",
    names_pattern = "med_pvalue_(.*)"
)

dat_long_plot <- dat_long %>%
    left_join(dat_long2, by = c("mass_feature", "core_metabolite", "data_origin", "sample_fraction", "treatment", "experiment_id", "timepoint")) %>%
    mutate(p_value_flag = ifelse(q_value < 0.05, "significant", "not significant")) %>%
    mutate(enrichment = ifelse(enrichment > 2^15, 2^15, enrichment)) %>%
    filter(core_metabolite != "core") %>%
    bind_rows(dat_long_core)

# Modify treatment and experiment_id columns for plotting
dat_long_plot <- dat_long_plot %>%
    mutate(experiment_id_to_plot = case_when(
        data_origin == "rpom" ~ "RPom",
        data_origin == "obi1" ~ "Obi1",
        data_origin == "rc104" ~ "RC104",
        experiment_id == "G4_1" ~ "G4 (1)",
        experiment_id == "G4_2" ~ "G4 (2)",
        experiment_id == "G51" ~ "G5 (1)",
        experiment_id == "G52" ~ "G5 (2)"
    )) %>%
    mutate(mass_feature_oi = ifelse(is.na(mass_feature_oi), mass_feature, mass_feature_oi)) %>%
    filter(!(data_origin != "rpom" & data_origin != "obi1" & treatment == "hNH4")) %>%
    mutate(treatment = ifelse(
        data_origin %in% c("rpom", "obi1"), treatment, paste0(timepoint, "hr")
        ) %>% factor(levels = c("h", "hNH4", "0hr", "2hr", "6hr", "12hr", "24hr", "48hr", "96hr")))

# Drop unknown_C6H9NO2_3_neg
dat_long_plot <- dat_long_plot %>%
    filter(mass_feature_oi != "unknown_C6H9NO3_3_neg")



# Make a tile plot, with y as mass_feature, x as treatment, and color as enrichment
my_fill_scale <- scale_fill_scico(palette = 'vik',
                                midpoint = 0,
                                limits = c(-3.5, 15))
g_all <- ggplot() +
    geom_tile(
        data = dat_long_plot,
        aes(
            fill = log2(enrichment),
            y = mass_feature_oi,
            x = treatment
        ),
        color = "black",
        width = 0.9,
        height = 0.9) +
    # Add asterisks for significant p-values
    geom_text(
        data = dat_long_plot %>%
            filter(p_value_flag == "significant"),
        aes(
            y = mass_feature_oi,
            x = treatment,
            label = "*"
        ),
        color = "yellow",
        size = 4
    ) +
    facet_wrap(facets = vars(experiment_id_to_plot), scales = "free_x", ncol = 5) +
    my_fill_scale +
    theme_bw()
g_all


g_orgs <- ggplot() +
    geom_tile(
        data = dat_long_plot  %>%
            filter(data_origin %in% c("rpom", "obi1")),
        aes(
            fill = log2(enrichment),
            y = mass_feature,
            x = treatment
        ),
        color = "black",
        width = 0.9,
        height = 0.9) +
    facet_wrap(facets = vars(data_origin)) +
    my_fill_scale +
    theme_bw()

g4 <- ggplot() +
    geom_tile(
        data = dat_long_plot  %>%
            filter(data_origin %in% c("g4")) %>%
            filter(treatment == "h") %>%
            filter(!(mass_feature %in% c(
                "unknown_C6H9NO3_3_neg",
                "unknown_C6H9NO3_1_neg",
                "n-methyl glutamine",
                "n-methyl glutamic acid"
            ))),
        aes(
            fill = log2(enrichment),
            y = mass_feature,
            x = factor(timepoint)
        ),
        color = "black",
        width = 0.9,
        height = 0.9) +
    facet_wrap(facets = vars(experiment_id)) +
    my_fill_scale +
    theme_bw()

g5 <- ggplot() +
    geom_tile(
        data = dat_long_plot  %>%
            filter(data_origin %in% c("g5")) %>%
            filter(treatment == "h") %>%
            filter(!(mass_feature %in% c(
                "unknown_C6H9NO3_3_neg",
                "unknown_C6H9NO3_1_neg",
                "n-methyl glutamine",
                "n-methyl glutamic acid"
            ))),
        aes(
            fill = log2(enrichment),
            y = mass_feature,
            x = factor(timepoint)
        ),
        color = "black",
        width = 0.9,
        height = 0.9) +
    facet_wrap(facets = vars(experiment_id)) +
    my_fill_scale +
    theme_bw()

g_rc104 <- ggplot() +
    geom_tile(
        data = dat_long_plot  %>%
            filter(data_origin %in% c("rc104")) %>%
            filter(treatment == "h") %>%
            filter(!(mass_feature %in% c(
                "unknown_C6H9NO4_pos_early",
                "unknown_C6H9NO4_pos",
                "unknown_C6H9NO4_neg",
                "unknown_C6H9NO3_3_neg",
                "unknown_C6H9NO3_2_neg",
                "unknown_C6H9NO3_1_neg",
                "n-methyl glutamine",
                "n-methyl glutamic acid"
            ))),
        aes(
            fill = log2(enrichment),
            y = mass_feature,
            x = factor(timepoint)
        ),
        color = "black",
        width = 0.9,
        height = 0.9) +
    facet_wrap(facets = vars(experiment_id)) +
    my_fill_scale +
    theme_bw()

g_cmb <- g_orgs + g4 + g5 + g_rc104 +
    plot_layout(
        guides = "collect",
        axis_title = "collect",
        ncol = 4,
        widths = c(1, 1, 1, 1)
    )

g_cmb
ggsave(here(output_dir, "intermediates_tiles.pdf"), g_cmb, width = 16, height = 4)

