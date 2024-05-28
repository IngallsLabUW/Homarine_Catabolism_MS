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
        data_origin = "g4"
)
g5_dat <- read_csv(
    here(
        g5_dat_dir,
        "enrichment_summary.csv"
    ), show_col_types = FALSE
) %>%
    mutate(
        data_origin = "g5"
)

# Combine data -----
dat_cmb <- bind_rows(rpom_dat, obi1_dat) %>%
    bind_rows(g4_dat) %>%
    bind_rows(g5_dat)

# Pull out mass features of interest ------
mass_feature_oi <- c(
    #"Homarine",
    "N-Methyl-L-glutamic acid",
    "n-methyl glutamine",
    "n-methyl glutamine_D3",
    "n-methyl glutamine_13C_15N",
    "n-methyl glutamic acid",
    "n-methyl glutamic acid_D3",
    "n-methyl glutamic acid_13C_15N",
    "unknown_C6H9NO3_1_neg",
    "unknown_C6H9NO3_1_neg_D3",
    "unknown_C6H9NO3_1_neg_13C_15N",
    "unknown_C6H9NO3_2_neg",
    "unknown_C6H9NO3_2_neg_D3",
    "unknown_C6H9NO3_2_neg_13C_15N",
    "unknown_C6H9NO3_3_neg",
    "unknown_C6H9NO3_3_neg_D3",
    "unknown_C6H9NO3_3_neg_13C_15N",
    "unknown_C6H9NO4_neg",
    "unknown_C6H9NO4_neg_D3",
    "unknown_C6H9NO4_neg_13C_15N",
    "unknown_C6H9NO4_pos",
    "unknown_C6H9NO4_pos_D3",
    "unknown_C6H9NO4_pos_13C_15N",
    "unknown_C6H9NO4_pos_early",
    "unknown_C6H9NO4_pos_early_D3",
    "unknown_C6H9NO4_pos_early_13C_15N",
    "unknown_C7H9NO5_pos",
    "unknown_C7H9NO5_pos_D3",
    "unknown_C7H9NO5_pos_13C_15N"
    )

# Pull out particulate data and the core metabolites
dat_sub <- dat_cmb %>%
    filter((core_metabolite == "core" &
                sample_fraction == "Particulate") | (mass_feature %in% mass_feature_oi)) %>%
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
        select(mass_feature, core_metabolite, data_origin, sample_fraction, med_pvalue_h, med_pvalue_hNH4, experiment_id, timepoint),
    cols = contains("med_pvalue"),
    names_to = "treatment",
    values_to = "q_value",
    names_pattern = "med_pvalue_(.*)"
)

dat_long_plot <- dat_long %>%
    left_join(dat_long2, by = c("mass_feature", "core_metabolite", "data_origin", "sample_fraction", "treatment", "experiment_id", "timepoint")) %>%
    mutate(mass_feature = ifelse(mass_feature == "N-Methyl-L-glutamic acid", "n-methyl glutamic acid", mass_feature)) %>%
    mutate(p_value_flag = ifelse(q_value < 0.05, "significant", "not significant")) %>%
    mutate(enrichment = ifelse(enrichment > 2^10, 2^10, enrichment)) %>%
    filter(core_metabolite != "core") %>%
    bind_rows(dat_long_core)

# Make a tile plot, with y as mass_feature, x as treatment, and color as enrichment
my_fill_scale <- scale_fill_scico(palette = 'vik',
                                midpoint = 0,
                                limits = c(-3.5, 10))
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

g_cmb <- g_orgs + g4 + g5 +
    plot_layout(
        guides = "collect",
        axis_title = "collect",
        ncol = 3,
        widths = c(1, 1, 1.5)
    )
ggsave(here(output_dir, "intermediates_tiles.pdf"), g_cmb, width = 16, height = 4)

