# PURPOSE: combine and tidy the untargeted data with the targeted data
# 20240313  KRH moving from Obi1-specific repo

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(janitor)
library(viridis)
library(broom)
library(coin)
library(here)
library(gghighlight)

#TODO:
# Figure out where rpom particulate homarine data are?

# SET FILE LOCATIONS  -----
rpom_data_dir <- here("data", "intermediate", "metabolomics", "rpom")
obi1_dat_dir <- here("data", "intermediate", "metabolomics", "obi1")

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

# Combine data -----
dat_cmb <- bind_rows(rpom_dat, obi1_dat)

# Pull out mass features of interest ------
mass_feature_oi <- c(
    #"Homarine",
    "N-Methyl-L-glutamic acid",
    "n-methyl glutamine",
    "n-methyl glutamic acid",
    "unknown_C6H9NO3_1_neg",
    "unknown_C6H9NO3_2_neg",
    "unknown_C6H9NO3_3_neg",
    "unknown_C6H9NO4_neg",
    "unknown_C6H9NO4_pos",
    "unknown_C6H9NO4_pos_early",
    "unknown_C7H9NO5_pos"
    )

# Plot Homarine and core metabs' enrichments as an example
dat_sub <- dat_cmb %>%
    filter((core_metabolite == "core" &
                sample_fraction == "Particulate") | (mass_feature %in% mass_feature_oi))


# pivot longer, use _h or _hNH4 as the suffix
dat_long <- pivot_longer(
    dat_sub %>%
        select(mass_feature, core_metabolite, data_origin, sample_fraction, med_fc_h, med_fc_hNH4),
    cols = contains("med_fc"),
    names_to = "treatment",
    values_to = "enrichment",
    names_pattern = "med_fc_(.*)"
)

dat_long2 <- pivot_longer(
    dat_sub %>%
        select(mass_feature, core_metabolite, data_origin, sample_fraction, med_pvalue_h, med_pvalue_hNH4),
    cols = contains("med_pvalue"),
    names_to = "treatment",
    values_to = "q_value",
    names_pattern = "med_pvalue_(.*)"
)
dat_long <- left_join(dat_long, dat_long2, by = c("mass_feature", "core_metabolite", "data_origin", "sample_fraction", "treatment")) %>%
    # Rename "N-Methyl-L-glutamic acid" to "n-methyl glutamic acid"
    mutate(mass_feature = ifelse(mass_feature == "N-Methyl-L-glutamic acid", "n-methyl glutamic acid", mass_feature)) %>%
    mutate(p_value_flag = ifelse(q_value < 0.05, "significant", "not significant"))

g <- ggplot() +
    geom_point(
        data = dat_long %>%
            filter(core_metabolite == "core"),
        aes(
            shape = sample_fraction,
            y = log2(enrichment),
            x = treatment
        ),
        color = "grey",
        alpha = 0.5,
        width = 0.1,
        position = position_dodge2(width = 0.1, preserve = "single")) +
    geom_point(
        data = dat_long %>%
            filter(core_metabolite == "non-core"),
        aes(
            shape = sample_fraction,
            color = mass_feature,
            y = log2(enrichment),
            x = treatment
        ),
        size = 3,
        position = position_dodge2(width = 0.2, preserve = "single")) +
    facet_wrap(facets = vars(data_origin)) +
    # use good colors for discrete variables
    scale_color_brewer(palette = "Set1")+
    theme_bw()
g

# Make a tile plot, with y as mass_feature, x as treatment, and color as enrichment
# use a color palette that is centered at 0 (grey for 0)
g2 <- ggplot() +
    geom_tile(
        data = dat_long ,
        aes(
            fill = log2(enrichment),
            y = mass_feature,
            x = interaction(treatment, sample_fraction)
        ),
        color = "black",
        width = 0.9,
        height = 0.9) +
    facet_wrap(facets = vars(data_origin)) +
    scale_fill_viridis_c(option = "magma", na.value = "grey") +
    theme_bw()
g2
