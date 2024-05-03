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
    "Homarine",
    "N-Methyl-L-glutamic acid",
    "n-methyl glutamine",
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
                sample_fraction == "Particulate") | (mass_feature %in% mass_feature_oi)) %>%
    select(mass_feature, core_metabolite, data_origin, sample_fraction, med_fc_h, med_fc_hNH4)

# pivot longer, use _h or _hNH4 as the suffix
dat_long <- pivot_longer(
    dat_sub,
    cols = contains("med_fc"),
    names_to = "treatment",
    values_to = "enrichment",
    names_pattern = "med_fc_(.*)"
)

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

