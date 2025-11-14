# PURPOSE: make tile plots of the enrichment of core metabolites and conserved intermediates
# 20240528  KRH Making a tile plot of the enrichment of core metabolites and conserved intermediates
# TODO: species names should be R. pomeroyi and Cobetia sp. OBi1

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)
library(here)
library(patchwork)
library(ggtext)

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
    ),
    show_col_types = FALSE
) %>%
    mutate(
        data_origin = "rpom"
    )
obi1_dat <- read_csv(
    here(
        obi1_dat_dir,
        "enrichment_summary.csv"
    ),
    show_col_types = FALSE
) %>%
    mutate(
        data_origin = "obi1"
    )
g4_dat <- read_csv(
    here(
        g4_dat_dir,
        "enrichment_summary.csv"
    ),
    show_col_types = FALSE
) %>%
    mutate(
        data_origin = "g4",
        tag = "_D3"
    )
g5_dat <- read_csv(
    here(
        g5_dat_dir,
        "enrichment_summary.csv"
    ),
    show_col_types = FALSE
) %>%
    mutate(
        data_origin = "g5",
        tag = "_13C_15N"
    )
rc104_dat <- read_csv(
    here(
        rc104_dat_dir,
        "enrichment_summary.csv"
    ),
    show_col_types = FALSE
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

# Harmonize naming of metabolites without and with isotopes----
dat_cmb <- dat_cmb %>%
    mutate(
        mass_feature = ifelse(mass_feature == "N-Methyl-L-glutamic acid", "n-methyl glutamic acid", mass_feature),
        mass_feature = ifelse(mass_feature == "Homarine", "homarine", mass_feature),
        mass_feature = ifelse(mass_feature == "L-Glutamic acid", "glutamic_acid", mass_feature),
    ) %>%
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

# Summarize core metabolite enrichment signal ------
# Pull out core metabolites and mass features of interest
dat_sub <- dat_cmb %>%
    filter((core_metabolite == "core" &
                sample_fraction == "Particulate") | !is.na(mass_feature_oi)) %>%
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

# Summarize pvalues for all metabolites of interest ------
dat_long2 <- pivot_longer(
    dat_sub %>%
        select(
            mass_feature,
            mass_feature_short,
            mass_feature_oi, core_metabolite, data_origin, sample_fraction, med_pvalue_h, med_pvalue_hNH4, experiment_id, timepoint
        ),
    cols = contains("med_pvalue"),
    names_to = "treatment",
    values_to = "q_value",
    names_pattern = "med_pvalue_(.*)"
)

# Prepare data for plotting=====
dat_long_plot <- dat_long %>%
    left_join(
        dat_long2,
        by = c(
            "mass_feature", "core_metabolite", "data_origin",
            "sample_fraction", "treatment", "experiment_id", "timepoint"
        )
    ) %>%
    mutate(p_value_flag = ifelse(q_value < 0.05, "significant", "not significant")) %>%
    mutate(enrichment = ifelse(enrichment > 2^15, 2^15, enrichment)) %>%
    # Convert glutamic_acid to non-core for plotting
    mutate(
        core_metabolite = case_when(
            mass_feature == "glutamic_acid" & data_origin %in% c("rpom", "obi1") ~ "non-core",
            TRUE ~ core_metabolite
        )
    ) %>%
    filter(core_metabolite != "core") %>%
    bind_rows(dat_long_core) %>%
    # if mass_feature_short is na, fill in with mass_feature_oi
    mutate(
        mass_feature_short = case_when(
            str_detect(mass_feature, "core") ~ mass_feature,
            is.na(mass_feature_short) & !is.na(mass_feature_oi) ~ mass_feature_oi,
            TRUE ~ mass_feature_short
        )
    )

## Modify treatment and experiment_id columns for plotting =====
dat_long_plot <- dat_long_plot %>%
    mutate(
        experiment_id_to_plot = case_when(
            data_origin == "rpom" ~ "DSS-3",
            data_origin == "obi1" ~ "OBi1",
            data_origin == "rc104" ~ "RC104 St2",
            experiment_id == "G4_1" ~ "TN397 St2",
            experiment_id == "G4_2" ~ "TN397 St13",
            experiment_id == "G51" ~ "TN412 St2",
            experiment_id == "G52" ~ "TN412 St7"
        ) %>%
            factor(
                levels = c("OBi1", "DSS-3", "RC104 St2", "TN397 St2", "TN397 St13", "TN412 St2", "TN412 St7"),
                labels = c("OBi1", "DSS-3", "RC104 St2", "TN397 St2", "TN397 St13", "TN412 St2", "TN412 St7")
            ),
        mass_feature_oi = ifelse(is.na(mass_feature_oi), mass_feature_short, mass_feature_oi),
        mass_feature_oi = ifelse(mass_feature_oi == "", mass_feature_short, mass_feature_oi)
    ) %>%
    mutate(
        mass_feature_oi = ifelse(is.na(mass_feature_oi), mass_feature, mass_feature_oi)
    ) %>%
    filter(!(data_origin != "rpom" & data_origin != "obi1" & treatment == "hNH4")) %>%
    mutate(
        treatment = ifelse(
            data_origin %in% c("rpom", "obi1"),
            treatment, paste0(timepoint, "hr")
        ) %>%
            factor(
                levels = c(
                    "h", "hNH4", "0hr", "2hr", "6hr", "12hr", "24hr", "48hr", "96hr"
                ),
                labels = c(
                    "homarine",
                    "glucose + homarine",
                    "<1hr", "2hr", "6hr", "12hr", "24hr", "48hr", "96hr"
                )
            )
    )

# Drop unknown_C6H9NO2_3_neg and reorder the levels of mass_feature_oi
dat_long_plot2 <- dat_long_plot %>%
    filter(mass_feature_short != "unknown_C6H9NO3_3_neg") %>%
    mutate(
        mass_feature_short = factor(
            mass_feature_short,
            levels = c(
                "homarine",
                "n-methyl glutamic acid",
                "n-methyl glutamine",
                "glutamic_acid",
                "unknown_C6H9NO3_1_neg",
                "unknown_C6H9NO3_2_neg",
                "unknown_C6H9NO4_neg",
                "unknown_C6H9NO4_pos",
                "unknown_C7H9NO5_pos",
                "core_max",
                "core_median",
                "core_min"
            ),
            labels = c(
                "Homarine",
                "N-methylglutamic Acid",
                "N-methylglutamine",
                "Glutamic Acid",
                "C<sub>6</sub>H<sub>9</sub>NO<sub>3</sub> (-@rt10)",
                "C<sub>6</sub>H<sub>9</sub>NO<sub>3</sub> (-@rt11)",
                "C<sub>6</sub>H<sub>9</sub>NO<sub>4</sub> (-@rt3.5)",
                "C<sub>6</sub>H<sub>9</sub>NO<sub>4</sub> (+@rt12)",
                "C<sub>7</sub>H<sub>9</sub>NO<sub>5</sub> (+@rt12)",
                "Core metabs (max)",
                "Core metabs (med)",
                "Core metabs (min)"
            )
        )
    )

# PLOT -----
## Define themes, scales, etc to share between plots-------
# Make a tile plot, with y as mass_feature, x as treatment, and color as enrichment
# Dark blue #005F99 for low and #FF9896 for high, with white as 0
my_fill_scale <- scale_fill_gradient2(
    low = "#005F99",
    mid = "white",
    high = "#FF7876",
    midpoint = 0,
    limits = c(-3.5, 15),
    na.value = "grey80",
    breaks = c(-3, 0, 5, 10, 15)
)
my_theme <- theme_bw() +
    # remove x- and y-axis labels
    theme(
        axis.title.x = element_blank(),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8)
    ) +
    # Make all margins 0
    theme(plot.margin = margin(0, 0, 0, 0))

## Function to plot tile enrichment data -----
plot_tiles <- function(data_input) {
    ggplot(data = data_input) +
        geom_tile(
            aes(
                fill = log2(enrichment),
                y = mass_feature_short,
                x = treatment
            ),
            color = "black",
            width = 0.9,
            height = 0.9
        ) +
        geom_point(
            data = data_input %>%
                filter(p_value_flag == "significant"),
            aes(
                y = mass_feature_short,
                x = treatment,
            ),
            color = "black",
            shape = 8,
            size = 0.7
        ) +
        facet_wrap(
            facets = vars(experiment_id_to_plot)
        ) +
        coord_equal() +
        my_fill_scale +
        my_theme +
        # Add thick horizontal line above "Core metabs" group
        geom_hline(yintercept = 3.5, color = "black", linewidth = 1) +
        scale_y_discrete(limits = rev) +
        theme(axis.text.y = element_markdown(size = 7))
}

## Plot for RPom and Obi1 data  -----
dat_long_plot2_org <- dat_long_plot2 %>%
    filter(data_origin %in% c("rpom", "obi1")) 
levels(dat_long_plot2_org$mass_feature_short)[
    levels(dat_long_plot2_org$mass_feature_short) == "Homarine"
] <-
    "Homarine (C<sub>7</sub>H<sub>7</sub>NO<sub>3</sub>)"
levels(dat_long_plot2_org$mass_feature_short)[
    levels(dat_long_plot2_org$mass_feature_short) == "N-methylglutamic Acid"
] <-
    "N-methylglutamic Acid (C<sub>6</sub>H<sub>11</sub>NO<sub>4</sub>)"
levels(dat_long_plot2_org$mass_feature_short)[
    levels(dat_long_plot2_org$mass_feature_short) == "N-methylglutamine"
] <-
    "N-methylglutamine (C<sub>6</sub>H<sub>12</sub>N<sub>2</sub>O<sub>3</sub>)"
levels(dat_long_plot2_org$mass_feature_short)[
    levels(dat_long_plot2_org$mass_feature_short) == "Glutamic Acid"
] <-
    "Glutamic Acid (C<sub>5</sub>H<sub>9</sub>NO<sub>4</sub>)"

g_orgs <- plot_tiles(
    data_input = dat_long_plot2_org
) +
    # Fix the caption in the manuscript
    theme(
        legend.position = "bottom",  # Keeps the legend at the bottom
        legend.justification = "left",
        legend.direction = "horizontal", 
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 0, -60, unit = "pt"), # Adjusts the position
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        # rotate the x axis labels 45 degrees
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.key.size = unit(0.4, "cm"), # Reduces size of legend keys
        legend.spacing = unit(0.1, "cm") # Reduces spacing between legend items
    ) +
    labs(fill = "enrichment \n factor (log2)") +
    theme(axis.title.y = element_blank())

## Plot for G4 data -----
g4 <- plot_tiles(
    data_input = dat_long_plot2 %>%
        filter(data_origin == "g4")
) +
    theme(
        legend.position = "none",
        axis.title.y = element_markdown(
           # margin = margin(t = -50),
            size = 8
        )
    ) +
    ylab("<sup>2</sup>H<sub>3</sub> isotopologues")


## Plot for RC104 and G5 data -----
g5 <- plot_tiles(
    data_input = dat_long_plot2 %>%
        filter(data_origin %in% c("rc104", "g5"))) +
    theme(
      legend.position = "none",
        axis.title.y = element_markdown(
        #    margin = margin(t = -1500),
            size = 8
        )
  ) +
    # Add y axis label of "3H2 isotopologues" to the left of the plot
    ylab("<sup>13</sup>C<sub>5-7</sub>,<sup>15</sup>N isotopologues") 

## Map for homarine fate experiments -----
# this adds the object "hom_fate_map' to the environment
source(here("figures", "maps", "homarine_fate_experiments_map.R"))
hom_fate_map <- hom_fate_map

## Structures ----
structures_file <- list(
    here("figures", "structures", "homarine_v2.png"),
    here("figures", "structures", "nmethyl_glutamic_acid_isotopologue.png"),
    here("figures", "structures", "nmethyl_glutamine_isotopologue.png")
)
library(figpatch)
# Annotate titles directly onto the homarine figure
homarine <- fig(structures_file[[1]])
homarine <- fig_tag(homarine, "Homarine", fontsize = 6)

n_methyl_glutamatic_acid <- fig(structures_file[[2]]) 
n_methyl_glutamatic_acid <- fig_tag(n_methyl_glutamatic_acid, "N-methylglutamic acid", fontsize = 6)

n_methyl_glutamine <- fig(structures_file[[3]])
n_methyl_glutamine <- fig_tag(n_methyl_glutamine, "N-methylglutamine", fontsize = 6)

# Combine figures side-by-side in one row
structures <- wrap_plots(homarine,  n_methyl_glutamatic_acid, n_methyl_glutamine) +
    plot_layout(ncol = 3, widths = c(1,1,1)) # Combine everything into one row

# Add border around the entire combined plot
structured_with_border <- wrap_elements(plot = structures) +
    theme(
        plot.background = element_rect(
            color = "grey",  # Border color
            linewidth = 1       # Border thickness
        ),
        plot.margin = margin(2, 2, 2, 2, "pt")  # Add small margin to prevent border clipping
    )
structured_with_border

## Chromatogram for obi1 and rpom homarine control v +homarine -----
# this sources eic_acid and eic_amine
source(here("figures", "metabolomics", "chromats_subplot.R"))
g_chromat <- eic_acid / eic_amine +
    plot_layout(axis = "collect") +
    plot_annotation(tag_levels = list(c("B", "")))


## Combine plots -----
layout <- "
AABBBCCC
AABBBCCC
AABBBCCC
AABBBCCC
DEEEEEFF
DEEEEEFF
GGHHHHHH
GGHHHHHH
GGHHHHHH
"

g_cmb <- (free(g_orgs) +
              g_chromat +
              hom_fate_map +
              plot_spacer() + structured_with_border + plot_spacer() + 
              g4 +
              free(g5) +
              plot_layout(design = layout)) +
    plot_annotation(
        tag_levels = 
            list(c("     A", "     B", "     C", "     D", "     E", "     F", "     G"))) &
    theme(
        #plot.margin = unit(c(0, 0.1, 0, -0.6), "cm"),
        plot.margin = unit(c(0.05, 0.05, 0.05, -0.5), "cm"), #top, right, bottom, left
        
        panel.spacing = unit(0.3, "cm")
    )


g_cmb

#### Save plot -----
ggsave(
    here(
        "figures",
        "metabolomics", "homarine_fate_experiments.pdf"
    ),
    width = 8, height = 8, dpi = 300
)
ggsave(
    here(
        "figures",
        "metabolomics", "homarine_fate_experiments.png"
    ),
    width = 8, height = 8, dpi = 300
)
