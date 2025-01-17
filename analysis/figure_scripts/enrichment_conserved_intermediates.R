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
    mass_feature = ifelse(mass_feature == "Homarine", "homarine", mass_feature)
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
  mutate(experiment_id_to_plot = case_when(
    data_origin == "rpom" ~ "RPom",
    data_origin == "obi1" ~ "Obi1",
    data_origin == "rc104" ~ "RC104",
    experiment_id == "G4_1" ~ "G4 (1)",
    experiment_id == "G4_2" ~ "G4 (2)",
    experiment_id == "G51" ~ "G5 (1)",
    experiment_id == "G52" ~ "G5 (2)"
  )) %>%
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
          "+homarine",
          "+homarine, \nNH4, and glucose",
          # TODO KRH: Check this
          "0hr", "2hr", "6hr", "12hr", "24hr", "48hr", "96hr"
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
        "n-Methyl Glutamic Acid",
        "n-Methyl Glutamine",
        # TODO KRH: double check - are these ion or molecular formula?
        # TODO KRH: fix the retention times
        "C6H9NO3 (-@rt = 3.5)",
        "C6H9NO3 (-@rt = 4.5)",
        "C6H9NO4 (-@rt = 5.5)",
        "C6H9NO4 (+@rt = 5.5)",
        "C7H9NO5 (+@rt = 6.5)",
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
  high = "#FF9896",
  midpoint = 0,
  limits = c(-3.5, 15),
  na.value = "grey80",
  breaks = c(-3, 0, 5, 10, 15)
)
my_theme <- theme_bw() +
  # remove x- and y-axis labels
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

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
    geom_text(
      data = data_input %>%
        filter(p_value_flag == "significant"),
      aes(
        y = mass_feature_short,
        x = treatment,
        label = "*"
      ),
      color = "black",
      size = 5
    ) +
    facet_wrap(
      facets = vars(experiment_id_to_plot)
    ) +
    my_fill_scale +
    my_theme +
    # Add thick horizontal line above "Core metabs" group
    geom_hline(yintercept = 3.5, color = "black", linewidth = 1) +
    scale_y_discrete(limits = rev) +
    theme(
      legend.position = "right"
    )
}

## Plot for RPom and Obi1 data  -----
g_orgs <- plot_tiles(
  data_input = dat_long_plot2 %>%
    filter(data_origin %in% c("rpom", "obi1"))
)

## Plot for G4 data -----
# TODO:Why do these start at 2 hr and the other start at 0?
g4 <- plot_tiles(
  data_input = dat_long_plot2 %>%
    filter(data_origin == "g4")
) +
  theme(legend.position = "none")

## Plot for RC104 and G5 data -----
g5 <- plot_tiles(
  data_input = dat_long_plot2 %>%
    filter(data_origin %in% c("rc104", "g5"))
) +
  theme(legend.position = "none")


## Combine plots -----
g_cmb <- g_orgs + g4 + g5 +
  plot_layout(
    guides = "collect",
    axis_title = "collect"
  ) +
  plot_annotation(tag_levels = 'A')

g_cmb



ggsave(here(output_dir, "intermediates_tiles.pdf"), g_cmb, width = 16, height = 4)
