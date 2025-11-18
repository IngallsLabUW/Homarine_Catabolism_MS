# Load required libraries ----
library(here)
library(tidyverse)
library(ggplot2)
library(janitor)

fnt_size <- 5
# Load the data from Oscar of ncbi genomes with full pathway and assocaited taxonomy
dat <- read_delim(here("data", "intermediate", "ncbi_genome_counts", "complete_homarine_catabolic_operons_second_round_genome_metadata_and_taxonomy.tsv"),
  delim = "\t",
  show_col_types = FALSE
) %>%
  clean_names()

# Filter the data to only keep one nucleotide_accession per assembly
dat2 <- dat %>%
  group_by(assembly) %>%
  filter(nucleotide_accession == first(nucleotide_accession)) %>%
  ungroup()

# Clean up data for next steps
dat3 <- dat2 %>%
  select(
    assembly, nucleotide_accession, gene_name, start, end, gene_length,
    group_name, kingdom, phylum, class, order, family, genus, species
  ) %>%
  filter(!gene_name %in% c("BCCT"))


df2 <- dat3 %>%
  filter(!is.na(phylum)) %>%
  group_by(phylum, class, order) %>%
  summarise(value = n_distinct(assembly)) %>%
  ungroup() %>%
  mutate(class = case_when(
    phylum != "Pseudomonadota" ~ "Non-proteobacteria",
    TRUE ~ class
  ))

# Fill all NAs with "Unclassified"
df2[is.na(df2)] <- "Unclassified"

lvl0 <- tibble(name = "Parent", value = 0, level = 0, fill = "test")

lvl1 <- df2 %>%
  group_by(class) %>%
  summarise(value = sum(value)) %>%
  ungroup() %>%
  mutate(level = 1) %>%
  mutate(fill = class) %>%
  rename(name = class)

lvl2 <- df2 %>%
  group_by(class, order) %>%
  summarise(value = sum(value)) %>%
  ungroup() %>%
  mutate(fill = class) %>%
  select(name = order, value, fill) %>%
  mutate(level = 2)

final_df <- bind_rows(lvl1, lvl2) %>%
  mutate(name = as.factor(name) %>% fct_reorder2(fill, value)) %>%
  arrange(fill, name) %>%
  mutate(level = as.factor(level))

to_label <- final_df %>% filter(level == 2) %>%
    group_by(fill) %>%
    mutate(row_number = row_number()) %>%
    ungroup() %>%
    mutate(name = case_when(
    fill %in% c("Unclassified", "Non-proteobacteria", "Betaproteobacteria") ~ "",
      row_number < 6 & fill == "Gammaproteobacteria" ~ name,
      row_number < 4 & fill == "Alphaproteobacteria" ~ name,
      TRUE ~ ""
    ))

# TODO make points and labels separate than the columns?
sunburst <- ggplot(
    data = final_df,
    aes(x = level, y = value, fill = fill, alpha = level)
) +
    geom_col(
        width = 1, color = "black", size = 0.125, position = position_stack()
    ) +
    geom_text(
        data = to_label,
        aes(label = name), size = 1.5, position = position_stack(vjust = 0.5),
        alpha = 1
    ) +
    coord_polar(theta = "y") +
    scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 1), guide = FALSE) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    scale_fill_manual(values = c("#dd8fb9", "#6FA26F", "#16afca", "white", "grey")) +
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = fnt_size - 1),
        legend.key.size = unit(0.2, "cm"),  # Smaller key size
        legend.margin = margin(0, 0, 0, 0),  # Remove extra legend margin
        legend.spacing.y = unit(0.05, 'cm'),  # Reduce vertical spacing in legend
        legend.spacing.x = unit(0.05, 'cm'),  # Reduce horizontal spacing in legend
        plot.margin = unit(rep(-0.8, 4), "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.box.margin = margin(0, 0, 0, 0)
    )

