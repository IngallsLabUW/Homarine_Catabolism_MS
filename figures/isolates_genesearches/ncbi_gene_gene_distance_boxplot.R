# Load required libraries ----
library(here)
library(tidyverse)
library(ggplot2)
library(janitor)
library(patchwork)

closest_distance <- function (start_x, end_x, start_y, end_y) {

    # Calculate the distance between the start and end of x and y
    distance_start = abs(start_x - start_y)
    distance_end = abs(end_x - start_y)
    distance_start_end = abs(start_x - end_y)
    distance_end_start = abs(end_x - start_y)
    
    # Return the minimum distance
    return(min(distance_start, distance_end, distance_start_end, distance_end_start))
}

# Load the data from Oscar of ncbi genomes with full pathway and assocaited taxonomy
dat <- read_delim(here("data", "intermediate", "ncbi_genome_counts", "complete_homarine_catabolic_operons_second_round_genome_metadata_and_taxonomy.tsv"), delim = "\t",
                         show_col_types = FALSE) %>%
    clean_names()
# Load data with genome size
dat_genome_size <- read_csv(here("analysis/ncbi_genomes/genome_sizes_datasets.csv"),
                         show_col_types = FALSE) %>%
    clean_names() %>%
    select(assembly_accession, genome_size_bp) %>%
    rename(assembly = assembly_accession) %>%
    distinct()

# Filter the data to only keep one nucleotide_accession per assembly
dat2 <- dat %>%
    group_by(assembly) %>%
    filter(nucleotide_accession == first(nucleotide_accession)) %>%
    ungroup() 

# Clean up data for next steps
dat3 <- dat2 %>%
    select(assembly, nucleotide_accession, gene_name, start, end, gene_length, 
           group_name, kingdom, phylum, class, order, family, genus, species) %>%
    filter(!gene_name %in% c("BCCT"))

# Calculate the distance from homA for each nucleotide_accession
dat4 <- dat3 %>%
    group_by(nucleotide_accession) %>%
    mutate(start_homA = start[gene_name == "homA"],
           end_homA = end[gene_name == "homA"]) %>%
    ungroup() %>%
    mutate(distance_from_homA = NA)
for (i in 1:nrow(dat4)) {
    # Calculate the distance from homA for each gene
    dat4$distance_from_homA[i] <- closest_distance(dat4$start[i], dat4$end[i], dat4$start_homA[i], dat4$end_homA[i])
}

# Pull out just the alphas and gammas, factor the order level by the number of genomes with > 10 count of assembly, all other as "other"
dat5 <- dat4 %>%
    filter(class %in% c("Alphaproteobacteria", "Gammaproteobacteria")) %>%
    group_by(order) %>%
    mutate(count = n_distinct(assembly)) %>%
    ungroup() %>%
    mutate(grouped_order = case_when(
        is.na(order) ~ "All others",
        count >119 ~ order,
        TRUE ~ "All others")) %>%
    mutate(grouped_order = factor(grouped_order, levels = c(unique(order), "All others"))) %>%
    filter(gene_name != "homA") %>%
    filter(distance_from_homA < 10000) 

dat_count <- dat5 %>%
    group_by(grouped_order) %>%
    summarise(count = n_distinct(assembly)) %>%
    ungroup() %>%
    arrange(desc(count))

counts <- dat5 %>%
    group_by(class, grouped_order) %>%
    summarise(n = n_distinct(assembly),
              per_homeE = sum(gene_name == "homE")/n*100)%>%
    ungroup() 

facet_labels <- counts %>%
    # first add gamma or alpha symbol after order name, then add count
    mutate(label = case_when(
        grouped_order == "All others" ~ paste0("All others\n(n = ", n, ", ", round(per_homeE, 0), "% w/homE)"),
        class == "Alphaproteobacteria" ~ paste0(grouped_order, "\n(α, n = ", n,  ", ", round(per_homeE, 0), "% w/homE)"),
        class == "Gammaproteobacteria" ~ paste0(grouped_order, "\n(γ, n = ", n, ", ", round(per_homeE, 0), "% w/homE)"))) %>%
    select(grouped_order, label) %>%
    deframe()


# Plot box plot, with x-axis as gene_name, y-axis as distance_from_homA, and facet as group_name
g <- ggplot(dat5 %>%
                filter(class %in% c("Alphaproteobacteria", "Gammaproteobacteria")), 
            aes(x = gene_name, y = distance_from_homA)) +
    geom_boxplot(outliers = TRUE, outlier.size = 0.1, width = 0.5) +
    facet_wrap(~ grouped_order, ncol = 5, 
               labeller = labeller(grouped_order = facet_labels)) +
    theme_bw() +
    labs(y = "Distance from homA (bp)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 7),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") 
g

# Plot a box plot of genome size and presence/absense of homE by genome
dat_genome_homE <- dat5 %>%
    select(assembly, nucleotide_accession, gene_name) %>%
    group_by(assembly, nucleotide_accession) %>%
    summarise(homE = sum(gene_name == "homE")) %>%
    ungroup() %>%
    mutate(homE = case_when(homE > 0 ~ "homE present",
                             homE == 0 ~ "homE absent") %>%
               factor(levels = c("homE present", "homE absent"))) %>%
    left_join(dat_genome_size, by = c("assembly"))

t_test_result <- t.test(genome_size_bp ~ homE, data = dat_genome_homE)

g2 <- ggplot(dat_genome_homE, aes(x = factor(homE), y = genome_size_bp/1000000)) +
    geom_boxplot(outliers = FALSE, outlier.size = 0.1, width = 0.5) +
    geom_jitter(size = 0.1, alpha = 0.1, width = 0.2) +
    # Add annotation of number of genomes in each group
    geom_text(data = dat_genome_homE %>%
                  group_by(homE) %>%
                  summarise(n = n_distinct(assembly)) %>%
                  ungroup() %>%
                  mutate(label = paste0("n = ", n)),
              aes(label = label, y = 10.0, x = factor(homE)), 
              size = 2.5, hjust = 0.5) +
    # Add annotaiton of results of t-test
    geom_text(data = data.frame(x = 1.5, y = 10.5, label = 
                                    ifelse(t_test_result$p.value < 0.001, 
                                           "t-test: p < 0.001",
                                        paste0("t-test: p = ", round(t_test_result$p.value, 5)))),
              aes(x = x, y = y, label = label), 
              size = 2.5, hjust = 0.5) +
    theme_bw() +
    labs(y = "Genome size (Mb)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 7),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")    
# Combine the two plots
g_combined <- g + g2 + plot_layout(ncol = 2, widths = c(3, 1)) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 7))
g_combined

# Save the plot
ggsave("figures/homE/Figure_SX8_combined_boxplots.pdf", 
       g_combined, 
       width = 10, height = 4, units = "in", dpi = 300,
       device = cairo_pdf)
