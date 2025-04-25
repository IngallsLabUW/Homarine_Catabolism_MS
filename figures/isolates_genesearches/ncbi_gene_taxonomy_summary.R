# Load required libraries ----
library(here)
library(tidyverse)
library(ggplot2)
library(janitor)

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

# Sunburst plot?
df <- "name  type    value
    foo   all     444
    foo   type1   123
    foo   type2   321
    bar   all     111
    bar   type3   111
    baz   all     999
    baz   type1   456
    baz   type3   543" %>% 
    read_table2() %>%
    filter(type != "all") %>%
    mutate(name = as.factor(name) %>% fct_reorder(value, sum)) %>%
    arrange(name, value) %>%
    mutate(type = as.factor(type) %>% fct_reorder2(name, value))

df2 <- dat4 %>%
    filter(!is.na(phylum)) %>%
    group_by(phylum, class, order) %>%
    summarise(value = n_distinct(assembly)) %>%
    ungroup() %>%
    mutate(class = case_when(
        phylum != "Pseudomonadota" ~ "Non-proteobacteria",
        TRUE ~ class)) 
    
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
    mutate(name = as.factor(name) %>% fct_reorder2(fill, value))%>%
    arrange(fill, name) %>%
    mutate(level = as.factor(level))

#TODO make points and labels separate than the columns?
library(ggrepel)
g <- final_df %>%
    ggplot(aes(x = level, y = value, fill = fill, alpha = level)) +
    geom_col(
        width = 1, color = "gray90", size = 0.25, position = position_stack()
        ) +
    geom_point(
        data = . %>% filter(level == 2), size = 0.5, position = position_stack(vjust = 0.5),
        alpha = 1,
    ) +
    geom_text(
        data = . %>% filter(level == 2),
        aes(label = name), size = 2.5, position = position_stack(vjust = 0.5),
                    alpha = 1,
                    ) +
    coord_polar(theta = "y") +
    scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 0.7), guide = F) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    #scale_fill_brewer(palette = "Dark2", na.translate = F) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    # Take out the ring around the plot
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(rep(-0.5, 4), "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
    ) 
g
