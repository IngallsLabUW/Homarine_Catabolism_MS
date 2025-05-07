library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(forcats)

#Step 1: Load and merge metadata with operon cluster info
ncbi <- read_csv("ncbi_genome_summary.csv", skip = 1, show_col_types = FALSE)
operons <- read_csv("homarine_operon_cluster_summary_with_true_canonical.csv", show_col_types = FALSE)
genome_sizes <- read_csv("genome_sizes_datasets.csv", show_col_types = FALSE)

#Remove problematic genomes
samples_to_remove <- c("GCA_003942065.1", "GCA_003942045.1", "GCA_000987305.1")

genome_sizes <- genome_sizes %>%
  filter(!Assembly_Accession %in% samples_to_remove)

#Standardize column names
colnames(ncbi)[1] <- "Assembly"
colnames(ncbi)[3:7] <- c("homA_location", "homB_location", "homC_location", "homD_location", "homE_location")
colnames(genome_sizes)[1] <- "Assembly"

#Merge all data
df <- ncbi %>%
  left_join(operons %>% select(Assembly, gene_cluster, canonical_cluster), by = "Assembly") %>%
  left_join(genome_sizes, by = "Assembly")

#Canonicalize fallback if canonical_cluster is missing
canonicalize_cluster <- function(cluster) {
  if (is.na(cluster)) return(NA)
  cluster <- gsub("hom", "", cluster)
  if (nchar(cluster) == 0) return(NA)
  genes <- strsplit(cluster, "")[[1]]
  forward <- paste(genes, collapse = "")
  reverse <- paste(rev(genes), collapse = "")
  return(paste0("hom", min(forward, reverse)))
}

df <- df %>%
  mutate(canonical_cluster = ifelse(is.na(canonical_cluster),
                                    sapply(gene_cluster, canonicalize_cluster),
                                    canonical_cluster))

#Extract gene start and end coordinates
get_start <- function(x) {
  if (is.na(x) || !grepl("-", x)) return(NA_real_)
  suppressWarnings(as.numeric(sub("-.*", "", x)))
}
get_end <- function(x) {
  if (is.na(x) || !grepl("-", x)) return(NA_real_)
  suppressWarnings(as.numeric(sub(".*-", "", x)))
}

df <- df %>%
  mutate(
    homA_start = sapply(homA_location, get_start),
    homA_end   = sapply(homA_location, get_end),
    homB_start = sapply(homB_location, get_start),
    homB_end   = sapply(homB_location, get_end),
    homC_start = sapply(homC_location, get_start),
    homD_start = sapply(homD_location, get_start),
    homD_end   = sapply(homD_location, get_end),
    homE_start = sapply(homE_location, get_start),
    homE_end   = sapply(homE_location, get_end)
  )

#Rebuild gene cluster based on gene order
rebuild_cluster <- function(A, B, C, D, E) {
  genes <- list(A = A, B = B, C = C, D = D, E = E)
  genes <- genes[!sapply(genes, is.na)]
  if (length(genes) == 0) return(NA)
  ordered <- names(sort(unlist(genes)))
  paste(ordered, collapse = "")
}

df <- df %>%
  mutate(gene_cluster_reconstructed = mapply(rebuild_cluster,
                                             homA_start, homB_start,
                                             homC_start, homD_start,
                                             homE_start))

#Compute intergenic gap between homA and homB
df <- df %>%
  mutate(
    AB_gap = case_when(
      !is.na(homA_end) & !is.na(homB_start) & homA_end < homB_start ~ homB_start - homA_end,
      !is.na(homB_end) & !is.na(homA_start) & homB_end < homA_start ~ homA_start - homB_end,
      TRUE ~ 0
    )
  )

#Coordinate-based adjacency
adjacency_threshold <- 900
df_coord <- df %>%
  filter(!is.na(AB_gap)) %>%
  mutate(adjacent = AB_gap < adjacency_threshold)

coord_stats <- df_coord %>%
  summarise(
    method = "Coordinate-based",
    total = n(),
    adjacent = sum(adjacent),
    percent_adjacent = round(100 * sum(adjacent) / n(), 1)
  )

#Operon-based adjacency from reconstructed order
is_AB_adjacent <- function(cluster) {
  if (is.na(cluster)) return(FALSE)
  genes <- strsplit(cluster, "")[[1]]
  any(sapply(seq_along(genes[-1]), function(i) {
    setequal(c(genes[i], genes[i + 1]), c("A", "B"))
  }))
}

operon_stats <- df %>%
  filter(grepl("A", gene_cluster_reconstructed), grepl("B", gene_cluster_reconstructed)) %>%
  mutate(adjacent = sapply(gene_cluster_reconstructed, is_AB_adjacent)) %>%
  summarise(
    method = "Operon-order",
    total = n(),
    adjacent = sum(adjacent),
    percent_adjacent = round(100 * sum(adjacent) / n(), 1)
  )

comparison <- bind_rows(coord_stats, operon_stats)
print(comparison)

#Compare methods and highlight discrepancies
df <- df %>%
  mutate(coord_adjacent = AB_gap < adjacency_threshold,
         operon_adjacent = sapply(gene_cluster_reconstructed, is_AB_adjacent))

discrepancies <- df %>%
  filter(coord_adjacent != operon_adjacent) %>%
  select(Assembly, gene_cluster, canonical_cluster, gene_cluster_reconstructed,
         AB_gap, coord_adjacent, operon_adjacent)

discrepancies_filtered <- discrepancies %>%
  filter(!(coord_adjacent == TRUE & operon_adjacent == FALSE & AB_gap < 1000),
         !(coord_adjacent == FALSE & operon_adjacent == TRUE & AB_gap <= 1000))

print(discrepancies_filtered)

original_coord_adjacent_count <- sum(df$coord_adjacent, na.rm = TRUE)
n_removed <- sum(discrepancies$coord_adjacent & !discrepancies$operon_adjacent & discrepancies$AB_gap < 1000, na.rm = TRUE)
n_added   <- sum(!discrepancies$coord_adjacent & discrepancies$operon_adjacent & discrepancies$AB_gap > 1000, na.rm = TRUE)
adjusted_total <- original_coord_adjacent_count - n_removed + n_added
cat("\nAdjusted total count of homA-B adjacent genomes:", adjusted_total, "\n")

#Strand-aware logic for homE relative to homA and homD
df <- df %>%
  mutate(
    strand = case_when(
      !is.na(homA_start) & !is.na(homB_start) & homA_start < homB_start ~ "forward",
      !is.na(homA_start) & !is.na(homB_start) & homA_start > homB_start ~ "reverse",
      TRUE ~ NA_character_
    )
  )

df <- df %>%
  mutate(
    homE_position = case_when(
      strand == "forward" & !is.na(homE_end) & !is.na(homA_start) & homE_end < homA_start ~ "Before homA",
      strand == "forward" & !is.na(homE_start) & !is.na(homD_end) & homE_start > homD_end ~ "After homD",
      strand == "reverse" & !is.na(homE_start) & !is.na(homA_end) & homE_start > homA_end ~ "Before homA",
      strand == "reverse" & !is.na(homE_end) & !is.na(homD_start) & homE_end < homD_start ~ "After homD",
      !is.na(homE_start) ~ "Other position",
      TRUE ~ "Not present"
    )
  )

#Count and print homE position summary
homE_position_summary <- df %>%
  count(homE_position, sort = TRUE)

print(homE_position_summary)

# #Subset the data into three categories based on homE position
# df_before_homA <- df %>% filter(homE_position == "Before homA")
# df_after_homD  <- df %>% filter(homE_position == "After homD")
# df_not_present <- df %>% filter(homE_position == "Not present")

#Subset data into two groups: present and not present
df <- df %>%
  mutate(homE_grouped = ifelse(homE_position == "Not present", "Not present", "Present"))
df_present     <- df %>% filter(homE_grouped == "Present")
df_not_present <- df %>% filter(homE_grouped == "Not present")


#Bar plot of homE positions
ggplot(homE_position_summary, aes(x = homE_position, y = n, fill = homE_position)) +
  geom_col() +
  labs(
    title = "Location of homE Relative to Operon",
    x = "homE Position",
    y = "Number of Genomes"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#Taxonomic distribution of homE positions
homE_taxa_summary <- df %>%
  filter(homE_position %in% c("Before homA", "After homD", "Not present")) %>%
  group_by(homE_position, Class, Order) %>%
  tally(name = "Count") %>%
  arrange(homE_position, desc(Count))

print(homE_taxa_summary)

ggplot(homE_taxa_summary, aes(x = fct_infreq(Class), y = Count, fill = homE_position)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Taxonomic Distribution of homE Positions",
    x = "Class",
    y = "Number of Genomes",
    fill = "homE Position"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Perform t-test to compare genome sizes
 df <- df %>%
   mutate(homE_present = homE_position != "Not present")

t_test_result <- t.test(Genome_Size_bp ~ homE_present, data = df)

print(t_test_result)

df %>%
  group_by(homE_present) %>%
  summarise(mean_genome_size = mean(Genome_Size_bp, na.rm = TRUE),
            sd_genome_size = sd(Genome_Size_bp, na.rm = TRUE),
            n = n())


# Add gene_arrangement and homE_present to the original ncbi table
ncbi_export <- ncbi %>%
  left_join(df %>% select(Assembly, canonical_cluster), by = "Assembly") %>%
  mutate(
    canonical_cluster = ifelse(!is.na(canonical_cluster), paste0("hom", canonical_cluster), NA)
  ) %>%
  relocate(canonical_cluster, .after = homE_location) %>%
  rename(gene_arrangement = canonical_cluster)

# Write to CSV
write_csv(ncbi_export, "Table_SX8_ncbi_genome_summary_mod.csv")

#Plot with t-test added
pval <- signif(t_test_result$p.value, digits = 3)
mean_diff <- round(diff(t_test_result$estimate), 0)

ggplot(df, aes(x = fct_infreq(homE_grouped), y = Genome_Size_bp, fill = homE_grouped)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Genome Size by homE Presence",
    x = "homE Presence",
    y = "Genome Size (bp)"
  ) +
  annotate("text", x = 1.5, y = max(df$Genome_Size_bp, na.rm = TRUE) * 0.95,
           label = paste0("t-test p = ", pval, "\nΔ mean = ", mean_diff, " bp"),
           size = 4, hjust = 0) +
  theme_minimal() +
  theme(legend.position = "none")
