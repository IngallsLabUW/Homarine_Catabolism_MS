library(tidyverse)
library(patchwork)

# Read in data ------
growth.class.file <- "data/intermediate/isolategrowth/Isolates_hom_growth_classification.csv"
isolate.gene.file <- "data/isolate_platereader_data/Homarine_isolate_gene_info_042525.csv"

######## Pull in homarine catabolism gene information
iso.gene.info <- read_csv(isolate.gene.file) %>%
  select(-Gene) %>%
  pivot_longer(names_to = "gene", values_to = "presence", cols = `homE`:`soxB`) %>%
  mutate(presence = replace_na(presence, 0)) %>%
  mutate(type = "gene") %>%
  rename("isolate" = Isolate) %>%
  filter(!isolate == "G4i35")

### Pull in isolate growth data:
gene.iso.growth.class <- read_csv(growth.class.file) %>%
  mutate(isolate = case_when(
    isolate == "DSS-3" ~ "DSS3",
    TRUE ~ isolate
  )) %>%
  filter(isolate %in% iso.gene.info$isolate) %>%
  mutate(var = "growth") %>%
  mutate(hom_growth_status = case_when(
    hom_growth_status == "maybe_grower" ~ "non_grower",
    TRUE ~ hom_growth_status
  ))


### Combine datasets:
gene.plot.dat <- iso.gene.info %>%
  left_join(., gene.iso.growth.class %>%
    select(isolate, hom_growth_status)) %>%
  mutate(Family = as.factor(str_trim(Family, "left"))) %>%
  mutate(Class = as.factor(str_trim(Class, "left"))) %>%
  mutate(Family = fct_relevel(Family, c(
    "Micrococcales", "Mycobacteriales",
    "Roseobacteraceae", "Phyllobacteriaceae",
    "Alteromonadaceae", "Halomonadaceae",
    "Moraxellaceae", "Pseudoalteromonadaceae",
    "Pseudomonadaceae", "Vibrionaceae"
  )))

iso.order <- gene.plot.dat %>%
  select(isolate, hom_growth_status, Class, Family) %>%
  unique() %>%
  arrange(hom_growth_status, Class, Family) %>%
  mutate(order = row_number())

gene.plot.dat.order <- left_join(gene.plot.dat, iso.order)

# Growth Classification Heatmap --------------------------------------

### growth classification heatmap:
hom.growth.hm.dat <- gene.plot.dat.order %>%
  select(isolate, hom_growth_status, order)

growth.plot <- ggplot(hom.growth.hm.dat, aes(x = "Growth", y = reorder(isolate, -order))) +
  geom_tile(color = "black", width = 1, fill = "white") +
    # Add + sign to indicate growth
    geom_text(aes(label = ifelse(hom_growth_status == "grower", "+", "")), size = 5, color = "black") +
  #scale_fill_manual(values = c("#f8b767", "#fffcf9")) +
  theme_minimal() +
  ylab("Isolate") +
  theme(
    plot.margin = margin(
      t = 0, # Top margin
      r = 0, # Right margin
      b = 0, # Bottom margin
      l = 0
    ),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.8),
    axis.title.x = element_blank()
  ) +
  coord_fixed() +
  labs(fill = "Homarine Growth")
growth.plot


###

# Make Gene Presence Heatmaps --------------------------------------


# homarine genes
gene.plot.dat.order.hom <- gene.plot.dat.order %>%
  filter(str_detect(gene, "hom"))

hom.gene.plot <- ggplot(gene.plot.dat.order.hom, aes(x = gene, y = reorder(isolate, -order))) +
  geom_tile(color = "black", width = 1, fill="white") +
    geom_text(aes(label = ifelse(presence == 1, "+", "")), size = 5, color = "black") +
  #scale_fill_manual(values = c("white", "#b0a5cc")) +
  theme_minimal() +
  ylab("Isolate") +
  xlab("Gene") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(
      t = 0, # Top margin
      r = 0, # Right margin
      b = 0, # Bottom margin
      l = 0
    ),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.8),
    axis.title.x = element_blank()
  ) +
  coord_fixed() +
  labs(fill = "Gene Presence")
hom.gene.plot


# mgd genes
gene.plot.dat.order.mgd <- gene.plot.dat.order %>%
  filter(str_detect(gene, "mgd"))

mgd.gene.plot <- ggplot(gene.plot.dat.order.mgd, aes(x = gene, y = reorder(isolate, -order))) +
  geom_tile(color = "black", width = 1, fill = "white") +
    geom_text(aes(label = ifelse(presence == 1, "+", "")), size = 5, color = "black") +
  scale_fill_manual(values = c("white", "#b0a5cc")) +
  theme_minimal() +
  ylab("Isolate") +
  xlab("Gene") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(
      t = 0, # Top margin
      r = 0, # Right margin
      b = 0, # Bottom margin
      l = 0
    ),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.8)
  ) +
  coord_fixed() +
  labs(fill = "Gene Presence")
mgd.gene.plot



# sox genes
gene.plot.dat.order.sox <- gene.plot.dat.order %>%
  filter(str_detect(gene, "sox"))

sox.gene.plot <- ggplot(gene.plot.dat.order.sox, aes(x = gene, y = reorder(isolate, -order))) +
  geom_tile(color = "black", width = 1, fill="white") +
    geom_text(aes(label = ifelse(presence == 1, "+", "")), size = 5, color = "black") +
  scale_fill_manual(values = c("white", "#b0a5cc")) +
  theme_minimal() +
  ylab("Isolate") +
  xlab("Gene") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(
      t = 0, # Top margin
      r = 0, # Right margin
      b = 0, # Bottom margin
      l = 0
    ),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.8)
  ) +
  coord_fixed() +
  labs(fill = "Gene Presence")
sox.gene.plot








# Make Taxonomy Heatmaps --------------------------------------

# Class data:
class.plot <- ggplot(gene.plot.dat.order, aes(x = "Class", y = reorder(isolate, -order), fill = Class)) +
  geom_tile(color = "black", width = 1) +
  #  scale_fill_continuous(low = "white", high = "gray40") +
  scale_fill_manual(values = c("#f9c482", "#ecb1cf", "#6cccdc")) +
  theme_minimal() +
  ylab("Isolate") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(
      t = 0, # Top margin
      r = 0, # Right margin
      b = 0, # Bottom margin
      l = 0
    ),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.8),
    axis.title.x = element_blank()
  ) +
  coord_fixed()
class.plot

# Genus data:
genus.plot <- ggplot(gene.plot.dat.order, aes(x = "Genus", y = reorder(isolate, -order), fill = Family)) +
  geom_tile(color = "black", width = 1) +
  scale_fill_manual(values = c(
    "#f8b767", "#fcdcb5", "#dd8fb9", "#f5d6e5", "#16afca",
    "#26a6f3", "#42c7da", "#71caf8", "#a5dffb", "#d9f5f9"
  )) +
  #  scale_fill_tableau(palette = "Tableau 20") +
  #  scale_fill_continuous(low = "white", high = "gray40") +
  # scale_fill_manual(values = c("gray20", "white")) +
  theme_minimal() +
  ylab("Isolate") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(
      t = 0, # Top margin
      r = 0, # Right margin
      b = 0, # Bottom margin
      l = 0
    ),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.8),
    axis.title.x = element_blank()
  ) +
  coord_fixed()

genus.plot

####

# combine plots
combined.plot <-   (growth.plot | hom.gene.plot | mgd.gene.plot | sox.gene.plot | class.plot | genus.plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
# widths = c(1,-1 , 4, -1.5, 4,-1.5,
#                                4,-1.5 ,1, -0.8,1))

combined.plot


# save
#ggsave(combined.plot,
#  filename = "figures/isolates_genesearches/outputs/isolate_growth_gene_hm.png",
#  dpi = 600, height = 7.5, width = 7, scale = 1.25
#)



## Make supplemental Table:
iso.supp.table <- gene.plot.dat.order %>%
  select(isolate, hom_growth_status, Kingdom, Phylum, Class, Order, Family, Genus, Species, gene, presence) %>%
  pivot_wider(
    id_cols = c(isolate, hom_growth_status, Kingdom, Phylum, Class, Order, Family, Genus, Species),
    names_from = gene, values_from = presence
  ) %>%
  arrange(hom_growth_status, Class, Family)

write_csv(iso.supp.table, file = "tables/isolate_homarinegrowthstatus_and_genes.csv")


## Old Code:
# presence/absence of genes
