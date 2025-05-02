library(readxl)
library(tidyverse)
library(patchwork)

# Set font size
fnt_size <- 5

# Utility: clean all columns in a dataframe of non-breaking spaces
clean_nbsp <- function(df) {
  df %>%
    mutate(across(where(is.character), ~str_replace_all(., "\u00A0", " ") |> str_squish()))
}

# Read in data ------
growth.class.file <- "Isolates_hom_growth_classification.csv"
isolate.gene.file <- "Homarine_isolate_gene_info_042525.csv"
ncbi.accession.file <- "2025_5_1_homarine_isolates_NCBI_accession_IDs.csv"
isolate.locations.file <- "Isolation_Locations_01142025.csv"

# Read gene IDs 
gene.id.df <- read_excel("isolates_homi_nmgda_sox_catabolic_genes_final.xlsx", skip = 1) %>%
  clean_nbsp()

# Set proper column name
names(gene.id.df)[1] <- "Isolate"

gene.id.trim <- gene.id.df %>% select(Isolate, homE:soxB)

gene.id.long <- gene.id.trim %>%
  pivot_longer(cols = -Isolate, names_to = "gene", values_to = "gene_id") %>%
  rename(isolate = Isolate)

######## Pull in homarine catabolism gene information
iso.gene.info <- read_csv(isolate.gene.file) %>%
  clean_nbsp() %>%
  select(-Gene) %>%
  pivot_longer(names_to = "gene", values_to = "presence", cols = homE:soxB) %>%
  mutate(presence = replace_na(presence, 0)) %>%
  mutate(type = "gene") %>%
  rename("isolate" = Isolate) %>%
  filter(isolate != "G4i35")

### Pull in isolate growth data:
gene.iso.growth.class <- read_csv(growth.class.file) %>%
  clean_nbsp() %>%
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
  left_join(gene.iso.growth.class %>% select(isolate, hom_growth_status), by = "isolate") %>%
  mutate(Family = as.factor(str_trim(Family, "left")),
         Class = as.factor(str_trim(Class, "left")),
         Family = fct_relevel(Family, c(
           "Micrococcales", "Mycobacteriales",
           "Roseobacteraceae", "Phyllobacteriaceae",
           "Alteromonadaceae", "Halomonadaceae",
           "Moraxellaceae", "Pseudoalteromonadaceae",
           "Pseudomonadaceae", "Vibrionaceae"
         )))

iso.order <- gene.plot.dat %>%
  select(isolate, hom_growth_status, Class, Family) %>%
  distinct() %>%
  arrange(hom_growth_status, Class, Family) %>%
  mutate(order = row_number())

gene.plot.dat.order <- left_join(gene.plot.dat, iso.order,
                                 by = c("isolate", "Class", "Family", "hom_growth_status"))

# Growth Classification Heatmap --------------------------------------
hom.growth.hm.dat <- gene.plot.dat.order %>%
  select(isolate, hom_growth_status, order)

growth.plot <- ggplot(hom.growth.hm.dat, aes(x = "Growth", y = reorder(isolate, -order))) +
  geom_tile(color = "black", width = 1, fill = "white") +
  geom_text(aes(label = ifelse(hom_growth_status == "grower", "+", "")), size = 3, color = "black") +
  theme_minimal() +
  ylab("Isolate") +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 0, l = -20),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.8, size = fnt_size),
    axis.text.y = element_text(size = fnt_size),
    axis.title.y = element_text(size = fnt_size),
    axis.title.x = element_blank()
  ) +
  coord_fixed() +
  labs(fill = "Homarine Growth")

# Combine plots
combined.plot <- (growth.plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save plot (optional)
# ggsave(combined.plot, filename = "figures/isolates_genesearches/outputs/isolate_growth_gene_hm.png",
#        dpi = 600, height = 7.5, width = 7, scale = 1.25)

# Create supplemental table --------------------------------------
iso.supp.table <- gene.plot.dat.order %>%
  select(isolate, hom_growth_status, Kingdom, Phylum, Class, Order, Family, Genus, Species, gene, presence) %>%
  pivot_wider(
    id_cols = c(isolate, hom_growth_status, Kingdom, Phylum, Class, Order, Family, Genus, Species),
    names_from = gene, values_from = presence
  ) %>%
  arrange(hom_growth_status, Class, Family)

# Read accession IDs
accession.df <- read_csv(ncbi.accession.file) %>%
  clean_nbsp()

# Add accession ID to supplemental table
iso.supp.table <- iso.supp.table %>%
  left_join(accession.df %>% select(isolate = Isolate, Genome_Accession), by = "isolate") %>%
  relocate(Genome_Accession, .after = isolate)

# Merge gene IDs into final output
gene.values.long <- iso.supp.table %>%
  pivot_longer(cols = homE:soxB, names_to = "gene", values_to = "presence")

gene.values.updated <- gene.values.long %>%
  left_join(gene.id.long, by = c("isolate", "gene")) %>%
  mutate(presence = case_when(
    presence == 1 ~ gene_id,
    presence == 0 ~ "0",
    TRUE ~ as.character(presence)
  )) %>%
  select(-gene_id)

iso.supp.table <- gene.values.updated %>%
  pivot_wider(names_from = gene, values_from = presence) %>%
  relocate(Genome_Accession, .after = isolate)

# Read and clean isolation locations
location.df <- read_csv(isolate.locations.file) %>%
  clean_nbsp() %>%
  select(isolate = isolate, Lat, Lon, Depth_Feature, Depth_m)

# Join location data to final table
iso.supp.table <- iso.supp.table %>%
  left_join(location.df, by = "isolate")

# Move location columns after gene IDs
location_cols <- c("Lat", "Lon", "Depth_Feature", "Depth_m")
iso.supp.table <- iso.supp.table %>%
  select(-all_of(location_cols), everything(), all_of(location_cols))

iso.supp.table <- iso.supp.table %>%
  mutate(across(everything(), ~ ifelse(. %in% c("0", NA), "", .)))

# Save final output
write_csv(iso.supp.table, file = "isolate_homarinegrowthstatus_and_genes_mod.csv")

