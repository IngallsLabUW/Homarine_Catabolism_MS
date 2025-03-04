



library(tidyverse)
library(patchwork)





#Inputs
growth.class.file <- "data/intermediate/isolategrowth/Isolates_hom_growth_classification.csv"
isolate.gene.file <- "data/isolate_platereader_data/Homarine_isolate_info_112024.csv"





########Pull in homarine catabolism gene information 
iso.gene.info <- read_csv(isolate.gene.file) %>%
  pivot_longer(names_to = "var", values_to = "presence", cols = `IEKGNDPE_02314`:`SPO3192/IEKGNDPE_02320`) %>%
  mutate(presence = replace_na(presence, 0)) %>%
  mutate(type = "gene") %>%
  rename("isolate" = Isolate)








# Make Growth Classification Heatmap --------------------------------------

### growth classification heatmap:
gene.iso.growth.class <- read_csv(growth.class.file) %>%
  mutate(isolate = case_when(isolate == "DSS-3" ~ "DSS3",
                             TRUE ~ isolate)) %>%
  filter(isolate %in% iso.gene.info$isolate) %>%
  mutate(var = "growth")

growth.plot <- ggplot(gene.iso.growth.class, aes(x = "Growth", y = reorder(isolate, -as.numeric(as.factor(hom_growth_status))), fill = hom_growth_status)) +
  geom_tile(color = "black", width = 1) +
  scale_fill_manual(values = c("#f8b767", "#fcdcb5", "#fffcf9")) +
  theme_minimal() +
  ylab("Isolate") +
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8),
        axis.title.x = element_blank()) +
  coord_fixed()
growth.plot


###

##Growth data:
gene.plot.dat <- iso.gene.info %>%
  left_join(., gene.iso.growth.class %>%
              select(isolate, hom_growth_status)) %>%
  mutate(Family = as.factor(str_trim(Family, "left"))) %>%
  mutate(Class = as.factor(str_trim(Class, "left"))) %>%
  mutate(Family = fct_relevel(Family, c("Micrococcales", "Mycobacteriales",
                                        "Roseobacteraceae", "Phyllobacteriaceae",
                                        "Alteromonadaceae", "Halomonadaceae",
                                        "Moraxellaceae", "Pseudoalteromonadaceae",
                                        "Pseudomonadaceae", "Vibrionaceae"))) %>%
  mutate(var = case_when(var == "IEKGNDPE_02314" ~ "HomE",
                         var == "SPO3186/IEKGNDPE_02317" ~ "DddT",
                         var == "SPO3188/IEKGNDPE_02315" ~ "homA",
                         var == "SPO3189/IEKGNDPE_02316" ~ "homB",
                         var == "SPO3190/IEKGNDPE_02318" ~ "homC",
                         var == "SPO3191/IEKGNDPE_02319" ~ "homD",
                         var == "SPO3192/IEKGNDPE_02320" ~ "homR"))




# Make Gene Presence Heatmap --------------------------------------

#presence/absence of genes
gene.plot <- ggplot(gene.plot.dat, aes(x = var, y = reorder(isolate, -as.numeric(as.factor(hom_growth_status))), fill = as.factor(presence))) +
  geom_tile(color = "black", width = 1) +
  scale_fill_manual(values = c("white", "#b0a5cc")) +
  # scale_fill_continuous(low = "white", high = "gray40") +
  # scale_fill_manual(values = c("gray20", "white")) +
  theme_minimal() +
  ylab("Isolate")  +
  xlab("Gene") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8)) +
  coord_fixed()
gene.plot




# Make Taxonomy Heatmaps --------------------------------------

#Class data:
class.plot <- ggplot(gene.plot.dat, aes(x = "Class", y = reorder(isolate, -as.numeric(as.factor(hom_growth_status))), fill = Class)) +
  geom_tile(color = "black", width = 1) +
  #  scale_fill_continuous(low = "white", high = "gray40") +
  scale_fill_manual(values = c("#f8b767", "#ecb1cf", "#6cccdc")) +
  theme_minimal() +
  ylab("Isolate") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8),
        axis.title.x = element_blank()) +
  coord_fixed()
class.plot

#Genus data:
genus.plot <- ggplot(gene.plot.dat, aes(x = "Genus", y = reorder(isolate, -as.numeric(as.factor(hom_growth_status))), 
                                        fill = Family)) +
  geom_tile(color = "black", width = 1) +
  scale_fill_manual(values = c("#f9c482", "#fcdcb5", "#dd8fb9", "#f5d6e5", "#16afca", 
                               "#26a6f3","#42c7da", "#71caf8", "#a5dffb", "#d9f5f9")) +
  #  scale_fill_tableau(palette = "Tableau 20") +
  #  scale_fill_continuous(low = "white", high = "gray40") +
  # scale_fill_manual(values = c("gray20", "white")) +
  theme_minimal() +
  ylab("Isolate") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8),
        axis.title.x = element_blank()) +
  coord_fixed()

genus.plot

####

#work on comining plots:




combined.plot <- (growth.plot | plot_spacer() | gene.plot | plot_spacer() | class.plot | plot_spacer() | genus.plot) +
  plot_layout(guides = 'collect', widths = c(1,-1.7,4,-1.7,1,-0.8,1))

combined.plot























































