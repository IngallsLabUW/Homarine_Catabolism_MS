#Load packages:
library(tidyverse)
library(patchwork)
library(forcats)

#Define inputs:
## This file is "results from organisms from which I only found co-located in a their genome three out of the six known proteins involved in homarine uptake and degradation in OBi1." (from Oscar)
ncbi.hom.file <- "data/intermediate/ncbi_genome_counts/ncbi_homarine_operons_loci_with_taxonomy.txt"
ncbi.order.counts.file <- "data/intermediate/ncbi_genome_counts/ncbi_genomes_per_order.csv"
ncbi.class.counts.file <- "data/intermediate/ncbi_genome_counts/ncbi_genomes_per_class.csv"
 

#Organize data:

#Load data and rename homarine genes:
ncbi_dat_og <- read_tsv(ncbi.hom.file,show_col_types = FALSE) 

ncbi_dat1 <- ncbi_dat_og %>%
    mutate(gene = case_when(query %in% c("OBi1_02315", "SPO3188") ~ "homA",
                            query %in% c("OBi1_02316") ~ "homB",
                            # Obi1_02371 is the transporter
                            query %in% c("OBi1_02318", "SPO3190") ~ "homC",
                            query %in% c("OBi1_02319", "SPO3191") ~ "homD",
                            #query %in% c("OBi1_02320") ~ "homR",
                            query %in% c("OBi1_02314") ~ "homE")) %>%
    filter(!is.na(gene)) %>%
    mutate(present = TRUE) %>%
    select(gene, present, Organism, class, order, family, genus, Assembly) %>%
    rename(assembly = Assembly) %>%
    unique() %>%
    filter(!is.na(class))



############make dataframe of all organisms and all genes:

#make dataframe of genes
genes <- ncbi_dat1 %>%
    select(gene) %>%
    unique()

#make dataframe of organisms 
orgs <- ncbi_dat1 %>%
    select(assembly, Organism, class, order, gene) %>% 
    unique() %>%
    rename(gene.present = gene) 

#create presences/absence gene data
gene.pa.dat <- cross_join(orgs, genes) %>%
    mutate(gene.status = case_when(gene.present == gene ~ 1,
                                   TRUE ~ 0)) %>%
    group_by(assembly, Organism, class, order, gene) %>%
    reframe(gene.status = sum(gene.status)) 


#identify organisms with core 4 catabolic genes
core.catab <- gene.pa.dat %>%
    filter(gene %in% c("homA", "homB", "homC", "homD")) %>%
    group_by(assembly) %>%
    reframe(core.path.score = sum(gene.status)) %>%
    mutate(core.catab.path = case_when(core.path.score == 4 ~ "Yes",
                                       TRUE ~ "No")) %>%
    select(assembly, core.catab.path)


###combine presence/absence, core 4, and total ncbi genomes data:
final.pa.dat <- left_join(gene.pa.dat, core.catab, by = join_by(assembly)) %>%
    full_join(., read_csv(ncbi.order.counts.file, show_col_types = FALSE),
              by = join_by(order))


#get class counts:
class.counts <- read_csv(ncbi.class.counts.file)


####summarize selected groups:

###selected groups:

#Gammas:

#identify groups
gamma.groups <- final.pa.dat %>%
    filter(class == "Gammaproteobacteria") %>%
    select(assembly, core.catab.path, order, genomes) %>%
    distinct() %>%
    filter(!is.na(order)) %>%
    group_by(order) %>%
    reframe(order.count = n(),
            genomes = genomes, 
            percent = order.count/genomes*100) %>%
    unique() %>%
    mutate(fig.group = case_when(order.count > 20 ~ order,
                                 TRUE ~ "Other_Gammas")) %>%
    select(order, fig.group) %>%
    unique()

#organize figure data:
gamma.fig.dat <- final.pa.dat %>%
    filter(!is.na(order)) %>%
    filter(gene.status == 1) %>%
    filter(class == "Gammaproteobacteria") %>%
    left_join(., gamma.groups) %>%
    group_by(gene, order, fig.group, genomes) %>%
    reframe(group.count = n()) %>%
    group_by(gene, fig.group) %>%
    reframe(tot.count = sum(group.count),
            tot.genomes = sum(genomes)) %>%
    cross_join(class.counts %>%
                   filter(class == "Gammaproteobacteria") %>% 
                   select(genomes) %>% 
                   rename(tot.gamma.genomes = genomes)) %>%
    group_by(gene) %>%
    mutate(fig.group.genomes = case_when(fig.group == "Other_Gammas" ~ tot.gamma.genomes - 
                                             sum(tot.genomes[!fig.group == "Other_Gammas"]),
                                         TRUE ~ tot.genomes)) %>%
    group_by(gene, fig.group, tot.count, fig.group.genomes) %>%
    reframe(percent = tot.count/fig.group.genomes*100)

    
#Alphas:
#identify groups
alpha.groups <- final.pa.dat %>%
    filter(class == "Alphaproteobacteria") %>%
    select(assembly, core.catab.path, order, genomes) %>%
    distinct() %>%
    filter(!is.na(order)) %>%
    group_by(order) %>%
    reframe(order.count = n(),
            genomes = genomes, 
            percent = order.count/genomes*100) %>%
    unique() %>%
    mutate(fig.group = case_when(order.count > 20 ~ order,
                                 TRUE ~ "Other_Alphas")) %>%
    select(order, fig.group) %>%
    unique()

#organize figure data:
alpha.fig.dat <- final.pa.dat %>%
    filter(!is.na(order)) %>%
    filter(gene.status == 1) %>%
    filter(class == "Alphaproteobacteria") %>%
    left_join(., alpha.groups) %>%
    group_by(gene, order, fig.group, genomes) %>%
    reframe(group.count = n()) %>%
    group_by(gene, fig.group) %>%
    reframe(tot.count = sum(group.count),
            tot.genomes = sum(genomes)) %>%
    cross_join(class.counts %>%
                   filter(class == "Alphaproteobacteria") %>% 
                   select(genomes) %>% 
                   rename(tot.gamma.genomes = genomes)) %>%
    group_by(gene) %>%
    mutate(fig.group.genomes = case_when(fig.group == "Other_Alphas" ~ tot.gamma.genomes - 
                                             sum(tot.genomes[!fig.group == "Other_Alphas"]),
                                         TRUE ~ tot.genomes)) %>%
    group_by(gene, fig.group, tot.count, fig.group.genomes) %>%
    reframe(percent = tot.count/fig.group.genomes*100)


##All Others:

#identify groups
other.groups <- final.pa.dat %>%
    filter(!class %in% c("Gammaproteobacteria", "Alphaproteobacteria")) %>%
    filter(gene.status == 1) %>%
    filter(gene == "homB") %>%
    filter(!is.na(order)) %>%
    group_by(order) %>%
    reframe(order.count = n(),
            genomes = genomes, 
            percent = order.count/genomes*100) %>%
    unique() %>%
    mutate(fig.group = case_when(percent > 1 & order.count > 20 ~ order,
                                 TRUE ~ "Other")) %>%
    select(order, fig.group) %>%
    unique()

#organize figure data:
other.fig.dat <- final.pa.dat %>%
    filter(!is.na(order)) %>%
    filter(gene.status == 1) %>%
    filter(!class %in% c("Gammaproteobacteria", "Alphaproteobacteria")) %>%
    left_join(., other.groups) %>%
    group_by(gene, order, fig.group, genomes) %>%
    reframe(group.count = sum(gene.status)) %>%
    group_by(gene, fig.group) %>%
    reframe(tot.count = sum(group.count),
            tot.genomes = sum(genomes)) %>%
    cross_join(class.counts %>%
                   filter(!class %in% c("Gammaproteobacteria", "Alphaproteobacteria")) %>% 
                   select(genomes) %>% 
                   rename(tot.other.genomes = genomes)) %>%
    filter(!is.na(fig.group)) %>% 
    group_by(gene) %>%
    mutate(fig.group.genomes = case_when(fig.group == "Other" ~ sum(tot.other.genomes) - 
                                             sum(tot.genomes[!fig.group == "Other"]),
                                         TRUE ~ tot.genomes)) %>%
    group_by(gene, fig.group, tot.count, fig.group.genomes) %>%
    reframe(percent = tot.count/fig.group.genomes*100) %>%
    unique()





### combine all group data:
all.fig.dat <- rbind(gamma.fig.dat, alpha.fig.dat, other.fig.dat)


#Make heatmap figure:
library(viridis)

group.hm <- ggplot(all.fig.dat, aes(x = gene, y = reorder(fig.group, percent), fill = percent)) +
    geom_tile(color = "black") +
    scale_fill_viridis() +
   # scale_fill_gradient(low = "white", high = "darkblue") +
    labs(fill = "Percent of Genomes") +
    ylab("Order") +
    theme_test() +
    coord_fixed()
group.hm



#








 
    
    
alphas.select.dat <- final.pa.dat %>%
    filter(class == "Alphaproteobacteria") %>%
    filter(gene == "homB") %>%
    filter(!is.na(order)) %>%
    group_by(order) %>%
    reframe(order.count = n(),
            genomes = genomes, 
            percent = order.count/genomes*100) %>%
    unique() %>%
    mutate(fig.group = case_when(percent > 1 & order.count > 10 ~ order,
                                 TRUE ~ "Other_Alphas"))%>%
    group_by(fig.group) %>%
    reframe(tot.count = sum(order.count),
            tot.genomes = sum(genomes)) %>%
    cross_join(class.counts %>%
                   filter(class == "Alphaproteobacteria") %>% 
                   select(genomes) %>% 
                   rename(tot.gamma.genomes = genomes)) %>%
    mutate(fig.group.genomes = case_when(fig.group == "Other_Alphas" ~ tot.gamma.genomes - 
                                             sum(tot.genomes[!fig.group == "Other_Alphas"]),
                                         TRUE ~ tot.genomes)) %>%
    group_by(fig.group, tot.count, fig.group.genomes) %>%
    reframe(percent = tot.count/fig.group.genomes*100)


#Alphas



gamma.fig.dat <- final.pa.dat %>%
    filter(!is.na(order)) %>%
    filter(gene.status == 1) %>%
    filter(class == "Gammaproteobacteria") %>%
    left_join(., gamma.groups) %>%
    group_by(gene, order, fig.group, genomes) %>%
    reframe(group.count = n()) %>%
    group_by(gene, fig.group) %>%
    reframe(tot.count = sum(group.count),
            tot.genomes = sum(genomes)) %>%
    cross_join(class.counts %>%
                   filter(class == "Gammaproteobacteria") %>% 
                   select(genomes) %>% 
                   rename(tot.gamma.genomes = genomes)) %>%
    group_by(gene) %>%
    mutate(fig.group.genomes = case_when(fig.group == "Other_Gammas" ~ tot.gamma.genomes - 
                                             sum(tot.genomes[!fig.group == "Other_Gammas"]),
                                         TRUE ~ tot.genomes)) %>%
    group_by(gene, fig.group, tot.count, fig.group.genomes) %>%
    reframe(percent = tot.count/fig.group.genomes*100)











others.select.dat <- final.pa.dat %>%
    filter(!class %in% c("Gammaproteobacteria", "Alphaproteobacteria")) %>%
    filter(gene == "homB")  %>%
    filter(!is.na(order)) %>%
    group_by(order) %>%
    reframe(order.count = n(),
            genomes = genomes, 
            percent = order.count/genomes*100) %>%
    unique() %>%
    mutate(fig.group = case_when(percent > 1 & order.count > 10 ~ order,
                                 TRUE ~ "Other")) %>%
    group_by(fig.group) %>%
    reframe(tot.count = sum(order.count),
            tot.genomes = sum(genomes)) %>%
    cross_join(class.counts %>%
                   filter(!class %in% c("Gammaproteobacteria", "Alphaproteobacteria")) %>% 
                   select(genomes) %>% 
                   rename(tot.group.genomes = genomes)) %>%
    mutate(fig.group.genomes = case_when(fig.group == "Other" ~ sum(tot.group.genomes) - 
                                             sum(tot.genomes[!fig.group == "Other"]),
                                         TRUE ~ tot.genomes)) %>%
    group_by(fig.group, tot.count, fig.group.genomes) %>%
    reframe(percent = tot.count/fig.group.genomes*100) %>%
    unique()
        



full.fig.dat<- rbind(gamma.select.dat, alphas.select.dat, others.select.dat) %>%
    filter(!is.na(fig.group))






