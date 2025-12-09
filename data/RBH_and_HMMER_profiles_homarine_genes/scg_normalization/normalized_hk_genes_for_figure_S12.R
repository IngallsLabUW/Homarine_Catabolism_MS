library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

setwd("C:/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/hmmr_profiles/big_tables")

# Load Data
hk <- read_csv("BigTable_G3_NS_subset_alpha_gamma_hk_genes.csv")
spogenes <- read_csv("SPO3188_SPO3189_combined_G3_revised_annotations_with_metadata_final_revised_copy.csv")

rep_cols <- c("raw_counts_A", "raw_counts_B", "raw_counts_C")

# SCG KOfams
scgs <- c(
  "K06942","K01889","K01887","K01875","K01883",
  "K01869","K01873","K01409","K03106","K03110"
)

target_classes <- c("Alphaproteobacteria","Gammaproteobacteria")


# 1. Filter to Alpha + Gamma ONLY

hk_filt <- hk %>% filter(class %in% target_classes)
spogenes_filt <- spogenes %>% filter(class %in% target_classes)


# 2. Prepare HK SCGs — keep REPLICATES SEPARATE
hk_long <- hk_filt %>%
  pivot_longer(
    cols = all_of(rep_cols),
    names_to = "Replicate",
    values_to = "ReadCount"
  ) %>%
  group_by(station, KOfam, Replicate) %>%
  summarize(ReadCount = sum(ReadCount), .groups = "drop")

# Keep only SCG KOfams
hk_scg <- hk_long %>% filter(KOfam %in% scgs)

# Compute SCG normalization PER replicate PER station
scg_factor <- hk_scg %>%
  group_by(station, Replicate) %>%
  summarize(SCG = median(ReadCount), .groups = "drop")


# 3. Prepare SPO counts — keep REPLICATES separate
spo_long <- spogenes_filt %>%
  pivot_longer(
    cols = all_of(rep_cols),
    names_to = "Replicate",
    values_to = "ReadCount"
  ) %>%
  group_by(station, hmm_profile, Replicate) %>%
  summarize(ReadCount = sum(ReadCount), .groups = "drop")


# 4. Normalize EACH REPLICATE independently
normalized <- spo_long %>%
  left_join(scg_factor, by = c("station", "Replicate")) %>%
  mutate(
    NormReadCount = ifelse(SCG > 0, ReadCount / SCG, NA_real_),
    Log2Norm = log2(NormReadCount)
  )


# 5. Add Latitude & Longitude
station_meta <- spogenes_filt %>%
  distinct(station, Latitude, Longitude)

normalized <- normalized %>%
  left_join(station_meta, by = "station")

# 6. Rename homA / homB
normalized <- normalized %>%
  mutate(
    Gene = case_when(
      hmm_profile == "OBi1_02315_SPO3188" ~ "homA",
      hmm_profile == "OBi1_02316_SPO3189" ~ "homB",
      TRUE ~ hmm_profile
    )
  )

# 7. Compute mean + SD across replicates
plot_df <- normalized %>%
  group_by(Gene, station, Latitude, Longitude) %>%
  summarize(
    MeanLog2 = mean(Log2Norm, na.rm = TRUE),
    SDLog2   = sd(Log2Norm, na.rm = TRUE),
    SEMLog2  = sd(Log2Norm, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(Latitude)

# 8. Plot with ERROR BARS (horizontal)
ggplot(plot_df, aes(x = MeanLog2, y = Latitude, color = Gene)) +
    geom_path(linewidth = 1.2, aes(group = Gene)) +
    geom_point(size = 3) +
    geom_errorbarh(
        aes(
            xmin = MeanLog2 - SDLog2,
            xmax = MeanLog2 + SDLog2
        ),
        height = 0.15,
        linewidth = 0.8
    ) +
    coord_flip() +
    theme_bw(base_size = 14) +
    labs(
        title = "Normalized Homarine Gene Expression",
        x = "Log2 Normalized Counts (per replicate, pooled Alpha+Gamma)",
        y = "Latitude"
    )


p <- ggplot(plot_df, aes(x = MeanLog2, y = Latitude, color = Gene)) +
    geom_path(linewidth = 1.2, aes(group = Gene)) +
    geom_point(size = 3) +
    geom_errorbarh(
        aes(
            xmin = MeanLog2 - SDLog2,
            xmax = MeanLog2 + SDLog2
        ),
        height = 0.15,
        linewidth = 0.8
    ) +
    coord_flip() +
    theme_bw(base_size = 6) +   # set font size here (5–7 recommended)
    labs(
        title = "Normalized Homarine Gene Expression",
        x = "Log2 Normalized Counts (per replicate)",
        y = "Latitude"
    )

# ggsave(
#     filename = "Homarine_Gene_Expression_plot.pdf",
#     plot = p,
#     width = 18,      # in cm
#     height = 12,     # adjust as needed
#     units = "cm",
#     dpi = 300
# )
