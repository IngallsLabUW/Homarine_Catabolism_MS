library(RaMS)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)

# GET FILES ------
## Set dirs ------
data_dir <- here("data", "intermediate", "metabolomics", "rpom")
meta_data_dir <- here("data", "raw", "metabolomics", "rpom")
raw_dat_dir <- here("data", "raw", "metabolomics", "rpom", "particulate", "mzml")
fig_dir <- here("figures", "exploratory", "metabolomics", "obi1", "chromats")


## Get data ------
# dat_long_filename <- here(data_dir, "combined_long_dat_scaled.csv")
# dat_MF_filename <- here(data_dir, "combined_MF_info.csv")
sample_key <- read_csv(here(meta_data_dir, "SampleList_EpicFate_Particulate.csv"), show_col_types = FALSE, col_names = TRUE)
# dat_long <- read_csv(dat_long_filename, show_col_types = FALSE)
# dat_MF <-
#   read_csv(dat_MF_filename, show_col_types = FALSE)

## Set mf_oi ------
mf_oi <- c(
    "Homarine",
    # "HILICPos_148", # n-methyl glutamine
    # "HILICPos_192", # unknown, but likely C7H9NO5 by mass and iso patt
    # "HILICPos_29", # unknown
    "HILICPos_153" # n-methyl glutamic acid
    )
bad_mf <- c("HILICNeg_123",  "HILICPos_180")
# Reference spectrum for n-methyl glutamic acid here: https://mona.fiehnlab.ucdavis.edu/spectra/display/CCMSLIB00005720398
# Reference spectrum for glutamine Same https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-BGC_Munich-RP001301

samples_oi <- sample_key %>%
  filter(
    # sample_fraction == "particulate",
    chromat == "HILIC",
    treatment == "Homarine" | treatment == "Glucose" | treatment == "Glucose and homarine" | treatment == "Homarine blank"
    |treatment == "Glucose blank" | treatment == "Glucose and homarine blank"
    # sample_set == "OBi1_set2" | replicate_name == "211202_Blk_BlankTube_3",
    # treatment == "Homarine" | treatment == "Glucose + NH4" | replicate_name == "211202_Blk_BlankTube_3"
  ) %>%
    mutate(filename = paste0(replicate_name, ".mzML")) %>%
    mutate(treatment = ifelse(is.na(treatment), "Blank", treatment))

# LOOP THROUGH EACH mf_oi_i-------
# for (mf_oi_i in mf_oi){
# mf_info <- dat_MF %>%
#   filter(mass_feature == mf_oi_i)
# if (mf_info$z == 1) {
#   data_subdir <- "HILIC_positive_mzml"
# } else {
#   data_subdir <- "HILIC_negative_mzml"
# }

data_subdir <- "HILIC_positive_mzml"

# all_ms_files <- list.files(here(raw_dat_dir, data_subdir))
all_ms_files <- list.files(here(raw_dat_dir, data_subdir))
msdata_files <- here(raw_dat_dir, data_subdir, all_ms_files[all_ms_files %in% samples_oi$filename])
dda_files <- here(raw_dat_dir, data_subdir, all_ms_files[str_detect(all_ms_files, "_DDA")])
msdata <- grabMSdata(files = msdata_files, grab_what = c("MS1"))

rt_dat <- msdata$MS1 %>% select(rt, filename) %>% distinct() %>%
    filter(rt %between% c(mf_info$RT-2.5, mf_info$RT+2.5))
ms1_data <- msdata$MS1[mz %between% pmppm(mf_info$mz, 20) & rt %between%c(mf_info$RT-3, mf_info$RT+3)] %>%
    full_join(rt_dat, by = join_by(rt, filename)) %>%
    left_join(samples_oi %>%
                  select(filename, treatment),
              by = join_by(filename)) %>%
    mutate(int = ifelse(is.na(int), 0, int))

g_ms1 <- ggplot(ms1_data) +
    geom_line(aes(x = rt, y = int, group  = filename)) +
    facet_wrap(~treatment, ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    geom_vline("vline", xintercept = mf_info$RT, linetype = "dashed") +
    scale_y_continuous(expand = c(0, NA)) +
    labs(y = "Intensity",
         x = "Retention Time (min)")

msdata_dda <- grabMSdata(files = dda_files, grab_what = c("MS2"))
ms2data <-  msdata_dda$MS2[premz %between% pmppm(mf_info$mz, 50) & rt %between%c(mf_info$RT-1, mf_info$RT+1)] %>%
    mutate(mz_round = round(fragmz, digits = 3)) %>%
    group_by(mz_round) %>%
    summarise(int = sum(int)) %>%
    ungroup() %>%
    mutate(int_norm = int/max(int)) %>%
    filter(int_norm > 0.01)
if (nrow(ms2data) > 0) {
  g_ms2 <- ggplot(ms2data) +
      geom_segment(
          aes(
              x = mz_round,
              xend = mz_round,
              y = 0,
              yend = int_norm
              ),
          alpha = 0.5,
          ) +
      geom_text_repel(
          aes(
              x = mz_round,
              y = int_norm,
              label = mz_round
              ),
          nudge_y = 0.05,
          size = 2.5
          ) +
      theme_bw()+
      labs(
           y = NULL,
           x = "m/z") +
      scale_y_continuous(expand = c(0, NA),
                         limits = c(0, 1.07)) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  g_save <- g_ms1 +
      g_ms2 +
      plot_layout(ncol = 2, widths = c(3, 1)) +
      plot_annotation(
          title = mf_oi_i,
          caption = paste0(
              "m/z: ",
              round(mf_info$mz, digits = 4),
              ", RT: ",
              round(mf_info$RT, digits = 2),
              " min, z: ",
              mf_info$z
              )
          )

  save_width <- 9
} else {
  g_save <- g_ms1 +
      plot_annotation(
          title = mf_oi_i,
          caption = paste0(
              "m/z: ",
              round(mf_info$mz, digits = 4),
              ", RT: ",
              round(mf_info$RT, digits = 2),
              " min, z: ",
              mf_info$z
          )
      )

  save_width <- 6
}
ggsave(here(fig_dir, paste0(mf_oi_i, ".pdf")), g_save, width = save_width, height = 5)


#RDisop for isotope envelopes
#plot(int ~ mz, type = "h", data = ms1_dat, ylab = "Intensity", xlab = "Fragment m/z")
#molecules <- decomposeMass(188.05644, ppm = 5)
