library(RaMS)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
source(here("analysis", "obi1_metabolomics", "functions", "plotting.R"))

# GET FILES ------
## Set dirs ------
data_dir <- here("data", "intermediate", "metabolomics", "obi1")
meta_data_dir <- here("data", "raw", "metabolomics", "obi1")
raw_dat_dir <- here("data", "raw", "metabolomics", "obi1", "particulate", "mzml")
fig_dir <- here("figures", "exploratory", "metabolomics", "obi1", "chromats")


## Get data ------
dat_long_filename <- here(data_dir, "combined_long_dat_scaled.csv")
dat_MF_filename <- here(data_dir, "combined_MF_info.csv")
sample_key <- read_csv(here(meta_data_dir, "sample_key.csv"), show_col_types = FALSE)
dat_long <- read_csv(dat_long_filename, show_col_types = FALSE)
dat_MF <-
  read_csv(dat_MF_filename, show_col_types = FALSE)

## Set mf_oi ------
mf_oi <- c(
    "Homarine",
    "HILICPos_148", # n-methyl glutamine
    "HILICPos_192", # unknown, but likely C7H9NO5 by mass and iso patt
    "HILICPos_29", # unknown
    "HILICPos_153" # n-methyl glutamic acid
    )
bad_mf <- c("HILICNeg_123",  "HILICPos_180")
# Reference spectrum for n-methyl glutamic acid here: https://mona.fiehnlab.ucdavis.edu/spectra/display/CCMSLIB00005720398
# Reference spectrum for glutamine Same https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-BGC_Munich-RP001301

samples_oi <- sample_key %>%
  filter(
    sample_fraction == "particulate",
    chromat == "HILIC",
    sample_set == "OBi1_set2" | replicate_name == "211202_Blk_BlankTube_3",
    treatment == "Homarine" | treatment == "Glucose + NH4" | replicate_name == "211202_Blk_BlankTube_3"
  ) %>%
    mutate(filename = paste0(replicate_name, ".mzML")) %>%
    mutate(treatment = ifelse(is.na(treatment), "Blank", treatment))

# LOOP THROUGH EACH mf_oi_i-------
for (mf_oi_i in mf_oi){
mf_info <- dat_MF %>%
  filter(mass_feature == mf_oi_i)
if (mf_info$z == 1) {
  data_subdir <- "HILIC_positive_mzml"
} else {
  data_subdir <- "HILIC_negative_mzml"
}
all_ms_files <- list.files(here(raw_dat_dir, data_subdir))
msdata_files <- here(raw_dat_dir, data_subdir, all_ms_files[all_ms_files %in% samples_oi$filename])
dda_files <- here(raw_dat_dir, data_subdir, all_ms_files[str_detect(all_ms_files, "_DDA")])
msdata <- grabMSdata(files = msdata_files, grab_what = c("MS1"))

g_ms1 <- plot_EIC(ms1data = msdata$MS1, m_z = mf_info$mz, r_t = mf_info$RT, samples_oi = samples_oi)

msdata_dda <- grabMSdata(files = dda_files, grab_what = c("MS2"))
ms2data <- pull_ms2_data(msdata_dda$MS2, m_z = mf_info$mz, r_t = mf_info$RT)

if (nrow(ms2data) > 0) {
  g_ms2 <- plot_spectrum (ms2data)
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
}

#RDisop for isotope envelopes
#plot(int ~ mz, type = "h", data = ms1_dat, ylab = "Intensity", xlab = "Fragment m/z")
#molecules <- decomposeMass(188.05644, ppm = 5)
