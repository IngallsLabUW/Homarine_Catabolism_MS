library(RaMS)
library(Rdisop)
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
sample_key <- read_csv(here(meta_data_dir, "sample_key.csv"), show_col_types = FALSE)
dat_long <- read_csv(dat_long_filename, show_col_types = FALSE)


## Set mf_oi ------
mf_oi <- c( #THESE ARE ION FORMULAS, NOT MOLECULAR FORMULAS
    #"C7H8NO3",
    #"C6H8NO2",
    #"C6H11NO2",
    #"C7H12NO3",
    #"C6H8NO3",
    #"C6H10NO3",  #insource
    "C6H10NO4",
    "C6H8NO4",
    "C6H8NO3"
    #"C7H9NO5Na"
)
rt_oi <- c(
    #3,
    #9.6,
    #10.6,
    #8,
    #11.9,
    #9.6,
    11.9,
    3,
    11
    #10
)
z_oi <- c(
    #1,
    #1,
    #1,
    #1,
    #1,
    #1,
    1,
    -1,
    -1
    #1
)

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
for (i in 1:length(mf_oi)){
    mf_oi_i <- mf_oi[i]
    z_oi_i <- z_oi[i]
    rt_oi_i <- rt_oi[i]
    mz_oi_i <- getMolecule(mf_oi_i)$isotopes[[1]][1,1]
    if (z_oi_i == 1) {
        data_subdir <- "HILIC_positive_mzml"
    } else {
        data_subdir <- "HILIC_negative_mzml"
    }
    all_ms_files <- list.files(here(raw_dat_dir, data_subdir))
    msdata_files <- here(raw_dat_dir, data_subdir, all_ms_files[all_ms_files %in% samples_oi$filename])
    dda_files <- here(raw_dat_dir, data_subdir, all_ms_files[str_detect(all_ms_files, "(_DDA)|(240319)")])
    msdata <- grabMSdata(files = msdata_files, grab_what = c("MS1"))

    g_ms1 <- plot_EIC(ms1data = msdata$MS1, m_z = mz_oi_i, r_t = rt_oi_i, samples_oi = samples_oi)
    if (mf_oi_i == "C7H9NO5Na"){
        g_ms1 <- plot_EIC(ms1data = msdata$MS1, m_z = mz_oi_i, r_t = rt_oi_i, rt_buffer = 10, samples_oi = samples_oi)
    }

    msdata_dda <- grabMSdata(files = dda_files, grab_what = c("MS2"))
    ms2data <- pull_ms2_data(msdata_dda$MS2, m_z = mz_oi_i, r_t = rt_oi_i, rt_buffer = 2)

    if (nrow(ms2data) > 0) {
        g_ms2 <- plot_spectrum (ms2data)
        g_save <- g_ms1 +
            g_ms2 +
            plot_layout(ncol = 2, widths = c(3, 1)) +
            plot_annotation(
                title = mf_oi_i,
                caption = paste0(
                    "m/z: ",
                    round(mz_oi_i, digits = 4),
                    ", RT: ",
                    round(rt_oi_i, digits = 2),
                    " min, z: ",
                    z_oi_i
                )
            )

        save_width <- 9
    } else {
        g_save <- g_ms1 +
            plot_annotation(
                title = mf_oi_i,
                caption = paste0(
                    "m/z: ",
                    round(mz_oi_i, digits = 4),
                    ", RT: ",
                    round(rt_oi_i, digits = 2),
                    " min, z: ",
                    z_oi_i
                )
            )

        save_width <- 6
    }

    if (z_oi_i == 1) {
        ggsave(here(fig_dir, paste0(mf_oi_i, ".pdf")), g_save, width = save_width, height = 5)
    } else {
        ggsave(here(fig_dir, paste0(mf_oi_i, "_neg", ".pdf")), g_save, width = save_width, height = 5)
    }

}
