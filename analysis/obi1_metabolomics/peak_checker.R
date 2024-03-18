library(RaMS)
library(here)



dat_long_filename <- here(data_dir, "combined_long_dat_scaled.csv")
dat_MF_filename <- here(data_dir, "combined_MF_info.csv")

dat_long <- read_csv(dat_long_filename, show_col_types = FALSE)
dat_MF <-
  read_csv(dat_MF_filename, show_col_types = FALSE)

# Positive files
dat_loc <- here("data", "raw", "metabolomics", "obi1", "particulate", "mzml", "HILIC_negative_mzml")
files_oi <- paste0(c(
  # Homarine
  "211202_Smp_OBi1_4",
  "211202_Smp_OBi1_5",
  "211202_Smp_OBi1_6",

  # Control
  "211202_Smp_OBi1_1",
  "211202_Smp_OBi1_2",
  "211202_Smp_OBi1_3",

  # Homarine+gluce
  "211202_Smp_OBi1_7",
  "211202_Smp_OBi1_8",
  "211202_Smp_OBi1_9",

  # Blank
  "211202_Smp_MQBlk2_A",

  # Standard
  "211202_Std_4uMStdsMix1InH2O_1",
  "211202_Std_4uMStdsMix2InH2O_1",

  # Standard_matrix
  "211202_Std_4uMStdsMix1InMatrix_1",
  "211202_Std_4uMStdsMix2InMatrix_1"
), ".mzML")

msdata_files <- list.files(dat_loc)
msdata_files <- here(dat_loc, msdata_files[msdata_files %in% files_oi])

msdata <- grabMSdata(files = msdata_files, grab_what = c("MS1"))

mass <- 188.05644 # unknown HILICPos_192 (11.85 min) (7-8 c by iso windwo, likely C7H10NO5 (as ion))
mass <- 109.02882 # unknown HILICpos_29 (10.68 minuts) (5-6 carbons by iso env, likely C6H5O2)
# mass = 285.15558    # This one goes down
# mass = 148.06021    # glutamic acid
# mass = 161.09213    # GORGEOUS! HILIC pos 148 [probably n-methyl glutamine!] by iso envelope
# mass = 162.07591    # n-methyl glutamic acid :)



mass <- 180.0294 # unknown HILICNeg_123 @ 5.5
mass <- 136.03932 # Homarine
# mass = 126.90396    # Iodine, pretty sure.

peak_data <- msdata$MS1[mz %between% pmppm(mass, 20)]

ggplot(peak_data) +
  geom_line(aes(x = rt, y = int, color = filename)) +
  facet_wrap(~filename, ncol = 1)


ms1_dat <- msdata$MS1[rt %between% c(11.8, 11.9) & mz %between% c(188, 200)] %>%
  filter(filename == "211202_Smp_OBi1_5.mzML") %>%
  mutate(mz = round(mz, digits = 3)) %>%
  group_by(mz) %>%
  summarise(int = sum(int))

plot(int ~ mz, type = "h", data = ms1_dat, ylab = "Intensity", xlab = "Fragment m/z")
molecules <- decomposeMass(188.05644, ppm = 5)
