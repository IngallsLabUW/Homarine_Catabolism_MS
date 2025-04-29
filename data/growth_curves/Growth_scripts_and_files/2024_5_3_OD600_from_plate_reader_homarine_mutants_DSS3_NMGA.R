# Load required libraries
library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
library(dplyr)
library(ggpubr)

# Read the dataset (ONLY ONCE)
OD600 <- fread('2024_4_12_growth_curves_NMGAD_SPO3191_Homi_NMGA.csv', header = TRUE)

# Convert Time to numeric (hours)
OD600[, Time := as.numeric(period_to_seconds(hms(Time))) / 3600]

# Verify conversion
print(str(OD600$Time))  # Should return "num"
print(head(OD600$Time))  # Should show numeric values

# Subset the first 192 rows
OD600_24 <- head(OD600, 192)

# Identify numeric columns except Time
numeric_cols <- names(OD600)[sapply(OD600, is.numeric) & !(names(OD600) %in% c("TimeS", "Time"))]

# Subtract first row from all numeric columns
OD600[, (numeric_cols) := lapply(.SD, function(x) x - x[1]), .SDcols = numeric_cols]

# Plot temperature over time
ggplot() +
  geom_line(data = OD600, aes(x = Time, y = Temp, group = 1, color = "#A54657")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Time (h)", y = "Temperature") +  
  ggtitle('Temperature versus Time') + theme(legend.position = "none")

# ============================ DSS3 Mutants in One Figure =============================

# Define function for generating ggplot graphs
generate_plot <- function(data, wells, title) {
  ggplot() +
    geom_line(data = data, aes(x = Time, y = get(wells[1]), group = 1), color = "#E87D94") +
    geom_line(data = data, aes(x = Time, y = get(wells[2]), group = 1), color = "#E87D94") +
    geom_line(data = data, aes(x = Time, y = get(wells[3]), group = 1), color = "#E87D94") +
    geom_line(data = data, aes(x = Time, y = get(wells[4]), group = 1), color = "#E87D94") +
    geom_line(data = data, aes(x = Time, y = get(wells[5]), group = 1), color = "#E87D94") +
    geom_line(data = data, aes(x = Time, y = get(wells[6]), group = 1), color = "#E87D94") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x = "Time (h)", y = "OD600") +  
    ggtitle(title) + scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") + ylim(0, 0.2)
}

# Create plots for each strain
SPO3191_Hom12 <- generate_plot(OD600, c("B1", "C1", "D1", "E1", "F1", "G1"), "SPO3191 HOM12")
SPO1588_Hom12 <- generate_plot(OD600, c("B2", "C2", "D2", "E2", "F2", "G2"), "SPO1588 HOM12")
SPO1585_Hom12 <- generate_plot(OD600, c("B3", "C3", "D3", "E3", "F3", "G3"), "SPO1585 HOM12")
DSS3_Hom12 <- generate_plot(OD600, c("B4", "C4", "D4", "E4", "F4", "G4"), "DSS3 HOM12")

SPO3191_Hom4 <- generate_plot(OD600, c("B5", "C5", "D5", "E5", "F5", "G5"), "SPO3191 HOM4")
SPO1588_Hom4 <- generate_plot(OD600, c("B6", "C6", "D6", "E6", "F6", "G6"), "SPO1588 HOM4")
SPO1585_Hom4 <- generate_plot(OD600, c("B7", "C7", "D7", "E7", "F7", "G7"), "SPO1585 HOM4")
DSS3_Hom4 <- generate_plot(OD600, c("B8", "C8", "D8", "E8", "F8", "G8"), "DSS3 HOM4")

SPO3191_NMGA <- generate_plot(OD600, c("B9", "C9", "D9", "E9", "F9", "G9"), "SPO3191 NMGA")
SPO1588_NMGA <- generate_plot(OD600, c("B10", "C10", "D10", "E10", "F10", "G10"), "SPO1588 NMGA")
SPO1585_NMGA <- generate_plot(OD600, c("B11", "C11", "D11", "E11", "F11", "G11"), "SPO1585 NMGA")
DSS3_NMGA <- generate_plot(OD600, c("B12", "C12", "D12", "E12", "F12", "G12"), "DSS3 NMGA")

# Arrange all plots in one figure
Gplot <- ggarrange(SPO3191_Hom12, SPO1588_Hom12, SPO1585_Hom12, DSS3_Hom12,
                   SPO3191_Hom4, SPO1588_Hom4, SPO1585_Hom4, DSS3_Hom4,
                   SPO3191_NMGA, SPO1588_NMGA, SPO1585_NMGA, DSS3_NMGA,
                   ncol = 2, nrow = 6)

# Display the arranged plot
print(Gplot)
