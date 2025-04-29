library(tidyverse)
library(ggpubr)
library(data.table)
library(lubridate)

# Load OD600 values from Gen5 plate reader
OD600 <- fread('2022_12_12_Plate_2_DSS3_homarine_mutants.csv', header=T)

# Convert time to seconds and hours
OD600$TimeS <- period_to_seconds(lubridate::hms(OD600$Time))
OD600$TimeH <- OD600$TimeS / 3600

# Define the blank/control well
control_col <- "B1"  # Background control well

# List of sample wells (ensuring correct placements based on your script)
sample_wells <- c("B2", "B3", "B4", "B5", "B6",  # DSS3 in H3
                  "C2", "D2", "E2", "F2", "G2",  # Mutant 1 in H3
                  "C3", "D3", "E3", "F3", "G3",  # Mutant 2 in H3
                  "C4", "D4", "E4", "F4", "G4",  # Mutant 3 in H3
                  "C5", "D5", "E5", "F5", "G5",  # Mutant 4 in H3
                  "C6", "D6", "E6", "F6", "G6",  # Mutant 5 in H3
                  "C7", "D7", "E7", "F7", "G7",  # Mutant 6 in H3
                  "C8", "D8", "E8", "F8", "G8",  # Mutant 7 in H3
                  "C9", "D9", "E9", "F9", "G9",  # Mutant 8 in H3
                  "C10", "D10", "E10", "F10", "G10",  # Mutant 9 in H3
                  "C11", "D11", "E11", "F11", "G11",  # Mutant 10 in H3
                  "B12", "C12", "D12", "E12", "F12")  # DSS3_Kan

# Ensure the control column exists
if (control_col %in% colnames(OD600)) {
  # Subtract the background OD600 from each sample well, ensuring non-negative values
  OD600[, (sample_wells) := lapply(.SD, function(x) pmax(x - get(control_col), 0)), .SDcols = sample_wells]
} else {
  stop("Control well (B1) not found in the dataset.")
}

# Function to create growth curve plot for each group
create_group_plot <- function(data, wells, title) {
  ggplot(data = data) +
    geom_line(aes(x = TimeH, y = get(wells[1])), color = "#E87D94", group = 1) +
    geom_line(aes(x = TimeH, y = get(wells[2])), color = "#E87D94", group = 1) +
    geom_line(aes(x = TimeH, y = get(wells[3])), color = "#E87D94", group = 1) +
    geom_line(aes(x = TimeH, y = get(wells[4])), color = "#E87D94", group = 1) +
    geom_line(aes(x = TimeH, y = get(wells[5])), color = "#E87D94", group = 1) +
    labs(x = "Time (h)", y = "OD600 (Blank Corrected)", title = title) +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"))
}

# Define groups for plotting
group_titles <- c(
  "DSS3 (Group H3: Homarine + TM + PO4)", "SPO3186_48115 (Transporter)", 
  "SPO3186_4607 (Transporter)", "SPO3187_108C1 (Sulfolactate dehydrogenase)",
  "SPO3187_37010 (Sulfolactate dehydrogenase)", "SPO3189_17N5 (Phenol hydroxylase)",
  "SPO3189_43B4 (Phenol hydroxylase)", "SPO3190 (Putative amidohydrolase)",
  "SPO3192_23A2 (Aminoaldehyde dehydrogenase)", "SPO3192_117L10 (Aminoaldehyde dehydrogenase)",
  "SPO3185_110715 (Transcriptional regulator)", "DSS3_Kan (Wild type + Antibiotic)"
)

# Create plots for each group
group_wells <- list(
  c("B2", "B3", "B4", "B5", "B6"),
  c("C2", "D2", "E2", "F2", "G2"),
  c("C3", "D3", "E3", "F3", "G3"),
  c("C4", "D4", "E4", "F4", "G4"),
  c("C5", "D5", "E5", "F5", "G5"),
  c("C6", "D6", "E6", "F6", "G6"),
  c("C7", "D7", "E7", "F7", "G7"),
  c("C8", "D8", "E8", "F8", "G8"),
  c("C9", "D9", "E9", "F9", "G9"),
  c("C10", "D10", "E10", "F10", "G10"),
  c("C11", "D11", "E11", "F11", "G11"),
  c("B12", "C12", "D12", "E12", "F12")
)

group_plots <- map2(group_wells, group_titles, ~create_group_plot(OD600, .x, .y))

# Arrange all plots in one figure
mutant_plot_grid <- ggarrange(plotlist = group_plots, ncol = 2, nrow = 6)

# Annotate the figure
annotate_figure(mutant_plot_grid,
                top = text_grob("Growth of DSS3 mutants in Homarine 12 mM C (Blank Corrected)", 
                                face = "bold", size = 18))
# 
# # Save the final figure
# ggsave('Growth_of_DSS3_mutants_Homarine_12mM_C_Corrected.pdf', plot = mutant_plot_grid, width = 12, height = 15)
