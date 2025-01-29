library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
library(dplyr)
library(ggpubr)
#Table with OD600 values from Gen5 plate reader
data <- fread('2023_8_25_Growth_of_SPO3186_on_Hom_GBT_DMSP_Carn_Nmetglut.csv', header=T)

# Exclude Time and Temp columns from normalization
columns_to_normalize <- setdiff(names(data), c("Time", "Temp"))

# Subtract the minimum value for each column
OD600 <- data %>%
  mutate(across(all_of(columns_to_normalize), ~ . - min(.), .names = "normalized_{col}"))

# Print the first few rows to check
head(OD600)


OD600$TimeS <- period_to_seconds(lubridate::hms(OD600$Time))
OD600$TimeH <- OD600$TimeS/3600
OD600_24 <- head(OD600, 192)


#temperature
ggplot() +
  geom_line(data = OD600, aes(x = Time, y = Temp, 
                              group = 1, color = "#A54657")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="Temperature") +  
  ggtitle('Temperature versus Time') + theme(legend.position = "none") 

#ggsave('2022_12_12_Temperature_versus_Time.pdf')


###########################   These are the different plots for the plate reader  ##############################


#DSS3 in Homarine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B1, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C1, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D1, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in Homarine 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 (wt) in Homarine 4 mM C.pdf')

#DSS3 in Homarine 4 mM C with Kan 50
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E1, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F1, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G1, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 in Homarine 4 mM C with Kan 50.pdf')

#Mutant SPO3186_48115 in Homarine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B2, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C2, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D2, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E2, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F2, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_48115 in Homarine 4 mM C with Kan 50.pdf')


#Mutant SPO3186_46D7 in Homarine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B3, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C3, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D3, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E3, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F3, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_46D7 in Homarine 4 mM C with Kan 50.pdf')

#DSS3 in GBT
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B4, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C4, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D4, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in Glycine betaine 4 mM C') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 (wt) in Glycine betaine 4 mM C.pdf')

#DSS3 in Glycine betaine 4 mM C with Kan 50
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E4, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F4, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G4, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in Glycine betaine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 in Glycine betaine 4 mM C with Kan 50.pdf')

#Mutant SPO3186_48115 in Glycine betaine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B5, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C5, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D5, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E5, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F5, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in Glycine betaine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_48115 in Glycine betaine 4 mM C with Kan 50.pdf')


#Mutant SPO3186_46D7 in Glycine betaine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B6, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C6, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D6, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E6, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F6, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in Glycine betaine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_46D7 in Glycine betaine 4 mM C with Kan 50.pdf')

#DSS3 in DMSP
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B7, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C7, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D7, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in DMSP 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 (wt) in DMSP 4 mM C.pdf')

#DSS3 in DMSP 4 mM C with Kan 50
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E7, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F7, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G7, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 in DMSP 4 mM C with Kan 50.pdf')

#Mutant SPO3186_48115 in DMSP
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B8, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C8, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D8, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E8, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F8, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_48115 in DMSP 4 mM C with Kan 50.pdf')


#Mutant SPO3186_46D7 in DMSP
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B9, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C9, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D9, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E9, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F9, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_46D7 in DMSP 4 mM C with Kan 50.pdf')


#DSS3 in Carnitine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B10, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C10, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D10, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in Carnitine 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 (wt) in Carnitine 4 mM C.pdf')

#DSS3 in Carnitine 4 mM C with Kan 50
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E10, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F10, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G10, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in Carnitine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 in Carnitine 4 mM C with Kan 50.pdf')

#Mutant SPO3186_48115 in Carnitine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B11, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C11, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D11, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E11, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F11, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in Carnitine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_48115 in Carnitine 4 mM C with Kan 50.pdf')


#Mutant SPO3186_46D7 in Carnitine
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B12, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C12, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D12, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E12, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F12, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in Carnitine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_46D7 in Carnitine 4 mM C with Kan 50.pdf')

#DSS3 in N-methyl-glutamic acid
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = H1, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = H2, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G2, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in N-methyl-glutamic acid 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 (wt) in N-methyl-glutamic acid 4 mM C.pdf')

#DSS3 in N-methyl-glutamic acid 4 mM C with Kan 50
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = G3, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = H3, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H4, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in N-methyl-glutamic acid 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('DSS3 in N-methyl-glutamic acid 4 mM C with Kan 50.pdf')

#Mutant SPO3186_48115 in N-methyl-glutamic acid
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = G5, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = G6, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H6, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H7, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = H5, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in N-methyl-glutamic acid 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_48115 in N-methyl-glutamic acid 4 mM C with Kan 50.pdf')


#Mutant SPO3186_46D7 in N-methyl-glutamic acid
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = G8, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = G9, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H8, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H9, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = H10, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in N-methyl-glutamic acid 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

#ggsave('SPO3186_46D7 in N-methyl-glutamic acid 4 mM C with Kan 50.pdf')

#################  SPO3186 mutants in one figure  ############################

DSS3_Hom <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B1, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C1, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D1, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in Homarine 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_Hom_Kan <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E1, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F1, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G1, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_Homarine <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B2, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C2, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D2, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E2, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F2, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_Homarine <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B3, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C3, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D3, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E3, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F3, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_GBT <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B4, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C4, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D4, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in Glycine betaine 4 mM C') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_GBT_Kan <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E4, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F4, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G4, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in Glycine betaine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_GBT <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B5, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C5, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D5, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E5, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F5, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in Glycine betaine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_GBT <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B6, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C6, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D6, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E6, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F6, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in Glycine betaine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_DMSP <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B7, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C7, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D7, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in DMSP 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_DMSP_Kan <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E7, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F7, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G7, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_DMSP <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B8, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C8, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D8, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E8, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F8, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_DMSP <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B9, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C9, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D9, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E9, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F9, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_Carn <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B10, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C10, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D10, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in Carnitine 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_Carn_Kan <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E10, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = F10, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G10, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in Carnitine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_Carn <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B11, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C11, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D11, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E11, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F11, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in Carnitine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_Carn <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B12, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = C12, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = D12, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = E12, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = F12, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in Carnitine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_N_methyl_glutamic_acid <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = H1, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = H2, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = G2, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 (wt) in N-methyl-glutamic acid 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_N_methyl_glutamic_acid_Kan <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = G3, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = H3, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H4, 
                              group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('DSS3 in N-methyl-glutamic acid 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_N_methyl_glutamic_acid <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = G5, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = G6, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H6, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H7, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = H5, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_48115 in N-methyl-glutamic acid 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_N_methyl_glutamic_acid <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = G8, 
                              group = 1, color = "#E87D94")) +
  geom_line(data = OD600, aes(x = TimeH, y = G9, 
                              group=1,), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H8, 
                              group=1), color = "#E87D94") +
  geom_line(data = OD600, aes(x = TimeH, y = H9, 
                              group=1), color = "#E87D94") +  
  geom_line(data = OD600, aes(x = TimeH, y = H10, 
                              group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('SPO3186_46D7 in N-methyl-glutamic acid 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

Gplot <- ggarrange(DSS3_Hom, DSS3_Hom_Kan, SPO3186_48115_Homarine, 
                   SPO3186_46D7_Homarine, DSS3_GBT, 
                   DSS3_GBT_Kan, SPO3186_48115_GBT, 
                   SPO3186_46D7_GBT, DSS3_DMSP, DSS3_DMSP_Kan, 
                   SPO3186_48115_DMSP, SPO3186_46D7_DMSP, DSS3_Carn, DSS3_Carn_Kan, 
                   SPO3186_48115_Carn, SPO3186_46D7_Carn, DSS3_N_methyl_glutamic_acid, 
                   DSS3_N_methyl_glutamic_acid_Kan, SPO3186_48115_N_methyl_glutamic_acid,
                   SPO3186_46D7_N_methyl_glutamic_acid, ncol =2, nrow =10) 

annotate_figure(Gplot,
                top = text_grob("Growth of SPO3186 mutants in Homarine, 
                                Glycine Betaine, DMSP, Carnitine, and N-methyl glutamic acid", 
                                face ='bold', size = 12))
ggsave('2023_8_25_SPO3186_Hom_GBT_DMSP_Carn_Nmet_glut.pdf', plot = Gplot, width = 10, 
       height = 20, units = "in", dpi = 800)

#only Homarine and N-met-glut
Gplot <- ggarrange(DSS3_Hom, DSS3_Hom_Kan, SPO3186_48115_Homarine, 
                   SPO3186_46D7_Homarine, DSS3_N_methyl_glutamic_acid, 
                   DSS3_N_methyl_glutamic_acid_Kan, SPO3186_48115_N_methyl_glutamic_acid,
                   SPO3186_46D7_N_methyl_glutamic_acid, ncol =2, nrow =4) 

annotate_figure(Gplot,
                top = text_grob("Growth of SPO3186 mutants in Homarine, 
                                and N-methyl glutamic acid", 
                                face ='bold', size = 12))
ggsave('2023_8_25_SPO3186_Hom_Nmet_glut.pdf', plot = Gplot, width = 12, 
       height = 10, units = "in", dpi = 600)


###################################### 24 hour ####################################


DSS3_Hom <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B1, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C1, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D1, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 (wt) in Homarine 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_Hom_Kan <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = E1, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = F1, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = G1, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_Homarine <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B2, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C2, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D2, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E2, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F2, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_48115 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_Homarine <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B3, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C3, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D3, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E3, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F3, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_46D7 in Homarine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_GBT <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B4, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C4, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D4, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 (wt) in Glycine betaine 4 mM C') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_GBT_Kan <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = E4, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = F4, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = G4, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 in Glycine betaine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_GBT <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B5, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C5, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D5, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E5, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F5, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_48115 in Glycine betaine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_GBT <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B6, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C6, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D6, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E6, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F6, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_46D7 in Glycine betaine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_DMSP <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B7, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C7, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D7, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 (wt) in DMSP 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_DMSP_Kan <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = E7, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = F7, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = G7, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_DMSP <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B8, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C8, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D8, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E8, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F8, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_48115 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_DMSP <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B9, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C9, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D9, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E9, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F9, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_46D7 in DMSP 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_Carn <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B10, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C10, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D10, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 (wt) in Carnitine 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_Carn_Kan <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = E10, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = F10, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = G10, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 in Carnitine 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_Carn <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B11, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C11, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D11, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E11, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F11, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_48115 in Carnitine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_Carn <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = B12, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = C12, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = D12, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = E12, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = F12, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_46D7 in Carnitine 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_N_methyl_glutamic_acid <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = H1, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = H2, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = G2, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 (wt) in N-methyl-glutamic acid 4 mM C') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

DSS3_N_methyl_glutamic_acid_Kan <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = G3, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = H3, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = H4, 
                                 group=1), color = "#E87D94") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('DSS3 in N-methyl-glutamic acid 4 mM C with Kan 50') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_48115_N_methyl_glutamic_acid <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = G5, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = G6, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = H6, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = H7, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = H5, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_48115 in N-methyl-glutamic acid 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

SPO3186_46D7_N_methyl_glutamic_acid <-
  ggplot() +
  geom_line(data = OD600_24, aes(x = TimeH, y = G8, 
                                 group = 1, color = "#E87D94")) +
  geom_line(data = OD600_24, aes(x = TimeH, y = G9, 
                                 group=1,), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = H8, 
                                 group=1), color = "#E87D94") +
  geom_line(data = OD600_24, aes(x = TimeH, y = H9, 
                                 group=1), color = "#E87D94") +  
  geom_line(data = OD600_24, aes(x = TimeH, y = H10, 
                                 group=1), color = "#E87D94") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600_24") +  
  ggtitle('SPO3186_46D7 in N-methyl-glutamic acid 4 mM C with Kan 50') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")  + ylim(NA,.075) 

Gplot <- ggarrange(DSS3_Hom, DSS3_Hom_Kan, SPO3186_48115_Homarine, 
                   SPO3186_46D7_Homarine, DSS3_GBT, 
                   DSS3_GBT_Kan, SPO3186_48115_GBT, 
                   SPO3186_46D7_GBT, DSS3_DMSP, DSS3_DMSP_Kan, 
                   SPO3186_48115_DMSP, SPO3186_46D7_DMSP, DSS3_Carn, DSS3_Carn_Kan, 
                   SPO3186_48115_Carn, SPO3186_46D7_Carn, DSS3_N_methyl_glutamic_acid, 
                   DSS3_N_methyl_glutamic_acid_Kan, SPO3186_48115_N_methyl_glutamic_acid,
                   SPO3186_46D7_N_methyl_glutamic_acid, ncol =2, nrow =10) 

annotate_figure(Gplot,
                top = text_grob("Growth of SPO3186 mutants in Homarine, 
                                Glycine Betaine, DMSP, Carnitine, and N-methyl glutamic acid", 
                                face ='bold', size = 12))
ggsave('2023_8_25_SPO3186_Hom_GBT_DMSP_Carn_Nmet_glut_24.pdf', plot = Gplot, width = 10, 
       height = 20, units = "in", dpi = 600)


