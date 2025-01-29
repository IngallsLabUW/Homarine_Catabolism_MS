library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
library(dplyr)
library(ggpubr)

#Table with OD600 values from Gen5 plate reader
OD600 <- fread('2022_12_9_Plate_1_Obi1_DSS3_raw_reads.csv', header=T)
OD600$TimeS <- period_to_seconds(lubridate::hms(OD600$Time))
OD600$TimeH <- OD600$TimeS/3600

#temperature
ggplot() +
  geom_line(data = OD600, aes(x = Time, y = Temp, 
                              group = 1, color = "#A54657")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="Temperature") +  
  ggtitle('Temperature versus Time') + theme(legend.position = "none") 


#These are the different plots for the plate reader


#Obi1 in G1
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B2, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C2, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D2, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, Obi1, Group G1: Glucose + TM + PO4 + V + NH4') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 

#Obi1 in G2
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B3, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C3, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D3, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, Obi1, Group G2: Glucose + TM + PO4 + V') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 

#Obi1 in G3
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B4, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C4, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D4, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, Obi1, Group G3: Glucose + TM + PO4') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 

#Obi1 in H1
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B5, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C5, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D5, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, Obi1, Group H1: Homarine + TM + PO4 + V + NH4') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 

#Obi1 in H2
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B6, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C6, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D6, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, Obi1, Group H2: Glucose + TM + PO4 + V') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 

#Obi1 in H3
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B7, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C7, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D7, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, Obi1, Group H3: Glucose + TM + PO4') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 

#DSS3 in G1
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E2, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F2, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G2, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, DSS3, Group G1: Glucose + TM + PO4 + V + NH4 ') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 


#DSS3 in G2
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E3, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F3, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G3, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, DSS3, Group G2: Glucose + TM + PO4 + V') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 


#DSS3 in G3
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E4, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F4, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G4, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, DSS3, Group G3: Glucose + TM + PO4') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 


#DSS3 in H1
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E5, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F5, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G5, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, DSS3, Group H1: Homarine + TM + PO4 + V + NH4 ') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))  + theme(legend.position = "none") 

#DSS3 in H2
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E6, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F6, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G6, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, DSS3, Group H2: Glucose + TM + PO4 + V ') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))  + theme(legend.position = "none") 

#DSS3 in H3
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E7, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F7, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G7, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('OD versus Time, DSS3, Group H3: Glucose + TM + PO4 ') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))  + theme(legend.position = "none") 


#################  Obi1 in each growth condition in one figure  ############################

Obi1G1 <- ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B2, 
                                       group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C2, 
                                       group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D2, 
                                       group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y="OD600") +  
  ggtitle('Group G1: Glucose + TM + PO4 + V + NH4') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.title = element_blank()) + theme(legend.position = "none")

Obi1G2 <- ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B3, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C3, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D3, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y="OD600") +  
  ggtitle('Group G2: Glucose + TM + PO4 + V') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")

Obi1G3 <- ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B4, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C4, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D4, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('Group G3: Glucose + TM + PO4') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")

Obi1H1 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B5, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C5, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D5, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y=" ") +  
  ggtitle('Group H1: Homarine + TM + PO4 + V + NH4') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")

Obi1H2 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B6, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C6, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D6, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y=" ") +  
  ggtitle('Group H2: Homarine + TM + PO4 + V') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")

Obi1H3 <-
  ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = B7, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = C7, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = D7, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y=" ") +  
  ggtitle('Group H3: Homarine + TM + PO4') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")

Gplot <- ggarrange(Obi1G1, Obi1H1, Obi1G2, Obi1H2, Obi1G3,   
                   Obi1H3, ncol =2, nrow =3) 

annotate_figure(Gplot,
                top = text_grob("Growth of Obi1 in Glucose and Homarine 12 mM C, 48 h", 
                                face ='bold', size = 14))



#################  DSS3 in each growth condition in one figure  ############################


DSS3G1 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E2, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F2, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G2, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y="OD600") +  
  ggtitle('Group G1: Glucose + TM + PO4 + V + NH4 ') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 

DSS3G2 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E3, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F3, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G3, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y="OD600") +  
  ggtitle('Group G2: Glucose + TM + PO4 + V') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 


DSS3G3 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E4, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F4, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G4, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y="OD600") +  
  ggtitle('Group G3: Glucose + TM + PO4 ') + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none") 


DSS3H1 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E5, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F5, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G5, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y=" ") +  
  ggtitle('Group H1: Homarine + TM + PO4 + V + NH4 ') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))  + theme(legend.position = "none") 

DSS3H2 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E6, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F6, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G6, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=" ", y=" ") +  
  ggtitle('Group H2: Homarine + TM + PO4 + V ') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))  + theme(legend.position = "none") 

DSS3H3 <-
ggplot() +
  geom_line(data = OD600, aes(x = TimeH, y = E7, 
                              group = 1, color = "#A54657")) +
  geom_line(data = OD600, aes(x = TimeH, y = F7, 
                              group=1,), color = "#582630") +
  geom_line(data = OD600, aes(x = TimeH, y = G7, 
                              group=1), color = "#F7EE7F") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Time (h)", y=" ") +  
  ggtitle('Group H3: Homarine + TM + PO4 ') + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))  + theme(legend.position = "none") 


Gplot <- ggarrange(DSS3G1, DSS3H1, DSS3G2, DSS3H2, DSS3G3,  
                    DSS3H3, ncol =2, nrow =3) 

annotate_figure(Gplot,
                top = text_grob("Growth of DSS3 in Glucose and Homarine 12 mM C, 48 h", 
                                face ='bold', size = 14))



