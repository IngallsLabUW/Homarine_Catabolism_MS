

#Load Packages
library(tidyverse)
library(rstatix)
#library(ggthemes)



#This script loads in plate reader data and tidys it for analysis and then 
# classifies each isolate as a grower, non-grower, or maybe grower based on several 
# tests and comparisons. Those comparisons are:
# 
# 1. a paired t-test comparing all data from the +homarine treatment to the control, paired by timepoint
#
# 2. a t-test at each timepoint comparing +homarine to the control (2 time points much be significantly different)
#
# 3. A max OD reading of at least 0.1  in the + homarine treatment to be reached 
#
# 4. A max OD in the +homarine treateatment to be at least 5x the max OD in the control
#
# These thresholds are adjustable and set below:

#Growth/No_Growth Thresholds:
paired_ttest_pval_threshold <- 0.05
timepoint_ttest_pval_threshold <- 0.05
num_signif_diff_timepoints_threshold <- 2
minimum_maxOD <- 0.1  ###look at controls, look at average of the growth of the measurment of the controls
maxOD_above_control_threshold <- 5




#Define Inputs:
###import data
growth.data <- read.table('data/isolate_platereader_data/Homarine_growth_test_PS_G4_isolates_17Nov2023.txt',fill=TRUE,skip="#",fileEncoding="UTF-16LE")
isolate.info <- read_csv('data/isolate_platereader_data/Homarine_growth_test_PS_G4_isolate_names_17Nov2023.csv')



# Tidy Growth Data --------------------------------------------------------

###Tidy Growth Data (From Oscar's Code)
plates <- growth.data$V2[grep("Plate:",growth.data$V1)]
remove <- grep("Original|Plate|Temperature|~End",growth.data$V1)
growth.data <- growth.data[-remove,]
row.names(growth.data) <- NULL
##
new.growth.data<-growth.data[,1:12]
for(j in which(growth.data$V1>10)){
  new.growth.data[j,1:12] <- growth.data[j,2:13]}
row.names(new.growth.data) <- NULL

###Convert data to long form 
long.data <- data.frame(plate=rep(plates,each=96),abs=as.numeric(unlist(t(new.growth.data))))
#row.names(long.data) <- NULL 
long.data$position <- rep(1:96,length(plates))

#Indicate condition, time point, and plate no.
long.data$condition <- sapply(stringr::str_split(long.data$plate,"_"),"[",2)
long.data$condition <- ifelse(grepl("Control",long.data$plate),stringr::str_sub(long.data$condition,1,7),stringr::str_sub(long.data$condition,1,8))
long.data$time.point <- sapply(stringr::str_split(long.data$plate,"_"),"[",3)
long.data$plate <- sapply(stringr::str_split(long.data$plate,"_"),"[",1)
long.data$plate <- as.numeric(stringr::str_sub(long.data$plate,-1))


########Combine isolate information with growth data
iso.growth.dat <- left_join(long.data, isolate.info) %>% 
  arrange(isolate,condition,time.point)


####_________Combine growth data with information on time points

#create times dataframe
times <- data.frame(time.point=c("T0","T1","T2","T3","T4"),time=c(0,5.3,20,47,98)) #Time units is hours.

#Merge long data with time points table
dat.platereader <- left_join(iso.growth.dat, times)










# Analyze Growth Data to Identify Growers and Nongrowers ------------------

###Perform paired-sample t-test to identify growth relative to the control:
paired.ttest.output <- dat.platereader %>%
  filter(!isolate %in% c("control", "empty")) %>%
  arrange(isolate, time.point) %>%
  group_by(isolate) %>%
  t_test(., abs~condition, paired = TRUE, detailed = TRUE) %>%
  select(isolate, p) %>%
  rename("paired_ttest_p" = p) %>%
  mutate(paired_ttest_sig = case_when(paired_ttest_p < paired_ttest_pval_threshold ~ 1,
                                      TRUE ~ 0))

##Perform t-test on each time point
timepoint.ttest.output <- dat.platereader %>%
  filter(!isolate %in% c("control", "empty")) %>%
  filter(!time == 0) %>%
  group_by(isolate, time.point) %>%
  t_test(., abs~condition, paired = FALSE, detailed = TRUE) %>%
  select(isolate, time.point, p) %>%
  rename("timepoint_ttest_p" = p) %>%
  group_by(isolate)# %>%
#  mutate(numb_sigdiff_timepoints = sum(timepoint_ttest_p < timepoint_ttest_pval_threshold))


####Identify if points are above or below control
timepoint.relative.OD <- dat.platereader %>%
  filter(!isolate %in% c("control", "empty")) %>%
  filter(!time == 0) %>%
  group_by(isolate, time.point, condition) %>%
  reframe(mean_OD = mean(abs)) %>%
  group_by(isolate, time.point) %>%
  reframe(HC_OD_ratio = mean_OD[condition == "Homarine"]/mean_OD[condition == "Control"]) %>%
  ungroup()

###Combine time point analysis data:
timepoint.analysis <- left_join(timepoint.ttest.output, timepoint.relative.OD) %>%
  mutate(H_sig_greater = case_when(timepoint_ttest_p < timepoint_ttest_pval_threshold & HC_OD_ratio > 1 ~ 1,
                                   TRUE ~ 0)) %>%
  mutate(H_sig_smaller = case_when(timepoint_ttest_p < timepoint_ttest_pval_threshold & HC_OD_ratio < 1 ~ 1,
                                   TRUE ~ 0))


####Combine overall and timepoint analyses:
comb.growth.analysis <- left_join(timepoint.analysis, paired.ttest.output) 



###Classify as grower, non_grower, and maybe_grower 
growth.class <- comb.growth.analysis %>%
  group_by(isolate, paired_ttest_p, paired_ttest_sig) %>%
  reframe(Sum_SigPosGrowth = sum(H_sig_greater),
          Sum_SigNegGrowth = sum(H_sig_smaller)) %>%
  mutate(hom_growth_status = case_when(paired_ttest_sig == 1 & Sum_SigPosGrowth > 0 & Sum_SigNegGrowth == 0 ~ "grower",
                                       paired_ttest_sig == 0 & Sum_SigPosGrowth == 0 & Sum_SigNegGrowth == 0 ~ "non_grower",
                                       paired_ttest_sig == 0 & Sum_SigPosGrowth == 0 & Sum_SigNegGrowth > 0 ~ "non_grower",
                                       paired_ttest_sig == 1 & Sum_SigPosGrowth == 0 & Sum_SigNegGrowth > 0 ~ "non_grower",
                                       TRUE ~ "maybe_grower"))


###export growth classification data:
write_csv(growth.class, file = "data/intermediate/isolategrowth/Isolates_hom_growth_classification.csv")

























































































































