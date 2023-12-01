# Libraries ----

library(LNCDR)
library(data.table)
library(dplyr)
library(factoextra)
library(ggplot2)
library(e1071)
library(caret)
attach(mtcars)
library(grid)
library(gridExtra)
library(plotrix)
library(mgcv)
library(readxl)
library(lme4)
library(lubridate)
library(checkmate)
library(lmerTest)
library(tidyr)
library(jtools)

# Define Outlier Function ----
outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

# Load in merge 7T EEG ----
merge7tEEG <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/merge7tEEG.csv')

# Select SNR data ----
SNR <- merge7tEEG[c("labels","Subject", "ERSP", "ITC" ,"BaselinePower","Induced","Evoked", "visitno","age")]

## DLPFC ----

SNR_bothDLPFCs <- rbind(SNR %>% filter(labels %in% c('F3', 'F5', 'F7'))%>% mutate(Region="LDLPFC"), 
                        SNR %>% filter(labels %in% c('F4', 'F6', 'F8'))%>% mutate(Region="RDLPFC"))

SNR_avgDLPFCs <- aggregate(cbind(ERSP, ITC, BaselinePower, Induced, Evoked, age, visitno) ~ Subject + Region, data = SNR_bothDLPFCs, mean) %>% separate(Subject,c("luna","vdate"), remove=F)

### Outlier Detection Subject Level 
SNR_avgDLPFCs_naout <- SNR_avgDLPFCs %>% 
  group_by(Region) %>% 
  mutate(across(c("ERSP", "ITC", "Induced", "Evoked","BaselinePower"), naoutlier)) %>% mutate(Freq ="40Hz")


write.csv(SNR_avgDLPFCs_naout, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectDLPFC_SNRMeasures_20231201.csv')


## Central Parietal for Auditory Response

SNR_CP <- SNR %>% filter(labels %in% c('Cz', 'C2', 'CP2', 'Cz', 'C1', 'CP1', 'CP2', 'CP4', 'P2', 'Pz', 'P1', 'CP1', 'CP3', 'P6', 'CP6', 'P4', 'PO4', 'POz', 'PO3', 'P3', 'PO3', 'P5', 'CP5', 'P8', 'PO8', 'O2', 'Oz', 'O1', 'P7', 'PO7'))
  

SNR_CP <- aggregate(cbind(ERSP, ITC, BaselinePower, Induced, age, visitno) ~ Subject, data = SNR_CP, mean) %>% separate(Subject,c("luna","vdate"), remove=F)


### Outlier Detection Subject Level 
# SNR_CPs_naout <- SNR_CP  %>% 
#   mutate(across(c("ERSP", "ITC", "Induced", "BaselinePower"), naoutlier)) %>% mutate(Freq ="40Hz") %>% separate(Subject,c("luna","vdate"), remove=F)

write.csv(SNR_CP, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectCentralParietal_SNRMeasures_20231114.csv')


# Select FOOOF data ----
fooof <- merge7tEEG[c("labels","Subject","Offset_eyesClosed","Offset_eyesOpen", "Exponent_eyesClosed", "Exponent_eyesOpen","visitno","age")]

## DLPFC ----

fooof_bothDLPFCs <- rbind(fooof %>% filter(labels %in% c('F3', 'F5', 'F7'))%>% mutate(Region="LDLPFC"), 
                          fooof %>% filter(labels %in% c('F4', 'F6', 'F8'))%>% mutate(Region="RDLPFC"))


fooof_avgDLPFCs <- aggregate(cbind(Offset_eyesClosed,Offset_eyesOpen, Exponent_eyesClosed, Exponent_eyesOpen, age, visitno) ~ Subject + Region, data = fooof_bothDLPFCs, mean) %>% separate(Subject,c("luna","vdate"), remove=F)

### Outlier Detection Subject Level 
fooof_avgDLPFCs_naout <- fooof_avgDLPFCs %>% 
  group_by(Region) %>% mutate(across(c("Offset_eyesClosed","Offset_eyesOpen", "Exponent_eyesClosed", "Exponent_eyesOpen"), naoutlier)) %>% ungroup

write.csv(fooof_avgDLPFCs_naout, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectDLPFCfooofMeasures_20231113.csv')




