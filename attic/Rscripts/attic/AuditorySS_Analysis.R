
# Auditory Steady State


library(LNCDR)
library(data.table)
require(dplyr)
library(factoextra)
library(tidyverse)
require(knitr)
library(ggplot2)
library(e1071)
library(caret)
library(sjPlot)
library(directlabels)
attach(mtcars)
library(grid)
library(gridExtra)
library(cowplot)
library(plotrix)
library(mgcv)
library(plotly)
library(readxl)
library(interactions)
library(eegUtils)
library(lme4)
library("ggpubr")
library(jtools)

## load in age file and evoked/spontaneous file

Load40hz <- function() {
agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
agefile$inverseAge <- 1/agefile$age

Amplitudes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/AudSS_EvokedSpontaneousActivity_method3_40hzBandpassed_20230206.csv')
Amplitudes_AdditionalSubjects <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/AudSS_EvokedSpontaneousActivity_method3_40hzBandpassed_AdditionalSubjects_20230206.csv')

CombineSubjectAmplitudes <- rbind(Amplitudes, Amplitudes_AdditionalSubjects)

Amplitudes_age <- merge(agefile, CombineSubjectAmplitudes, by = c("Subject"), all = T)

## Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 

outliers <- function(x) {
  
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
  
  
}

Amplitudes_age_new <- Amplitudes_age

cols = names(Amplitudes_age[7:8])
for (col in cols) {
  
  Amplitudes_age_grouped <- Amplitudes_age %>% group_by(Subject)
  
  indx <- outliers(Amplitudes_age_grouped[[col]])
  
  Amplitudes_age_new[[col]] <- Map(replace, Amplitudes_age_new[[col]], indx, NA)
  Amplitudes_age_new[[col]] <- as.numeric(Amplitudes_age_new[[col]])
  
}  
## Add behavior
sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/03.BehaviorAcrossAge.R", envir = knitr::knit_global(), chdir = TRUE)

Behavior <- Behavior_Sublevel_Maria()
Behavior <- Behavior %>% separate(Subject,c("luna","vdate"), remove=F)


Amplitudes_age_new <- merge(Amplitudes_age_new, Behavior, by = c("Subject", "age", "inverseAge"))

# log the evoked activity to create a more normal distrubution
Amplitudes_age_new$logAmp <- log10(Amplitudes_age_new$Amplitude)

#write.csv(Amplitudes_age_new, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/EvokedSpontaneousBehavior_method3_40hz.csv')

## Map evoked power on brain 
# Load channel locations
chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
chanLocs$Channel <- chanLocs$labels

Amplitudes_channels <- merge(Amplitudes_age_new, chanLocs, by = "Channel")
Amplitudes_channels$zscoredDifference <- zscore(Amplitudes_channels$logAmp) - zscore(Amplitudes_channels$Spontaneous)
Amplitudes_channels$Difference <- (Amplitudes_channels$logAmp) - (Amplitudes_channels$Spontaneous)


SNRCoef <- data.frame()
for (chan in unique(Amplitudes_channels$Channel)) {
  ageRes <- lm(z~invage , data = Amplitudes_channels %>% filter(Channel == chan) %>% mutate(z = scale(Difference)[,1], invage=1/age))   
  SNRCoef <- rbind(SNRCoef, data.frame(chan,ageRes$coefficients[2]))
}

SNRCoef$labels <- SNRCoef$chan
SNRCoef <- merge(SNRCoef, chanLocs, by = "labels")


ggplot(SNRCoef, aes(x = -Y, y = X, fill = ageRes.coefficients.2., z = ageRes.coefficients.2., label = Channel)) + geom_topo(chan_markers = "text") + scale_fill_distiller(palette = "RdBu")

lunaize(ggplot(Amplitudes_channels %>% mutate(isDLPFC = labels %in% c('F3','F5','F7','F4','F6','F8')), aes(x = age, y = Difference, color = isDLPFC)) + geom_smooth(aes(group = labels), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.01,size=1)) + scale_color_manual(values=c("lightgray", "blue4")) +ylab("SNR")


## average the evoked and spontaneous activity together to get a "whole brain" measure 
WholeBrainAmplitudes <- aggregate(.~ Subject+sex, Amplitudes_age_new , mean)
WholeBrainAmplitudes$Region <- "Whole Brain"
WholeBrainAmplitudes$Difference <- (WholeBrainAmplitudes$logAmp) - (WholeBrainAmplitudes$Spontaneous)

write.csv(WholeBrainAmplitudes, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/wholeBrainSNR.csv')

## Regional Data

regionDF <- Amplitudes_channels %>% filter(., Channel == 'F3' | Channel == 'F5') %>% group_by(Subject) %>% 
  summarise(Amplitude = mean(Amplitude, na.rm=T), 
            Spontaneous = mean(Spontaneous, na.rm=T), 
            logAmp = mean(logAmp, na.rm=T),
            zscoredDifference = mean(zscoredDifference, na.rm=T),
            Difference = mean(Difference, na.rm=T)) %>%
  mutate(Region = "LDLPFC") %>% 
  rbind(.,Amplitudes_channels %>% filter(.,  Channel == 'F4' | Channel == 'F6') %>% group_by(Subject) %>% 
          summarise(Amplitude = mean(Amplitude, na.rm=T), 
                    Spontaneous = mean(Spontaneous, na.rm=T), 
                    logAmp = mean(logAmp, na.rm=T),
                    zscoredDifference = mean(zscoredDifference, na.rm=T),
                    Difference = mean(Difference, na.rm=T)) %>%
          mutate(Region = "RDLPFC"))  %>%
  rbind(.,Amplitudes_channels %>% filter(., str_detect(Channel, "F")) %>% group_by(Subject) %>% 
          summarise(Amplitude = mean(Amplitude, na.rm=T), 
                    Spontaneous = mean(Spontaneous, na.rm=T), 
                    logAmp = mean(logAmp, na.rm=T),
                    zscoredDifference = mean(zscoredDifference, na.rm=T),
                    Difference = mean(Difference, na.rm=T)) %>%
          mutate(Region = "MPFC"))

regionDF_new_age <- regionDF %>% merge(., agefile, by = "Subject")

ASSR40 <- regionDF_new_age
ASSR40$Frequency <- '40 Hz'

ASSR40 <- ASSR40 %>% separate(Subject,c("luna","vdate"), remove=F)
write.csv(ASSR40, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/EvokedSpontaneousBehavior_40Hz_20221010.csv')
write.csv(regionDF, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_40hz_SNR_bandpassed.csv')



}

Load30hz <- function() {
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  agefile$inverseAge <- 1/agefile$age
  
  Amplitudes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/AudSS_EvokedSpontaneousActivity_method3_30hz.csv')
  
  Amplitudes_age <- merge(agefile, Amplitudes, by = "Subject")
  Amplitudes_age <- subset(Amplitudes_age, visitno == 1)
  
  ## Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  Amplitudes_age_new <- Amplitudes_age
  
  cols = names(Amplitudes_age[8:9])
  for (col in cols) {
    
    Amplitudes_age_grouped <- Amplitudes_age %>% group_by(Subject)
    
    indx <- outliers(Amplitudes_age_grouped[[col]])
    
    Amplitudes_age_new[[col]] <- Map(replace, Amplitudes_age_new[[col]], indx, NA)
    Amplitudes_age_new[[col]] <- as.numeric(Amplitudes_age_new[[col]])
    
  }  
  ## Add behavior
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/03.BehaviorAcrossAge.R", envir = knitr::knit_global(), chdir = TRUE)
  
  Behavior <- Behavior_Sublevel_Maria()
  
  Amplitudes_age_new <- merge(Amplitudes_age_new, Behavior, by = c("Subject", "age", "inverseAge"))
  
  # log the evoked activity to create a more normal distrubution
  Amplitudes_age_new$logAmp <- log10(Amplitudes_age_new$Amplitude)
  
  write.csv(Amplitudes_age_new, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/EvokedSpontaneousBehavior_30hz.csv')
  
  
  ## Regional Data
  #stat_smooth(method = "lm", formula = 'y ~ I(1/x)')
  
  
  # one channel
  AmplitudesC4 <- filter(Amplitudes_age_new, Amplitudes_age_new$Channel == 'C4')
  AmplitudesC4 <- aggregate(.~Subject, AmplitudesC4, mean)
  AmplitudesC4$Region <- "C4"
  
  lunaize(ggplot(data = AmplitudesC4, 
                 aes(x = age, y = logAmp)) + geom_point() + stat_smooth(method='gam', alpha = 0.1))
  
  
  gam.model <- gam(logAmp ~ s(age, k=4),
                   data = AmplitudesC4 , 
                   random = ~(1|subj))
  summary(gam.model)
  
  
  # one channel
  AmplitudesC3 <- filter(Amplitudes_age_new, Amplitudes_age_new$Channel == 'C3')
  AmplitudesC3 <- aggregate(.~Subject, AmplitudesC3, mean)
  AmplitudesC3$Region <- "C3"
  
  
  lunaize(ggplot(data = AmplitudesC3, 
                 aes(x = age, y = logAmp)) + geom_point() + stat_smooth(method='gam', alpha = 0.1))
  
  
  gam.model <- gam(logAmp ~ s(age, k=4),
                   data = AmplitudesC3 , 
                   random = ~(1|subj))
  summary(gam.model)
  
  # Left DLPFC
  Amplitudes_LDLPFC <- filter(Amplitudes_age_new, Amplitudes_age_new$Channel == 'F3')
  Amplitudes_LDLPFC <- aggregate(.~Subject, Amplitudes_LDLPFC, mean)
  Amplitudes_LDLPFC$Region <- "Left DLPFC"
  Amplitudes_LDLPFC <- filter(Amplitudes_LDLPFC, Amplitudes_LDLPFC$logAmp > -3) 
  
  # Right DLPFC
  Amplitudes_RDLPFC <- filter(Amplitudes_age_new, Amplitudes_age_new$Channel == 'F4')
  Amplitudes_RDLPFC <- aggregate(.~Subject, Amplitudes_RDLPFC, mean)
  Amplitudes_RDLPFC$Region <- "Right DLPFC"
  Amplitudes_RDLPFC <- filter(Amplitudes_RDLPFC, Amplitudes_RDLPFC$logAmp > -3) 
  
  
  # frontal
  Amplitudes_frontal <- filter(Amplitudes_age_new, str_detect(Amplitudes_age_new$Channel, "F"))
  Amplitudes_frontal<- aggregate(.~Subject, Amplitudes_frontal, mean)
  Amplitudes_frontal$Region <- "Frontal"
  
  # occipical
  Amplitudes_occipital <- filter(Amplitudes_age_new, str_detect(Amplitudes_age_new$Channel, "O"))
  Amplitudes_occipital<- aggregate(.~Subject, Amplitudes_occipital, mean)
  Amplitudes_occipital$Region <- "Occipital"
  
  # parietal
  Amplitudes_parietal <- filter(Amplitudes_age_new, str_detect(Amplitudes_age_new$Channel, "P"))
  Amplitudes_parietal<- aggregate(.~Subject, Amplitudes_parietal, mean)
  Amplitudes_parietal$Region <- "Parietal"
  
  Amplitudes_allRegions <- rbind(Amplitudes_LDLPFC, Amplitudes_RDLPFC) %>% rbind (., Amplitudes_frontal) %>% rbind(., Amplitudes_occipital) %>% rbind (., Amplitudes_parietal)%>% rbind (., AmplitudesC4)  %>% rbind (., AmplitudesC3)
  
  ## Compute SNR by z scoring and then taking the difference
  Amplitudes_allRegions$zscoredDifference <- zscore(Amplitudes_allRegions$logAmp) - zscore(Amplitudes_allRegions$Spontaneous)
  
  ASSR30 <- Amplitudes_allRegions
  ASSR30$Frequency <- '30 Hz'
  
}

Load20hz <- function() {
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  agefile$inverseAge <- 1/agefile$age
  
  Amplitudes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/AudSS_EvokedSpontaneousActivity_method3_20hzBandpassed_20230206.csv')
  Amplitudes_AdditionalSubjects <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/AudSS_EvokedSpontaneousActivity_method3_20hzBandpassed_AdditionalSubjects_20230206.csv')
  
  CombineSubjectAmplitudes <- rbind(Amplitudes, Amplitudes_AdditionalSubjects)
  
  Amplitudes_age <- merge(agefile, CombineSubjectAmplitudes, by = c("Subject"), all = T)
  
  ## Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  Amplitudes_age_new <- Amplitudes_age
  
  cols = names(Amplitudes_age[7:8])
  for (col in cols) {
    
    Amplitudes_age_grouped <- Amplitudes_age %>% group_by(Subject)
    
    indx <- outliers(Amplitudes_age_grouped[[col]])
    
    Amplitudes_age_new[[col]] <- Map(replace, Amplitudes_age_new[[col]], indx, NA)
    Amplitudes_age_new[[col]] <- as.numeric(Amplitudes_age_new[[col]])
    
  }  
  
  # log the evoked activity to create a more normal distrubution
  Amplitudes_age_new$logAmp <- log10(Amplitudes_age_new$Amplitude)
  
  #write.csv(Amplitudes_age_new, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/EvokedSpontaneousBehavior_method3_20hz.csv')
  
  ## Map evoked power on brain 
  # Load channel locations
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs$Channel <- chanLocs$labels
  
  Amplitudes_channels <- merge(Amplitudes_age_new, chanLocs, by = "Channel")
  Amplitudes_channels$zscoredDifference <- zscore(Amplitudes_channels$logAmp) - zscore(Amplitudes_channels$Spontaneous)
  Amplitudes_channels$Difference <- (Amplitudes_channels$logAmp) - (Amplitudes_channels$Spontaneous)
  
  
  ggplot(Amplitudes_channels, aes(x = -Y, y = X, fill = Difference, z = Difference, label = Channel)) + geom_topo(chan_markers = "text") + scale_fill_distiller(palette = "RdBu") + facet_wrap(.~Adult)
  
  ## average the evoked and spontaneous activity together to get a "whole brain" measure 
  WholeBrainAmplitudes <- aggregate(.~ Subject+sex, Amplitudes_age_new , mean)
  WholeBrainAmplitudes$Region <- "Whole Brain"
  WholeBrainAmplitudes$Difference <- (WholeBrainAmplitudes$logAmp) - (WholeBrainAmplitudes$Spontaneous)
  
  write.csv(WholeBrainAmplitudes, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/wholeBrainSNR.csv')
  
  ## Regional Data
  
  regionDF <- Amplitudes_channels %>% filter(., Channel == 'F3' | Channel == 'F5') %>% group_by(Subject) %>% 
    summarise(Amplitude = mean(Amplitude, na.rm=T), 
              Spontaneous = mean(Spontaneous, na.rm=T), 
              logAmp = mean(logAmp, na.rm=T),
              zscoredDifference = mean(zscoredDifference, na.rm=T),
              Difference = mean(Difference, na.rm=T)) %>%
    mutate(Region = "LDLPFC") %>% 
    rbind(.,Amplitudes_channels %>% filter(.,  Channel == 'F4' | Channel == 'F6') %>% group_by(Subject) %>% 
            summarise(Amplitude = mean(Amplitude, na.rm=T), 
                      Spontaneous = mean(Spontaneous, na.rm=T), 
                      logAmp = mean(logAmp, na.rm=T),
                      zscoredDifference = mean(zscoredDifference, na.rm=T),
                      Difference = mean(Difference, na.rm=T)) %>%
            mutate(Region = "RDLPFC"))  %>%
    rbind(.,Amplitudes_channels %>% filter(., str_detect(Channel, "F")) %>% group_by(Subject) %>% 
            summarise(Amplitude = mean(Amplitude, na.rm=T), 
                      Spontaneous = mean(Spontaneous, na.rm=T), 
                      logAmp = mean(logAmp, na.rm=T),
                      zscoredDifference = mean(zscoredDifference, na.rm=T),
                      Difference = mean(Difference, na.rm=T)) %>%
            mutate(Region = "MPFC"))
  
  regionDF_new_age <- regionDF %>% merge(., agefile, by = "Subject")
  
  ASSR20 <- regionDF_new_age
  ASSR20$Frequency <- '20 Hz'
  
  ASSR20 <- ASSR20 %>% separate(Subject,c("luna","vdate"), remove=F)
  write.csv(ASSR20, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/EvokedSpontaneousBehavior_20Hz_20230207.csv')
  write.csv(regionDF, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_20hz_SNR_bandpassed.csv')
  
}

#combine all ASSR dataframes
allASSRfreqs <- rbind(ASSR40, ASSR30) %>% rbind(., ASSR20)
allASSRfreqs <- allASSRfreqs %>% separate(Subject,c("luna","vdate"), remove=F)


EvokedActivity <- function() {
  
regionDF_40 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_40hz_SNR_bandpassed.csv') %>% dplyr::select(2:8)
regionDF_20 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_20hz_SNR_bandpassed.csv') %>% dplyr::select(2:8)

regionDF_20$Frequency <- '20 Hz'
regionDF_40$Frequency <- '40 Hz'

rbind40and20 <- rbind(regionDF_20, regionDF_40) %>% merge(., agefile, by = "Subject")%>% separate(Subject,c("luna","vdate"), remove=F)
write.csv(rbind40and20, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/EvokedSpontaneousBehavior_40and20Hz_20230209.csv')


# evoked
lunaize(ggplot(data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = logAmp, by = luna))+ geom_line(aes(group = interaction(luna, Region, Frequency)), alpha = 0.1) + geom_point(aes(color=Region),alpha=.5)+ geom_smooth(aes(group = Frequency, linetype = Frequency), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1) + scale_color_manual(values=c("gold3", "blue4")) + ylab("Evoked") +xlab("Age"))
        
gam.model <-  gam(logAmp ~ s(age, k = 3) + Region + Frequency, data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
summary(gam.model)

rbind40and20$oFreq <- ordered(rbind40and20$Frequency, levels = c("40 Hz","20 Hz")) # eyes open will be the reference group
model_formula <- as.formula("logAmp ~ oFreq + s(age, k = 4, fx = T) + s(age, by = oFreq, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)
        

# spontaneous vs age for 40 and 20
lunaize(ggplot(data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Spontaneous)) + geom_point(aes(color=Region),alpha=.8)+ geom_smooth(aes(group = Frequency, linetype = Frequency), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + scale_color_manual(values=c("gold3", "blue4")) + ylab("Spontaneous") +xlab("Age"))

gam.model <-  gam(Spontaneous ~ s(age, k = 3) + Region + Frequency, data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
summary(gam.model)

rbind40and20$oFreq <- ordered(rbind40and20$Frequency, levels = c("40 Hz","20 Hz")) # eyes open will be the reference group
model_formula <- as.formula("Spontaneous ~ oFreq + s(age, k = 4, fx = T) + s(age, by = oFreq, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)
      


# SNR vs age for 40 and 20
lunaize(ggplot(data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Difference)) + geom_point(aes(color=Region),alpha=.8)+ geom_smooth(aes(group = Frequency, linetype = Frequency), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + scale_color_manual(values=c("gold3", "blue4")) + ylab("SNR") +xlab("Age"))

gam.model <-  gam(Difference ~ s(age, k = 3) + Region + Frequency, data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
summary(gam.model)

rbind40and20$oFreq <- ordered(rbind40and20$Frequency, levels = c("40 Hz","20 Hz")) # eyes open will be the reference group
model_formula <- as.formula("Difference ~ oFreq + s(age, k = 4, fx = T) + s(age, by = oFreq, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = rbind40and20 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)


## CREATE 40-20 MEASURES

combine40and20 <- merge(regionDF_20, regionDF_40, by = c("Subject", "Region"), suffixes = c("_20", "_40"))

subtract20from40 <- combine40and20 %>% group_by(Subject, Region) %>%
  summarise(logAmp = logAmp_40 - logAmp_20, 
            Spontaneous = Spontaneous_40 - Spontaneous_20, 
            Difference = Difference_40 - Difference_20,
            zscoredDifference = zscoredDifference_40 - zscoredDifference_20)

regionDF_new_age <- subtract20from40 %>% merge(., agefile, by = "Subject")
regionDF_new_age <- regionDF_new_age %>% separate(Subject,c("luna","vdate"), remove=F)
write.csv(regionDF_new_age, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_40minus20_SNR_bandpassed.csv')


Amplitudes_age_new <- regionDF_new_age

cols = names(regionDF_new_age[5:8])
for (col in cols) {
  
  regionDF_new_age_grouped <- regionDF_new_age %>% group_by(Subject, Region)
  
  indx <- outliers(regionDF_new_age_grouped[[col]])
  
  Amplitudes_age_new[[col]] <- Map(replace, Amplitudes_age_new[[col]], indx, NA)
  Amplitudes_age_new[[col]] <- as.numeric(Amplitudes_age_new[[col]])
  
}  
regionDF_new_age <- Amplitudes_age_new

## EVOKED VS AGE 
lunaize(ggplot(data = regionDF_new_age %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult')))%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = logAmp)) + geom_line(aes(group=interaction(luna,Region), color =Region), alpha = 0.8) + geom_point(aes(color=Region),alpha=.8) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + scale_color_manual(values=c("gold3", "blue4")) + ylab("Evoked") +xlab("Age")

gam.model <-  gam(logAmp ~ s(age, k = 3) + Region, data = regionDF_new_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
summary(gam.model)

regionDF_new_age$oRegion <- ordered(regionDF_new_age$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
model_formula <- as.formula("logAmp ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = regionDF_new_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)

p.adjust((2e-16), method = "bonferroni", n = 3)




ggplot(Amplitudes_channels, aes(x = -Y, y = X, fill = logAmp, z = logAmp, label = Channel)) + geom_topo(chan_markers = "text") + scale_fill_distiller(palette = "RdBu") 


 }
 
Spontaneous <- function() {
 
  
  ## SPONTANEOUS VS AGE
  lunaize(ggplot(data = regionDF_new_age %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Spontaneous, by =luna)) + geom_line(aes(group=interaction(luna,Region), color =Region), alpha = 0.2) + geom_point(aes(color=Region),alpha=.8) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + scale_color_manual(values=c("gold3", "blue4"))+ ylab("Spontaneous") +xlab("Age")
  
  gam.model <-  gam(Spontaneous ~ s(age, k = 3) + Region, data = regionDF_new_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  regionDF_new_age$oRegion <- ordered(regionDF_new_age$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("Spontaneous ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = regionDF_new_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
p.adjust((2e-16), method = "bonferroni", n = 3)

 }
 
 

SNR <- function() {
  ## SNR VS AGE
  lunaize(ggplot(data = regionDF_new_age %>% filter(age <= 26) %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Difference, by =luna))  + geom_line(aes(group=interaction(luna,Region), color =Region), alpha = 0.2) + geom_point(aes(color=Region),alpha=.8) + geom_smooth(aes(group = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + scale_color_manual(values=c("gold3", "blue4"))+ ylab("SNR") +xlab("Age")
  
  gam.model <-  gam(Difference ~ s(age, k = 3) + Region, data = regionDF_new_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  regionDF_new_age$oRegion <- ordered(regionDF_new_age$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("Difference ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = regionDF_new_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
 }
 
 
SNRvsBestSaccade <- function() {
  
lunaize(ggplot(data = Amplitudes_Frontal %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))), 
               aes(y = mgsLatency_sd, x = zscoredDifference, color = Region)) + stat_smooth(method = 'lm', alpha = 0.1))

  
gam.model <- gam(mgsLatency_sd ~ s(zscoredDifference, k=4, by = as.factor(Region)) + s(age, k=4),
                 data = Amplitudes_Frontal %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))))

summary(gam.model)


lm.model <- lm(mgsLatency_sd ~ zscoredDifference*age + Region , data = Amplitudes_Frontal %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))))
anova(lm.model)

}


SNRvsBestSaccadeVar <- function() {

gam.model <- gam(absBestError_sd ~ s(zscoredDifference, k=3, by=ageGroup) + s(age, k=3),
                 data = ASSR40 %>% filter(Region == 'Left DLPFC') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))))

summary(gam.model)


lunaize(ggplot(data = ASSR40 %>% filter(Region == "Left DLPFC") %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))), 
               aes(y = absBestError_sd, x = zscoredDifference, color = ageGroup)) + geom_point()+geom_smooth(method = "gam", formula = y ~ splines::bs(x, 3)) + xlab("Signal to Noise Ratio") + ylab("Best Saccade Var (degs)") + facet_wrap( . ~ Region, ncol = 8)) 


lunaize(ggplot(data = ASSR40 %>% filter(Region == "Left DLPFC") %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))), 
               aes(y = absBestError_sd, x = zscoredDifference, color = ageGroup)) + geom_point() + stat_smooth(method = 'lm') + xlab("Signal to Noise Ratio") + ylab("Best Saccade Var (degs)") + facet_wrap( . ~ Region, ncol = 8)) 



# all subjects 
geom_smooth(method = "gam", formula = y ~ splines::bs(x, 3))

lunaize(ggplot(data = ASSR40 %>% filter(Region == "Left DLPFC") %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))), 
               aes(y = absBestError_sd, x = zscoredDifference)) + geom_point()+ stat_smooth(method='lm')
 + xlab("Signal to Noise Ratio") + ylab("Best Saccade Var (degs)") + facet_wrap( . ~ Region, ncol = 8)) 

lm.model <- lm(absBestError_sd ~ zscoredDifference + age, data = ASSR40 %>% filter(Region == 'Left DLPFC'))
anova(lm.model)

gam.model <- gam(absBestError_sd ~ s(zscoredDifference, k=3) + s(age, k=2),
                 data = ASSR40 %>% filter(Region == 'Left DLPFC') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))))
summary(gam.model)

AIC(lm.model, gam.model)

}

AvgEvokedSpontanous <- function() {

## average the evoked and spontaneous activity together to get a "whole brain" measure 
WholeBrainAmplitudes <- aggregate(.~ Subject, Amplitudes_age_new , mean)
WholeBrainAmplitudes_new <- WholeBrainAmplitudes


cols = names(WholeBrainAmplitudes[7:8])
for (col in cols) {
  
  
  indx <- outliers(WholeBrainAmplitudes[[col]])
  
  WholeBrainAmplitudes_new[[col]] <- Map(replace, WholeBrainAmplitudes_new[[col]], indx, NA)
  WholeBrainAmplitudes_new[[col]] <- as.numeric(WholeBrainAmplitudes_new[[col]])
  
}  
WholeBrainAmplitudes_new$SNR <- WholeBrainAmplitudes_new$Amplitude / WholeBrainAmplitudes_new$Spontaneous

cols = names(WholeBrainAmplitudes[9])
for (col in cols) {
  
  
  indx <- outliers(WholeBrainAmplitudes_new[[col]])
  
  WholeBrainAmplitudes_new[[col]] <- Map(replace, WholeBrainAmplitudes_new[[col]], indx, NA)
  WholeBrainAmplitudes_new[[col]] <- as.numeric(WholeBrainAmplitudes_new[[col]])
  
}  

lunaize(ggplot(data = WholeBrainAmplitudes_new[], aes(x = age, y = Spontaneous)) + geom_point() + stat_smooth(method = "lm", formula = 'y ~ I(1/x)') + ggtitle("") + xlab("") + ylab("") )

lunaize(ggplot(data = WholeBrainAmplitudes_new[], aes(x = age, y = SNR)) + geom_point() + stat_smooth(method = "gam") + ggtitle("") + xlab("") + ylab("") )


lm.model <- lm(data = WholeBrainAmplitudes_new[], (SNR) ~ (inverseAge))
summary(lm.model)
}


 
SNRvsMRS <- function() {
  
  regionDF <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_40hz_SNR.csv')
  regionDF_new_age <- regionDF %>% merge(., agefile, by = "Subject")
  regionDF_new_age <- regionDF_new_age %>% separate(Subject,c("luna","vdate"), remove=F)
  rbind40and20 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/EvokedSpontaneousBehavior_40and20Hz_20230209.csv')
  MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  
  idx <- which(!is.na(MRSregionsRes$GABA.Cr.adj) & !is.na(MRSregionsRes$Glu.Cr.adj))
  gabaglu.lm <- lm(Glu.Cr.adj ~ GABA.Cr.adj + Region, data = MRSregionsRes[idx,])
  MRSregionsRes$GluGABAimbalance <- NA
  MRSregionsRes[idx,]$GluGABAimbalance <- abs(gabaglu.lm$residuals)
  
  ## MERGE SNR AND FOOOF
  SNRMRS <- merge(MRSregionsRes, rbind40and20, by = c("luna", "Region", "visitno")) %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult')))
  
  
  ## gaba GLU imbalance VS SNR    scale_color_gradient2(midpoint=22, low="gold1", mid="blue2", high="blue4")
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Difference, x = GluGABAimbalance, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.8) + scale_color_manual(values=c("gold3", "blue4"))) + ylab("SNR")
  
  lm.model <- lmer(Difference ~ GluGABAimbalance + age.x  + Region + Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Difference ~ GluGABAimbalance*age.x + Region + Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## GLU VS Evoked    scale_color_gradient2(midpoint=22, low="gold1", mid="blue2", high="blue4")
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = logAmp, x = gluRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.8) + scale_color_manual(values=c("gold3", "blue4")))
  
  lm.model <- lmer(logAmp ~ gluRes.residuals + age.x  + Region + Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(logAmp ~ gluRes.residuals*age.x + Region + Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(logAmp ~ gluRes.residuals*Frequency+Frequency*age.x + Region +  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)

  ##if significant age interaction
  lm.model <- lmer(logAmp ~ gluRes.residuals + age.x  + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(ageGroup == "Adult"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## GLU VS spontaneous  
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Spontaneous, x = gluRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))
  
  lm.model <- lmer(Spontaneous ~ gluRes.residuals + age.x + Region + Frequency+ (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Spontaneous ~ gluRes.residuals * age.x +Region+ Frequency +  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(Spontaneous ~ gluRes.residuals*Frequency +Frequency*age.x + Region +  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  ##if significant main effect of region
  lm.model <- lm(logAmp ~ gluRes.residuals + age.x + Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## GLU VS SNR  
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Difference, x = gluRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))
  
  lm.model <- lmer(Difference ~ gluRes.residuals + age.x + Region + Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Difference ~ gluRes.residuals * age.x +Region + Frequency +(1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(Difference ~ gluRes.residuals*Frequency +Frequency*age.x + Region +  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ##if significant age interaction
  lm.model <- lmer(Difference ~ gluRes.residuals + age.x  + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(ageGroup == "Adol"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## GABA VS evoked 
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
  aes(y = logAmp, x = gabaRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + 
  geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))


  lm.model <- lmer(logAmp ~ gabaRes.residuals + age.x  + Frequency + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)

  lm.model <- lmer(logAmp ~ gabaRes.residuals * age.x + Region + Frequency+  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(logAmp ~ gabaRes.residuals *Frequency+Frequency*age.x + Region +   (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  ##if significant age interaction
  lm.model <- lmer(logAmp ~ gabaRes.residuals + age.x  + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(ageGroup == "Adult"))
  summ(lm.model)
  
  
  ## GABA VS spontaneous
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Spontaneous, x = gabaRes.residuals, by = luna, color = ageGroup, shape = Region)) + geom_line(aes(group=interaction(luna, Region, Frequency)), alpha = 0.1) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))
  
  lm.model <- lmer(Spontaneous ~ gabaRes.residuals + age.x + Frequency + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
    lm.model <- lmer(Spontaneous ~ gabaRes.residuals * age.x+ Frequency + Region+  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(Spontaneous ~ gabaRes.residuals *Frequency+Frequency*age.x + Region +   (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  ## GABA VS SNR
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Difference, x = gabaRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(Difference ~ gabaRes.residuals + age.x +Frequency + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  lm.model <- lmer(Difference ~ gabaRes.residuals * age.x + Region +Frequency+  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(Difference ~ gabaRes.residuals *Frequency+Frequency*age.x + Region +   (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
 
  ## Ratio VS Evoked
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = logAmp, x = ratioRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(logAmp ~ ratioRes.residuals + age.x +Frequency  + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  lm.model <- lmer(logAmp ~ ratioRes.residuals * age.x + Region +Frequency +  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(logAmp ~ ratioRes.residuals *Frequency+Frequency*age.x + Region +   (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
   ## Ratio VS Spontaneous
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Spontaneous, x = ratioRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(Spontaneous ~ ratioRes.residuals + age.x  + Region +Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  lm.model <- lmer(Spontaneous ~ ratioRes.residuals * age.x + Region +Frequency+  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(Spontaneous ~ ratioRes.residuals *Frequency+Frequency*age.x + Region +   (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  ##if significant age interaction
  lm.model <- lmer(Spontaneous ~ ratioRes.residuals + age.x  + Region + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(ageGroup == "Adult"))
  summ(lm.model)
  
  ## Ratio VS SNR
  lunaize(ggplot(data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Difference, x = ratioRes.residuals, by = luna, color = ageGroup)) + geom_line(aes(group=interaction(luna, Region)), alpha = 0.2) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Frequency, linetype = Frequency), method="lm", alpha = 0.5) + scale_color_manual(values=c("gold3", "blue4")))
  
  lm.model <- lmer(Difference ~ ratioRes.residuals + age.x  + Region +Frequency + (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  lm.model <- lmer(Difference ~ ratioRes.residuals * age.x + Region +Frequency+  (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  #mrs by freq interaction
  lm.model <- lmer(Difference ~ ratioRes.residuals *Frequency+Frequency*age.x + Region +   (1|luna), data = SNRMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  p.adjust((0.78), method = "bonferroni", n = 3)
  
  
  
    ## Mediation
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
  
  
  mediationMatrix <- fooofMRS %>% dplyr::select(gluRes.residuals, Offset, Condition, Region, luna, visitno, age.x) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on offset (c)
  model.0 <- lmer(Offset ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  summ(model.0)
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on Ratio  (a)
  model.M <- lmer(gluRes.residuals ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  summ(model.M)
  
  #the effect of MRS ratio on offset (b)
  model.Y <- lmer(Offset ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = 'gluRes.residuals', boot = FALSE, sims = 1000)
  (summary(results))
  
  
  
  
}



