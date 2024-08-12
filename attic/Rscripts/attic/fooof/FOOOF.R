
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

createDataFrames <- function() {
  
  wholeBrain <- function() {
    
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/PrepMRSdata.R", envir = knitr::knit_global(), chdir = TRUE)
    agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsFooofMeasures_20230111.csv')
    
    outliers <- function(x) {
      
      (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
      
      
    }
    
    fooof_new <- fooof
    
    #outliers per channel
    cols = names(fooof[4:5])
    for (col in cols) {
      
      fooof_grouped <- fooof %>% group_by(Subject) %>% group_by(Condition)
      
      indx <- outliers(fooof_grouped[[col]])
      
      fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
      fooof_new[[col]] <- as.numeric(fooof_new[[col]])
      
    }  
    
    fooof_age <- merge(fooof_new, agefile, by = 'Subject', all.x = T, all.y = T)
    #fooof_age<- fooof_age %>% subset(visitno ==1)
    
    lunaize(ggplot(data = fooof_age , aes(x = age, y = Offset, by =Channel )) + geom_smooth( method = "lm", formula = y ~ I(1/x))) + facet_wrap(.~Condition)
    
    
    # fooof_age_chanLocs <- merge(fooof_age, chanLocs, by = "Channel")
    
    fooof_age_wholeScalp <- aggregate(.~Subject+Condition+sex, fooof_age[c(1,4:9)], mean)
    fooof_age_wholeScalp <- fooof_age_wholeScalp %>% separate(Subject,c("luna","vdate"), remove=F)
    
    # outliers per subject 
    fooof_age_wholeScalp_new <- fooof_age_wholeScalp
    
    cols = names(fooof_age_wholeScalp[5:6])
    for (col in cols) {
      
      fooof_grouped <- fooof_age_wholeScalp %>% group_by(Subject) %>% group_by(Condition)
      
      indx <- outliers(fooof_grouped[[col]])
      
      fooof_age_wholeScalp_new[[col]] <- Map(replace, fooof_age_wholeScalp_new[[col]], indx, NA)
      fooof_age_wholeScalp_new[[col]] <- as.numeric(fooof_age_wholeScalp_new[[col]])
      
    }
  
    
    fooof_age_wholeScalp_new <- fooof_age_wholeScalp_new%>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult')))
    write.csv(fooof_age_wholeScalp_new, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/WholeBrainFOOOFMeasures_20230111.csv')
  }
  
  allChannels <- function (){
    
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/PrepMRSdata.R", envir = knitr::knit_global(), chdir = TRUE)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsFooofMeasures_20230111.csv')
    
    outliers <- function(x) {
      
      (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
      
      
    }
    
    fooof_new <- fooof
    
    #outliers per channel
    cols = names(fooof[4:5])
    for (col in cols) {
      
      fooof_grouped <- fooof %>% group_by(Subject) %>% group_by(Condition)
      
      indx <- outliers(fooof_grouped[[col]])
      
      fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
      fooof_new[[col]] <- as.numeric(fooof_new[[col]])
      
    }  
    
    fooof_age <- merge(fooof_new, agefile, by = 'Subject')
    #fooof_age<- fooof_age %>% subset(visitno ==1)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof_age <- subset(fooof_age, select = c(-X))
  
    fooof_channels <- merge(fooof_age, chanLocs, by = c("Channel"))
    
  
    write.csv(fooof_channels, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allChannelsFOOOFMeasures_20230214.csv')
    
    
  }
  
  
  regions <- function (){
    
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/PrepMRSdata.R", envir = knitr::knit_global(), chdir = TRUE)

    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsFooofMeasures_20230111.csv')
    
    outliers <- function(x) {
      
      (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
      
      
    }
    
    fooof_new <- fooof
    
    #outliers per channel
    cols = names(fooof[4:5])
    for (col in cols) {
      
      fooof_grouped <- fooof %>% group_by(Subject) %>% group_by(Condition)
      
      indx <- outliers(fooof_grouped[[col]])
      
      fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
      fooof_new[[col]] <- as.numeric(fooof_new[[col]])
      
    }  
    
    fooof_age <- merge(fooof_new, agefile, by = 'Subject')
    #fooof_age<- fooof_age %>% subset(visitno ==1)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels

    fooof_channels <- merge(fooof_age, chanLocs, by = c("Channel"))
    
    regionDF <- fooof_channels %>% filter(., Channel == 'F3' | Channel == 'F5'| Channel == 'F7') %>% group_by(Subject, Condition) %>% 
      summarise(Offset = mean(Offset, na.rm=T), 
                Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "LDLPFC") %>% 
      rbind(.,fooof_channels %>% filter(., Channel == 'F4' | Channel == 'F6'| Channel == 'F8') %>% group_by(Subject, Condition) %>% 
              summarise(Offset = mean(Offset, na.rm=T), 
                        Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "RDLPFC")) %>%
      rbind(.,fooof_channels %>% filter(., str_detect(fooof_channels$Channel, "F")) %>% group_by(Subject, Condition) %>% 
              summarise(Offset = mean(Offset, na.rm=T), 
                        Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "MPFC"))

    write.csv(regionDF, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_20230111.csv')
    
    
  }
  
  regions_30to50hz <- function (){
    
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/PrepMRSdata.R", envir = knitr::knit_global(), chdir = TRUE)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/30-50hz/allSubjectsFooofMeasures_30to50hz_20230123.csv')
    
    outliers <- function(x) {
      
      (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
      
      
    }
    
    fooof_new <- fooof
    
    #outliers per channel
    cols = names(fooof[4:5])
    for (col in cols) {
      
      fooof_grouped <- fooof %>% group_by(Subject) %>% group_by(Condition)
      
      indx <- outliers(fooof_grouped[[col]])
      
      fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
      fooof_new[[col]] <- as.numeric(fooof_new[[col]])
      
    }  
    
    fooof_age <- merge(fooof_new, agefile, by = 'Subject')
    #fooof_age<- fooof_age %>% subset(visitno ==1)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof_channels <- merge(fooof_age, chanLocs, by = c("Channel"))
    
    regionDF <- fooof_channels %>% filter(., Channel == 'F3' | Channel == 'F5'| Channel == 'F7') %>% group_by(Subject, Condition) %>% 
      summarise(Offset = mean(Offset, na.rm=T), 
                Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "LDLPFC") %>% 
      rbind(.,fooof_channels %>% filter(., Channel == 'F4' | Channel == 'F6'| Channel == 'F8') %>% group_by(Subject, Condition) %>% 
              summarise(Offset = mean(Offset, na.rm=T), 
                        Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "RDLPFC")) %>%
      rbind(.,fooof_channels %>% filter(., str_detect(fooof_channels$Channel, "F")) %>% group_by(Subject, Condition) %>% 
              summarise(Offset = mean(Offset, na.rm=T), 
                        Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "MPFC"))
    
    write.csv(regionDF, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_30to50hz_20230123.csv')
    
    
  }
  
  
  gammapeaks_regions <- function (){
    
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/PrepMRSdata.R", envir = knitr::knit_global(), chdir = TRUE)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooofGammaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/Peaks/allSubjectsGammaPeakMeasures_20230123.csv')
    
    outliers <- function(x) {
      
      (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
      
      
    }
    
    fooof_new <- fooofGammaPeaks
    
    #outliers per channel
    cols = names(fooofGammaPeaks[4:6])
    for (col in cols) {
      
      fooof_grouped <- fooofGammaPeaks %>% group_by(Subject) %>% group_by(Condition)
      
      indx <- outliers(fooof_grouped[[col]])
      
      fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
      fooof_new[[col]] <- as.numeric(fooof_new[[col]])
      
    }  
    
    fooof_age <- merge(fooof_new, agefile, by = 'Subject')
    #fooof_age<- fooof_age %>% subset(visitno ==1)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof_channels <- merge(fooof_age, chanLocs, by = c("Channel"))
    

    
    is.na(fooof_channels$Center.Frequency) <- fooof_channels$Center.Frequency==0 #if there isnt a peak, python puts a 0, make NA
    is.na(fooof_channels$Bandwidth) <- fooof_channels$Bandwidth==0 #if there isnt a peak, python puts a 0, make NA
    is.na(fooof_channels$Power) <- fooof_channels$Power==0 #if there isnt a peak, python puts a 0, make NA
    
    
    regionDF <- fooof_channels %>% filter(., Channel == 'F3' | Channel == 'F5'| Channel == 'F7') %>% group_by(Subject, Condition) %>% 
      summarise(Center.Frequency = mean(Center.Frequency, na.rm=T), 
                Bandwidth = mean(Bandwidth, na.rm=T),
                Power = mean(Power, na.rm=T)) %>% mutate(Region = "LDLPFC") %>% 
      rbind(.,fooof_channels %>% filter(., Channel == 'F4' | Channel == 'F6'| Channel == 'F8') %>% group_by(Subject, Condition) %>% 
              summarise(Center.Frequency = mean(Center.Frequency, na.rm=T), 
                        Bandwidth = mean(Bandwidth, na.rm=T),
                        Power = mean(Power, na.rm=T)) %>% mutate(Region = "RDLPFC")) %>%
      rbind(.,fooof_channels %>% filter(., str_detect(fooof_channels$Channel, "F")) %>% group_by(Subject, Condition) %>% 
              summarise(Center.Frequency = mean(Center.Frequency, na.rm=T), 
                        Bandwidth = mean(Bandwidth, na.rm=T),
                        Power = mean(Power, na.rm=T)) %>% mutate(Region = "MPFC"))
    
    write.csv(regionDF, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalGammaPeaks_20230123.csv')
    
    
  }
  
  allBandPeaks_regions <- function (){
    
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
    sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/PrepMRSdata.R", envir = knitr::knit_global(), chdir = TRUE)
    
    chanLocs <- read.csv(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv'))
    chanLocs$Channel <- chanLocs$labels
    
    fooofPeaks <- read.csv(hera('/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsPeakMeasures_20230710.csv'))
    fooofPeaks$urchan <- fooofPeaks$Channel+1
    fooofPeaks <- subset(fooofPeaks, select = c(-Channel))
    
    outliers <- function(x) {
      
      (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
      
      
    }
    
    fooof_new <- fooofPeaks
    
    #outliers per channel
    cols = names(fooofPeaks[3:5])
    for (col in cols) {
      
      fooof_grouped <- fooofPeaks %>% group_by(Subject) %>% group_by(Condition)
      
      indx <- outliers(fooof_grouped[[col]])
      
      fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
      fooof_new[[col]] <- as.numeric(fooof_new[[col]])
      
    }  
    
    fooof_age <- merge(fooof_new, agefile, by = 'Subject')
    #fooof_age<- fooof_age %>% subset(visitno ==1)
    
    chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
    chanLocs$Channel <- chanLocs$labels
    
    fooof_channels <- merge(fooof_age, chanLocs, by = c("urchan"))
    
    
    
    is.na(fooof_channels$Center.Frequency) <- fooof_channels$Center.Frequency==0 #if there isnt a peak, python puts a 0, make NA
    is.na(fooof_channels$Bandwidth) <- fooof_channels$Bandwidth==0 #if there isnt a peak, python puts a 0, make NA
    is.na(fooof_channels$Power) <- fooof_channels$Power==0 #if there isnt a peak, python puts a 0, make NA
    
    fooof_channels_bands <- fooof_channels %>% filter(Center.Frequency >= 13 & Center.Frequency < 30) %>% mutate(Band = "Beta") %>%
      rbind(.,fooof_channels %>% filter(Center.Frequency >=4 & Center.Frequency < 8 ) %>% mutate(Band = "Theta")) %>%
      rbind(.,fooof_channels %>% filter(Center.Frequency >= 8  & Center.Frequency < 13) %>% mutate(Band = "Alpha")) %>%
      rbind(.,fooof_channels %>% filter(Center.Frequency >= 30 & Center.Frequency < 75) %>% mutate(Band = "Gamma"))
      
    
    
    regionDF <- fooof_channels_bands %>% filter(., Channel == 'F3' | Channel == 'F5'| Channel == 'F7') %>% group_by(Subject, Condition, Band) %>% 
      summarise(Center.Frequency = mean(Center.Frequency, na.rm=T), 
                Bandwidth = mean(Bandwidth, na.rm=T),
                Power = mean(Power, na.rm=T)) %>% mutate(Region = "LDLPFC") %>% 
      rbind(.,fooof_channels_bands %>% filter(., Channel == 'F4' | Channel == 'F6'| Channel == 'F8') %>% group_by(Subject, Condition, Band) %>% 
              summarise(Center.Frequency = mean(Center.Frequency, na.rm=T), 
                        Bandwidth = mean(Bandwidth, na.rm=T),
                        Power = mean(Power, na.rm=T)) %>% mutate(Region = "RDLPFC")) %>%
      rbind(.,fooof_channels_bands %>% filter(., str_detect(fooof_channels_bands$Channel, "F")) %>% group_by(Subject, Condition, Band) %>% 
              summarise(Center.Frequency = mean(Center.Frequency, na.rm=T), 
                        Bandwidth = mean(Bandwidth, na.rm=T),
                        Power = mean(Power, na.rm=T)) %>% mutate(Region = "MPFC"))
    
    write.csv(regionDF, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalPeaks_allBands.csv')
    
    
  }
  
  
  MRSdf <- function() {
    
    LDLPFC_GABA <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/LDLPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
    LDLPFC_Glu <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/LDLPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
    LDLPFC_Ratio <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/LDLPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
    
    LDLPFC <- LDLPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
      merge(., LDLPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"),all.x = T, all.y = T) %>%
      merge(., LDLPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
      mutate(Region = "LDLPFC")
    
    
    RDLPFC_GABA <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/RDLPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
    RDLPFC_Glu <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/RDLPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
    RDLPFC_Ratio <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/RDLPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
    
    RDLPFC <- RDLPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
      merge(., RDLPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"),all.x = T, all.y = T) %>%
      merge(., RDLPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
      mutate(Region = "RDLPFC")
    
    
    MPFC_GABA <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/MPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
    MPFC_Glu <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/MPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
    MPFC_Ratio <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/MPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
    
    MPFC <- MPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
      merge(., MPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"), all.x = T, all.y = T) %>%
      merge(., MPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
      mutate(Region = "MPFC")
    
    MRSregions <- rbind(LDLPFC, RDLPFC) %>% rbind(., MPFC) %>% mutate(luna = subjID)%>% mutate(visitno = visitnum)
    write.csv(MRSregions, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj.csv')
    
    gluRes <- lm(Glu.Cr.adj~GMrat, data = MRSregions) 
    Subject <- dplyr::select(MRSregions[!is.na(MRSregions$GMrat & MRSregions$Glu.Cr.adj),],Subject,Region)
    gluResDF <- data.frame(Subject, gluRes$residuals)
    
    gabaRes <- lm(GABA.Cr.adj~GMrat, data = MRSregions) 
    Subject <- dplyr::select(MRSregions[!is.na(MRSregions$GMrat & MRSregions$GABA.Cr.adj),],Subject,Region)
    gabaResDF <- data.frame(Subject, gabaRes$residuals)
    
    ratioRes <- lm(Ratio.adj~GMrat, data = MRSregions) 
    Subject <- dplyr::select(MRSregions[!is.na(MRSregions$GMrat & MRSregions$Ratio.adj),],Subject,Region)
    ratioResDF <- data.frame(Subject, ratioRes$residuals)
    
    MRSregionsRes <- merge(MRSregions, gluResDF, by = c("Subject", "Region"), all = T) %>% merge(.,gabaResDF, by = c("Subject", "Region"), all = T) %>% merge(.,ratioResDF, by = c("Subject", "Region"), all = T)
    write.csv(MRSregionsRes, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  }
  
  
}
fooof_allChannels <- function() {
  
  fooofallChannels <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allChannelsFOOOFMeasures_20230214.csv')
  fooofallChannels <- fooofallChannels %>% separate(Subject,c("luna","vdate"), remove=F)
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs$Channel <- chanLocs$labels
  
  
  
  lunaize(ggplot(fooofallChannels, aes(x = -Y, y = X, fill =Offset, z = Offset, label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.9, low="navy", mid="lightblue2", high="lightblue1")) + ggtitle("Offset") + theme(text = element_text(size = 30))
  
  
  exponentAgeCoef <- data.frame()
  for (chan in unique(fooofallChannels$Channel)) {
    ageRes <- lm(z~invage + Condition, data = fooofallChannels %>% filter(Channel == chan) %>% mutate(z = scale(Exponent)[,1], invage=1/age) )
    exponentAgeCoef <- rbind(exponentAgeCoef, data.frame(chan,ageRes$coefficients[2]))
  } 
  
  exponentAgeCoef <- merge(exponentAgeCoef, chanLocs, by = "labels")
  exponentAgeCoef$Measure <- "Exponent"
  
  
  offsetAgeCoef <- data.frame()
  for (chan in unique(fooofallChannels$Channel)) {
    ageRes <- lm(z~invage + Condition, data = fooofallChannels %>% filter(Channel == chan) %>% mutate(z = scale(Offset)[,1], invage=1/age))   
                 offsetAgeCoef <- rbind(offsetAgeCoef, data.frame(chan,ageRes$coefficients[2]))
  }
  
  offsetAgeCoef$labels <- offsetAgeCoef$chan
  offsetAgeCoef <- merge(offsetAgeCoef, chanLocs, by = "labels")
  offsetAgeCoef$Measure <- "Offset"
  
  allMeasues <- rbind(exponentAgeCoef,offsetAgeCoef)
  
  lunaize(ggplot(allMeasues, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=-20, low="navy", mid="lightblue2", high="lightblue1")) + ggtitle("Age Effects") +facet_wrap(.~Measure)+ theme(text = element_text(size = 30))
  
  
  ggplot(offsetAgeCoef, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_distiller(palette = "RdBu")+ ggtitle("Offset Age Effects")
  
  lunaize(ggplot(fooofallChannels %>% mutate(isDLPFC = labels %in% c('F3','F5','F7','F4','F6','F8')), aes(x = age, y = Exponent, color = isDLPFC)) + geom_smooth(aes(group = labels), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.1,size=1)) + scale_color_manual(values=c("lightgray", "blue4"))
 
  
   lunaize(ggplot(fooofallChannels %>% mutate(isDLPFC = labels %in% c('F3','F5','F7','F4','F6','F8')), aes(x = age, y = Offset, color = isDLPFC)) + geom_smooth(aes(group = labels), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.1,size=1)) + scale_color_manual(values=c("lightgray", "blue4"))
}

fooof_wholeBrain <- function() {
  
  fooof_wholeBrain <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/WholeBrainFOOOFMeasures_20230111.csv')
  

#check age vs offset and age vs exponent 
lunaize(ggplot(data = fooof_wholeBrain , aes(x = age, y = Exponent, by =luna)) + geom_point() +geom_line(aes(group=luna))+ geom_smooth(aes(group = 1), method = "lm", formula = y ~ I(1/x))) + facet_wrap(.~Condition)

gam.model <-  gam(Offset ~ s(age, k = 3)  + Condition + s(age, k=3, by=Condition), data = fooof_wholeBrain, random = ~(1|luna))

summary(gam.model)

#check offset vs exponent they should be very correlated
lunaize(ggplot(data = fooof_wholeBrain , aes(x = Offset, y = Exponent, by =luna , color=age)) + geom_point() +geom_line(aes(group=luna))+ geom_smooth(aes(group = 1), method = "lm")) + facet_wrap(.~Condition) + scale_color_gradient2(midpoint=20, low="gold3", mid="white", high="blue4")

lm.model <- lmer(Exponent ~ Offset + Offset*age + (1|luna), data =fooof_wholeBrain )
car::Anova(lm.model)

#trying things with different types of plots 
ggplot(fooof_wholeBrain, aes(x=ageGroup, y=Exponent, color=ageGroup)) + geom_beeswarm(alpha=0.7, size=4, cex=1.5) + theme_light() + facet_wrap(.~Condition) + scale_color_manual(values=c("gold3", "blue4")) + stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = 'mean', geom = 'crossbar', color = 'firebrick', width=0.5)

}

waterfall_group <- function(d) {
  if (!all(c("luna", "age") %in% names(d))){
    stop("dataframe must have columns 'luna' and 'age'")
  }
  age_ranked <-
    d %>%
    select(luna, age) %>%
    group_by(luna) %>%
    summarise(minage=min(age)) %>%
    ungroup() %>%
    mutate(age_id = rank(minage, ties.method="first") ) %>%
    inner_join(d, by="luna")
}

waterfall_plot <- function(d, ...) {
  
  age_ranked <- waterfall_group(d)
  
  p <-
    ggplot(age_ranked) +
    aes(x=age, y=age_id, group=age_id, color = sex) +
    geom_line() +
    geom_point()
  
  p <- lunaize(p) +
    ylab("") +
    xlab("Age") +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y  = element_blank()
    )
  
  return(p)
}

fooof_regions <- function() {
 
  fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_20230111.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  
  fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F)
  
 
  
  ##Exponent vs Offset
  lunaize(ggplot(data = fooof_regional_age %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent, x = Offset, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.4) + 
            geom_point(alpha=.4) + geom_smooth(aes(group = Condition, linetype = Condition), method="lm", alpha = 0.9, size = 1) + 
            scale_color_manual(values=c("gold3", "blue4")))
 
  corDF <- fooof_regional_age %>% filter(Condition == "eyesOpen") %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  cor.test(corDF$Exponent , corDF$Offset, method = "pearson")
  
  lm.model <- lmer(Exponent ~ Offset + age + Condition + Region + (1|luna), data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model) 
  
  lm.model <- lmer(Exponent ~ Offset*age + Condition + Region + (1|luna), data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model) 
  
 ##Exponent vs age
  lunaize(ggplot(data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Exponent)) + geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + geom_point(aes(shape=Region),alpha=.5) + geom_smooth(aes(group = sex, color = sex, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + scale_color_manual(values=c("pink", "blue")) + facet_wrap(~Condition) + theme(text = element_text(size = 30)) + xlab("Age") + ylab("Exponent")
  
  gam.model <-  gam(Exponent ~ s(age, k = 3)  + Condition + Region + sex, data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  # condition interaction, controlling for region
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
               random = list(luna=~1),
               data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  # region interaction, controlling for condition
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Exponent ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  
  ##offset vs age
  lunaize(ggplot(data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Offset)) + geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + geom_point(aes(shape=Region),alpha=.5) + geom_smooth(aes(group = sex, color = sex, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + scale_color_manual(values=c("pink", "blue")) + facet_wrap(~Condition)+ theme(text = element_text(size = 30)) + xlab("Age") + ylab("Offset")
  
  gam.model <-  gam(Offset ~ s(age, k = 3)  + Condition + Region + sex, data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  
  # region interaction, controlling for condition
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Offset ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  #condition interaction, controlling for region
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Offset ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
}
fooof_30to50hz_regions <- function() {
  
  fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_30to50hz_20230123.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  
  fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F)
  
  
  ##Exponent vs age
  lunaize(ggplot(data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Exponent, by =luna)) + geom_line(alpha=.3,aes(group=luna)) + geom_point(alpha=.3) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=1) + facet_wrap(.~Region+Condition))
  
  
  gam.model <-  gam(Exponent ~ s(age, k = 3)  + Condition + Region, data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  
  ##offset vs age
  lunaize(ggplot(data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Offset, by =luna)) + geom_line(alpha=.3,aes(group=luna)) + geom_point(alpha=.3) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=1) + facet_wrap(.~Region+Condition))
  
  gam.model <-  gam(Offset ~ s(age, k = 3)  + Condition + Region, data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Offset ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
}


peaks <- function(){
  # Gamma Peaks 
  fooofPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalPeaks_allBands.csv')
  fooofPeaks_age <- merge(fooofPeaks, agefile, by = 'Subject')
  fooofPeaks_age <- fooofPeaks_age %>% separate(Subject,c("luna","vdate"), remove=F) 
  fooofPeaks_age$luna <- as.factor(fooofPeaks_age$luna)
  
  
  ##power vs age
  lunaize(ggplot(data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC")%>% filter(Band == "Beta"), aes(x = age, y = Power, by =luna, color=Band)) + geom_line(alpha=.1,aes(group=interaction(luna, Region, Band))) + geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.5,size=1) + facet_wrap(.~Region+Condition))
  
  
  gam.model <-  gam(Power ~ s(age, k = 3)  + Condition + Region, data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"), random = ~(1|luna))
  summary(gam.model)
  
  fooofPeaks_age$oCondition <- ordered(fooofPeaks_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Power ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region + Band")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  
  ##Bandwidth vs age
  lunaize(ggplot(data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Bandwidth, by =luna, color=Band)) + geom_line(alpha=.1,aes(group=interaction(luna, Region, Band))) + geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.5,size=1) + facet_wrap(.~Region+Condition))
  
  
  gam.model <-  gam(Bandwidth ~ s(age, k = 3)  + Condition + Region, data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC")%>% filter(Band == "Beta"), random = ~(1|luna))
  summary(gam.model)
  
  fooofPeaks_age$oCondition <- ordered(fooofPeaks_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Bandwidth ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region + Band")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  ##center frequency vs age
  lunaize(ggplot(data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Center.Frequency, by =luna, color=Band)) + geom_line(alpha=.1,aes(group=interaction(luna, Region, Band))) + geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.5,size=1) + facet_wrap(.~Region+Condition))
  
  
  gam.model <-  gam(Center.Frequency ~ s(age, k = 3)  + Condition + Region, data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"), random = ~(1|luna))
  summary(gam.model)
  
  fooofPeaks_age$oCondition <- ordered(fooofPeaks_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Center.Frequency ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region + Band")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooofPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
}

gammaPeaks <- function() {

# Gamma Peaks 
fooofGammaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalGammaPeaks_20230123.csv')
fooofGammaPeaks_age <- merge(fooofGammaPeaks, agefile, by = 'Subject')
fooofGammaPeaks_age <- fooofGammaPeaks_age %>% separate(Subject,c("luna","vdate"), remove=F) 
fooofGammaPeaks_age$luna <- as.factor(fooofGammaPeaks_age$luna)


##Center Frequency vs age
lunaize(ggplot(data = fooofGammaPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = Power, by =luna)) + geom_line(alpha=.3,aes(group=luna)) + geom_point(alpha=.3) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=1) + facet_wrap(.~Region+Condition))


gam.model <-  gam(Exponent ~ s(age, k = 3)  + Condition + Region, data = fooofGammaPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
summary(gam.model)

fooofGammaPeaks_age$oCondition <- ordered(fooofGammaPeaks_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofGammaPeaks_age %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)

}

topoGammaPeaks <- function() {
  
  fooofGammaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsGammaPeakMeasures.csv')
  
  fooofGammaPeaks_age <- merge(fooofGammaPeaks, agefile, by = 'Subject') %>% subset(., select = -c(2) )

  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs$Channel <- chanLocs$labels
  
  gammaPeaks_chan <- merge(fooofGammaPeaks_age,chanLocs, by = "Channel")
  
  
  ggplot(gammaPeaks_chan, aes(x = -Y, y = X, fill = Center.Frequency, z =  Center.Frequency, label = Channel)) + geom_topo() +  scale_fill_gradient(low = "White", high = "purple") + facet_wrap(.~Condition+Adult)
  
  scale_fill_distiller(palette = "PuOr")
  
  
 
}

FOOOFvsSNR_wholeBrain <- function() {
  
  WholeBrainAmplitudes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/wholeBrainSNR.csv')
  
  fooofSNR <- merge(fooof_age_wholeScalp_new, WholeBrainAmplitudes, by = c("Subject", "visitno", "sex", "age"), all.x = T)
  
  lunaize(ggplot(data = fooofSNR , aes(x = Offset, y = Difference, by =luna, color =age)) + geom_point() + geom_line(aes(group=luna))+ geom_smooth(aes(group = 1), method = "lm")) + facet_wrap(.~Condition)+ scale_color_gradient2(midpoint=20, low="gold3", mid="white", high="blue4")
  
  gam.model <-  gam(Difference ~ s(Offset, k=3)  + s(age, k = 3) + s(age, k=3, by=Offset) + Condition + s(age, k=3, by=Condition), data = fooofSNR  %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) , random = ~(1|luna))
  
  summary(gam.model)
  
  
  p.adjust(0.047, method = 'bonferroni', n = 2)
  
  
  ## mediation
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
  
  other_vars <- c("Exponent","Spontaneous")

  mediationAnalysis_longitud(fooofSNR, other_vars)
  

}

FOOOFvsSNR_regions <- function() {
  
  SNRregions <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_40hz_SNR.csv')

  fooofSNR <- merge(fooof_regional_age, SNRregions, by = c("Subject", "Region"), all.x = T)

  # Exponent vs evoked
  lunaize(ggplot(data = fooofSNR %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = logAmp, x = Exponent, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna, Condition)), alpha = 0.5) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(logAmp ~ Exponent + age + Condition + Region + (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(logAmp ~ Exponent *age +Region + Condition  +  (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  corDF <- fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  cor.test(corDF$Exponent , corDF$logAmp, method = "pearson")
  
  # Exponent vs spontaneous
  lunaize(ggplot(data = fooofSNR %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Spontaneous, x = Exponent, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna, Condition)), alpha = 0.5) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(Spontaneous ~ Exponent + age + Condition + Region + (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Spontaneous ~ Exponent *age +Region + Condition  +  (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  corDF <- fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  cor.test(corDF$Exponent , corDF$Spontaneous, method = "pearson")
  
  # Exponent vs snr
  lunaize(ggplot(data = fooofSNR %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Difference, x = Exponent, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna, Condition)), alpha = 0.5) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  
  lm.model <- lmer(Difference ~ Exponent + age + Condition + Region + (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Difference ~ Exponent *age +Region + Condition  +  (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  corDF <- fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  cor.test(corDF$Exponent , corDF$Difference, method = "pearson")
  
  # Offset vs Evoked
  lunaize(ggplot(data = fooofSNR %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = logAmp, x = Offset, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna, Condition)), alpha = 0.5) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(logAmp ~ Offset + age + Condition + Region + (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(logAmp ~ Offset *age +Region + Condition  +  (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  corDF <- fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  cor.test(corDF$Offset , corDF$logAmp, method = "pearson")
  
  # Offset vs spontaneous
  lunaize(ggplot(data = fooofSNR %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Spontaneous, x = Offset, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna, Condition)), alpha = 0.5) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  lm.model <- lmer(Spontaneous ~ Offset + age + Condition + Region + (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Spontaneous ~ Offset *age +Region + Condition  +  (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  corDF <- fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  cor.test(corDF$Offset , corDF$Spontaneous, method = "pearson")
  
  
  # Offset vs SNR
  lunaize(ggplot(data = fooofSNR %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Difference, x = Offset, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna, Condition)), alpha = 0.5) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(Difference ~ Offset + age + Condition + Region + (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Difference ~ Offset *age +Region + Condition  +  (1|luna), data = fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  corDF <- fooofSNR %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  cor.test(corDF$Offset , corDF$Difference, method = "pearson")
  
  
  p.adjust((0.05), method = "bonferroni", n = 4)
  
  
  
  
  ## mediation
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
  
  other_vars <- c("Exponent","Spontaneous")
  
  mediationAnalysis_longitud(fooofSNR, other_vars)
  
  
}

FOOOFvsMRS <- function() {
  fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_20230111.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  
  fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F)
  

  MRSregionsResOld <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  
  idx <- which(!is.na(MRSregionsResOld$GABA.Cr.adj) & !is.na(MRSregionsResOld$Glu.Cr.adj))
  gabagluOld.lm <- lm(Glu.Cr.adj ~ GABA.Cr.adj + Region, data = MRSregionsResOld[idx,])
  MRSregionsResOld$ratioRes.residuals <- NA
  MRSregionsResOld$ratioRes.residuals_abs <- NA
  
  MRSregionsResOld[idx,]$ratioRes.residuals <- (gabaglu.lm$residuals)
  MRSregionsResOld[idx,]$ratioRes.residuals_abs <- abs(gabaglu.lm$residuals)
  
  MRSregionsResOld$GluMinusGABA <- MRSregionsResOld$gluRes.residuals - MRSregionsResOld$gabaRes.residuals
  
  
  ## MERGE MRS AND FOOOF
  fooofMRS <- merge(MRSregionsRes, fooof_regional_age, by = c("luna", "Region", "visitno", "sex")) %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult')))
  

  
  ## GLU VS AGE 
   lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult')))
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = gluRes.residuals, by = luna, shape = Region)) + 
             geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
     scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glutamate")+ theme(text = element_text(size = 30))
   
   gam.model <-  gam(gluRes.residuals ~ s(age, k = 3) + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
   summary(gam.model)
   
   MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
   model_formula <- as.formula("gluRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
   # Note we keep fx = T for reliable p-values.
   model <- gamm(model_formula,
                 random = list(luna=~1),
                 data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summary(model$gam)
   
   
   gam.model <- gam(gluRes.residuals ~ s(age), data = MRSregionsRes)
   gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
   lunaize(gam_growthrate_plot(MRSregionsRes, gam.model, gam.growthrate, agevar = 'age', yvar = 'gluRes.residuals', draw_points = F))
   
   
   ## GABA VS AGE
   lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = gabaRes.residuals, by =luna, shape = Region)) + 
             geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
     scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("GABA") + theme(text = element_text(size = 30))
   
   gam.model <-  gam(gabaRes.residuals ~ s(age, k = 3) + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
   summary(gam.model)
   
   MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
   model_formula <- as.formula("GABA.Cr.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
   # Note we keep fx = T for reliable p-values.
   model <- gamm(model_formula,
                 random = list(luna=~1),
                 data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summary(model$gam)
   
   ## RATIO VS AGE
   lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = ratioRes.residuals, by =luna, shape = Region)) + 
             geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = sex, color = sex, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1)) + 
     scale_color_manual(values=c("gold3", "blue4"))+ xlab("Age") +ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))
   
   
   gam.model <-  gam(ratioRes.residuals ~ s(age, k = 3) + Region + sex, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
   summary(gam.model)
   
   MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
   model_formula <- as.formula("Ratio.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)+GMrat")
   # Note we keep fx = T for reliable p-values.
   model <- gamm(model_formula,
                 random = list(luna=~1),
                 data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summary(model$gam)
   
   
   gam.model <- gam(ratioRes.residuals ~ s(age) + Region, data = MRSregionsRes%>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
   lunaize(gam_growthrate_plot(MRSregionsRes%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), gam.model, gam.growthrate, agevar = 'age', yvar = 'ratioRes.residuals', draw_points = T))
   
   
   ## GLU GABA Imbalance VS AGE 
   lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult')))
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = ratioRes.residuals_abs, by = luna, shape = Region)) + 
             geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = sex, color = sex, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
             scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))
   
   gam.model <-  gam(ratioRes.residuals ~ s(age, k = 3) + GMrat + Region + sex, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
   summary(gam.model)
   
   MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
   model_formula <- as.formula("ratioRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
   # Note we keep fx = T for reliable p-values.
   model <- gamm(model_formula,
                 random = list(luna=~1),
                 data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summary(model$gam)
   
   gam.model <- gam(ratioRes.residuals ~ s(age, k = 3) + Region, data = MRSregionsRes%>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
   lunaize(gam_growthrate_plot(MRSregionsRes%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), gam.model, gam.growthrate, agevar = 'age', yvar = 'ratioRes.residuals', draw_points = T))
   
   ## GLU VS Exponent
   lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                  aes(x = Exponent, y = gluRes.residuals, by = luna, color = ageGroup))+ 
             geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
             geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
             scale_color_manual(values=c("gold3", "blue4"))) + ylab("Glutamate") + xlab("Exponent")+ theme(text = element_text(size = 30))
   
   
   lm.model <- lmer(gluRes.residuals ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   car::Anova(lm.model)
   summ(lm.model)
   
   
   lm.model <- lmer(gluRes.residuals ~ Exponent * age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   car::Anova(lm.model)
   summary(lm.model)
   summ(lm.model)
   
   
   ## GLU VS Offset + scale_color_gradient2(midpoint=20, low="gold2", mid="lightblue4", high="blue4") )
   lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                  aes(x = Offset, y = gluRes.residuals, by = luna, color = ageGroup))+ 
             geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
             geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
             scale_color_manual(values=c("gold3", "blue4"))) + ylab("Glutamte") + xlab("Offset")+ theme(text = element_text(size = 30))
   
   
   lm.model <- lmer( gluRes.residuals ~ Offset + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   car::Anova(lm.model)
   summary(lm.model)
   
   
   lm.model <- lmer(gluRes.residuals ~  Offset * age.x +Condition +Region +   (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   car::Anova(lm.model)
   summary(lm.model)
   summ(lm.model)
   

   ## GABA VS Exponent
   lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                  aes(x = Exponent, y = gabaRes.residuals, by = luna, color = ageGroup))+ 
             geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
             geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
             scale_color_manual(values=c("gold3", "blue4"))) + ylab("GABA") + xlab("Exponent")+ theme(text = element_text(size = 30))
   
   lm.model <- lmer(gabaRes.residuals ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summ(lm.model)
   
   
   lm.model <- lmer(gabaRes.residuals ~ Exponent *  age.x +Condition + Region+ (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summ(lm.model)

   
   ## GABA VS Offset
   lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                  aes(x = Offset, y = gabaRes.residuals, by = luna, color = ageGroup))+ 
             geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
             geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
             scale_color_manual(values=c("gold3", "blue4"))) + ylab("GABA") + xlab("Offset")+ theme(text = element_text(size = 30))
   
   lm.model <- lmer(gabaRes.residuals  ~ Offset + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summ(lm.model)
   
   
   lm.model <- lmer(gabaRes.residuals  ~ Offset * age.x + Condition +Region +   (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   car::Anova(lm.model)
   summ(lm.model)
   
   ## Ratio VS Exponent+ scale_color_manual(values=c("gold3", "blue4")))
             
   lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = Exponent, y = ratioRes.residuals, by = luna))+ 
             geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
             geom_point(alpha=.5) + geom_smooth(aes(group = sex, color = sex, alpha = 0.1), method="lm", alpha = 0.8)  + ylab("Glu/GABA Ratio") +xlab("Exponent")+ theme(text = element_text(size = 30)))
   
   
   lm.model <- lmer(ratioRes.residuals  ~ Exponent + age.x + Condition + Region +sex+ (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summ(lm.model)
   car::Anova(lm.model)
   
   
   lm.model <- lmer(ratioRes.residuals  ~ Exponent * sex + age.x +Condition +Region+   (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summ(lm.model)
   car::Anova(lm.model)
   
   
   ## Ratio VS Offset 
   lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")%>% filter(ageGroup == "Adult"), 
                  aes(x = Offset, y = ratioRes.residuals, by = luna))+ 
             geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
             geom_point(alpha=.5) + geom_smooth(aes(group = sex, color = sex, alpha = 0.1), method="lm", alpha = 0.8) + 
             scale_color_manual(values=c("pink", "blue"))) + ylab("Glu/GABA Ratio") + xlab("Offset")+ theme(text = element_text(size = 30))



   lm.model <- lmer(ratioRes.residuals  ~  Offset + age.x + Condition + sex + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")%>% filter(ageGroup == "Adult"))
   summ(lm.model)
   car::Anova(lm.model)
   
   
   lm.model <- lmer(ratioRes.residuals  ~ Offset *sex +age.x +Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
   summ(lm.model)
   car::Anova(lm.model)
   
   p.adjust((0.15), method = "bonferroni", n = 4)
   
   
   ## gaba glu imbalance VS offset
   lunaize(ggplot(data = fooofMRS 
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = Offset, y = ratioRes.residuals_abs, by =luna)) + 
             geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = sex, color = sex, alpha = 0.01), method="lm",alpha=.8,size=1)) + scale_color_manual(values=c("gold3", "blue4")) + 
     ylab("Glu GABA Imbalance") + xlab("Offset")+ theme(text = element_text(size = 30))
   
   summary(lmerTest::lmer(ratioRes.residuals_abs ~ Offset + age.x + GMrat + Condition + Region + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   summary(lmerTest::lmer(ratioRes.residuals ~ Offset*sex + age.x + GMrat + Region + Condition  +(1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   
   ## gaba glu imbalance VS exponent 
   lunaize(ggplot(data = fooofMRS 
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(y = ratioRes.residuals, x = Exponent, by =luna)) + 
             geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = sex, color = sex, alpha = 0.01), method="lm",alpha=.8,size=1)) + scale_color_manual(values=c("gold3", "blue4"))+ 
     ylab("Glu GABA Imbalance") + xlab("Exponent")+ theme(text = element_text(size = 30))
   
   summary(lmerTest::lmer(ratioRes.residuals ~ Exponent + age.x + GMrat +Condition+ Region + sex+ (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   summary(lmerTest::lmer(ratioRes.residuals ~ Exponent*sex+age.x + GMrat + Region + Condition + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   p.adjust((0.0046), method = "bonferroni", n = 2)
   
   
   ## gaba glu difference VS exponent 
   lunaize(ggplot(data = fooofMRS 
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(y = abs(GluMinusGABA), x = Exponent, by =luna, color = ageGroup)) + 
             geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = 1), method="lm",alpha=.8,size=1)) + scale_color_manual(values=c("gold3", "blue4"))+ 
     ylab("ABS Glu - GABA") + xlab("Exponent")+ theme(text = element_text(size = 30))
   
   summary(lmerTest::lmer(abs(GluMinusGABA) ~ Exponent + age.x +Condition+ Region + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   summary(lmerTest::lmer(abs(GluMinusGABA) ~ Exponent*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   
   
   ## gaba glu difference VS offset 
   lunaize(ggplot(data = fooofMRS 
                  %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(y = abs(GluMinusGABA), x = Offset, by =luna, color = ageGroup)) + 
             geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
             geom_smooth(aes(group = 1), method="lm",alpha=.8,size=1)) + scale_color_manual(values=c("gold3", "blue4"))+ 
     ylab("ABS Glu - GABA") + xlab("Offset")+ theme(text = element_text(size = 30))
   
   summary(lmerTest::lmer(abs(GluMinusGABA) ~ Offset + age.x +Condition+ Region + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   summary(lmerTest::lmer(GluMinusGABA ~ Offset*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
   
   
   
   #effect size of imbalance-exponent on brain electrode map
   
   fooofallChannels <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allChannelsFOOOFMeasures_20230214.csv')
   fooofallChannels <- fooofallChannels %>% separate(Subject,c("luna","vdate"), remove=F)
   
   fooofallChannels_MRS <- merge(fooofallChannels, MRSregionsRes, by =c("luna", "visitno")) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
   
   exponentImbalanceCoef <- data.frame()
   for (chan in unique(fooofallChannels_MRS$Channel)) {
     ageRes <- lm(z~age.x + e + age.x + GMrat + Condition+ Region , data = fooofallChannels_MRS %>% filter(Channel == chan) %>% mutate(z = scale(ratioRes.residuals)[,1],e = scale(Exponent)[,1] ,invage=1/age.x))   
     exponentImbalanceCoef <- rbind(exponentImbalanceCoef, data.frame(chan,ageRes$coefficients[2]))
   }
   
   exponentImbalanceCoef$labels <- exponentImbalanceCoef$chan
   exponentImbalanceCoef <- merge(exponentImbalanceCoef, chanLocs, by = "labels")
  
   lunaize(ggplot(exponentImbalanceCoef, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.03, low="white", mid="red", high="firebrick4") + ggtitle("Effect size of imbalance~exponent controlling for age, condition, region, and GMrat")  + theme(text = element_text(size = 30)))
   
   
   offsetImbalanceCoef <- data.frame()
   for (chan in unique(fooofallChannels_MRS$Channel)) {
     ageRes <- lm(z~age.x + e + age.x + GMrat + Condition+ Region , data = fooofallChannels_MRS %>% filter(Channel == chan) %>% mutate(z = scale(ratioRes.residuals)[,1],e = scale(Offset)[,1] ,invage=1/age.x))   
     offsetImbalanceCoef <- rbind(offsetImbalanceCoef, data.frame(chan,ageRes$coefficients[2]))
   }
   
   offsetImbalanceCoef$labels <- offsetImbalanceCoef$chan
   offsetImbalanceCoef <- merge(offsetImbalanceCoef, chanLocs, by = "labels")
   
   lunaize(ggplot(offsetImbalanceCoef, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.025, low="white", mid="red", high="firebrick4") + ggtitle("Effect size of imbalance~offset controlling for age, condition, region, and GMrat") + theme(text = element_text(size = 30))) 
   
   
   
   ## Johnson Neyman Plots & Condition == "eyesClosed"
   
   fit <- lmer(Offset ~ age*Glu.Cr.adj + Region + Condition + (1|luna),
             data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(visitno == 1))
   
   johnson_neyman(model = fit, pred = Glu.Cr.adj, modx = age, alpha = .05)
   
   
   ## Mediation
   sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
    
   # Glu on imblance 
   mediationMatrix <- fooofMRS %>% dplyr::select(ratioRes.residuals, gluRes.residuals, Condition, Region, luna, visitno, age.x, GMrat) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
   mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
   ratioRes.residualsRes <- lm(ratioRes.residuals~GMrat, data = mediationMatrix) 
   residuals <- ratioRes.residualsRes$residuals
   mediationMatrix$ratioRes.residualsRes <- residuals
   
   
   #the effect of age on offset (c)
   model.0 <- lme4::lmer(ratioRes.residuals ~ age.x + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.0))
   print(summary(model.0))
   summ(model.0)
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   #the effect of age on Ratio  (a)
   model.M <- lme4::lmer(gluRes.residuals ~ age.x + Region + (1|luna), data = mediationMatrix )
   print(car::Anova(model.M))
   print(summary(model.M))
   summ(model.M)
   
   #the effect of MRS ratio on offset (b)
   model.Y <- lme4::lmer(ratioRes.residualsRes ~ gluRes.residuals + age.x + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   summ(model.Y)
   
   
   results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "gluRes.residuals", boot = FALSE, sims = 1000)
   (summary(results))
   
   
   
    # Glu gaba imbalance on fooof
   mediationMatrix <- fooofMRS %>% dplyr::select(ratioRes.residuals, Offset, Condition, Region, luna, visitno, age.x, GMrat) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
   mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
   ratioRes.residualsRes <- lm(ratioRes.residuals~GMrat, data = mediationMatrix) 
   residuals <- ratioRes.residualsRes$residuals
   mediationMatrix$ratioRes.residualsRes <- residuals
   
   
   #the effect of age on offset (c)
   model.0 <- lme4::lmer(Offset ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.0))
   print(summary(model.0))
   summ(model.0)
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   #the effect of age on Ratio  (a)
   model.M <- lme4::lmer(ratioRes.residualsRes ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
   print(car::Anova(model.M))
   print(summary(model.M))
   summ(model.M)
   
   #the effect of MRS ratio on offset (b)
   model.Y <- lme4::lmer(Offset ~ ratioRes.residualsRes + age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   summ(model.Y)
   
   
   results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "ratioRes.residualsRes", boot = FALSE, sims = 1000)
   (summary(results))

   #glu res on fooof 
   mediationMatrix <- fooofMRS %>% dplyr::select(gluRes.residuals, Offset, Condition, Region, luna, visitno, age.x, GMrat) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
   mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
   
   #the effect of age on offset (c)
   model.0 <- lme4::lmer(Offset ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.0))
   print(summary(model.0))
   summ(model.0)
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   #the effect of age on Ratio  (a)
   model.M <- lme4::lmer(gluRes.residuals ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
   print(car::Anova(model.M))
   print(summary(model.M))
   summ(model.M)
   
   #the effect of MRS ratio on offset (b)
   model.Y <- lme4::lmer(Offset ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   summ(model.Y)
   
   results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "gluRes.residuals", boot = FALSE, sims = 1000)
   (summary(results))
  
   
   #gaba res on exponent 
   mediationMatrix <- fooofMRS %>% dplyr::select(gabaRes.residuals, Offset, Condition, Region, luna, visitno, age.x, GMrat) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
   mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
   
   #the effect of age on offset (c)
   model.0 <- lme4::lmer(Offset ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.0))
   print(summary(model.0))
   summ(model.0)
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   #the effect of age on Ratio  (a)
   model.M <- lme4::lmer(gabaRes.residuals ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
   print(car::Anova(model.M))
   print(summary(model.M))
   summ(model.M)
   
   #the effect of MRS ratio on offset (b)
   model.Y <- lme4::lmer(Offset ~ gabaRes.residuals + age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   summ(model.Y)
   
   results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "gabaRes.residuals", boot = FALSE, sims = 1000)
   (summary(results))
   
   # Glu gaba ratio on fooof
   mediationMatrix <- fooofMRS %>% dplyr::select(ratioRes.residuals, Exponent, Condition, Region, luna, visitno, age.x, GMrat) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
   mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
   
   #the effect of age on offset (c)
   model.0 <- lme4::lmer(Exponent ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.0))
   print(summary(model.0))
   summ(model.0)
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   #the effect of age on Ratio  (a)
   model.M <- lme4::lmer(ratioRes.residuals ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
   print(car::Anova(model.M))
   print(summary(model.M))
   summ(model.M)
   
   #the effect of MRS ratio on offset (b)
   model.Y <- lme4::lmer(Exponent ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = mediationMatrix)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   summ(model.Y)
   
   
   results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "ratioRes.residuals", boot = FALSE, sims = 1000)
   (summary(results))
   
}


FOOOF_30to50_vsMRS <- function() {
  fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_30to50hz_20230123.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  
  fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F)
 
  MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  idx <- which(!is.na(MRSregionsRes$GABA.Cr.adj) & !is.na(MRSregionsRes$Glu.Cr.adj))
  gabaglu.lm <- lm(Glu.Cr.adj ~ GABA.Cr.adj + Region, data = MRSregionsRes[idx,])
  MRSregionsRes$ratioRes.residuals <- NA
  MRSregionsRes[idx,]$ratioRes.residuals <- abs(gabaglu.lm$residuals)
  
  ## MERGE MRS AND FOOOF
  fooofMRS <- merge(MRSregionsRes, fooof_regional_age, by = c("luna", "Region", "visitno")) 

  
  ## GLU VS AGE
  lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = gluRes.residuals, by =luna)) + geom_line(aes(group=luna)) + geom_point(alpha=.3) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=1)) + facet_wrap(.~Region)
  
  gam.model <-  gam(gluRes.residuals ~ s(age, k = 3) + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("gluRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  ## GABA VS AGE
  lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = gabaRes.residuals, by =luna)) + geom_line(aes(group=luna)) + geom_point(alpha=.3) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=1)) + facet_wrap(.~Region)
  
  gam.model <-  gam(gabaRes.residuals ~ s(age, k = 3) + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("gabaRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  ## RATIO VS AGE
  lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))%>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = ratioRes.residuals, by =luna)) + geom_line(aes(group=luna)) + geom_point(alpha=.3) + geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=1)) + facet_wrap(.~Region)
  
  gam.model <-  gam(ratioRes.residuals ~ s(age, k = 3) + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("ratioRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  
  ## GLU VS Exponent
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(x = Exponent, y = gluRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + ylab("Glutamate") + xlab("Exponent")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lmer(gluRes.residuals ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(gluRes.residuals ~ Exponent * age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## GLU VS Offset + scale_color_gradient2(midpoint=20, low="gold2", mid="lightblue4", high="blue4") )
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(x = Offset, y = gluRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + ylab("Glutamte") + xlab("Offset")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lmer( gluRes.residuals ~ Offset + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(gluRes.residuals ~  Offset * age.x + Condition +Region +   (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## GABA VS Exponent
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(x = Exponent, y = gabaRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + ylab("GABA") + xlab("Exponent")+ theme(text = element_text(size = 30))
  
  summary(lmerTest::lmer(gabaRes.residuals ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
  summ(lm.model)
  
  
  lm.model <- lmer(gabaRes.residuals ~ Exponent *  age.x +Condition + Region+ (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  ## GABA VS Offset
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(x = Offset, y = gabaRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + ylab("GABA") + xlab("Offset")+ theme(text = element_text(size = 30))
  
  lm.model <- lmer(gabaRes.residuals  ~ Offset + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  lm.model <- lmer(gabaRes.residuals  ~ Offset * age.x + Condition +Region +   (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  ## Ratio VS Exponent
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(x = Exponent, y = ratioRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + ylab("Glu/GABA Ratio") +xlab("Exponent")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lmer(ratioRes.residuals  ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  
  lm.model <- lmer(ratioRes.residuals  ~ Exponent * age.x +Condition +Region+   (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  
  ## Ratio VS Offset 
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(x = Offset, y = ratioRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + ylab("Glu/GABA Ratio") + xlab("Offset")+ theme(text = element_text(size = 30))
  
  
  
  lm.model <- lmer(ratioRes.residuals  ~  Offset + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  
  lm.model <- lmer(ratioRes.residuals  ~ Offset *age.x +Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  p.adjust((0.0006905), method = "bonferroni", n = 2)
  
  
  ## gaba glu imbalance VS offset
  lunaize(ggplot(data = fooofMRS  %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) 
                 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = Offset, y = ratioRes.residuals, by =luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm",alpha=.8,size=1)) + scale_color_manual(values=c("gold3", "blue4")) + 
    ylab("Glu GABA Imbalance") + xlab("Offset")+ theme(text = element_text(size = 30))
  
  summary(lmerTest::lmer(ratioRes.residuals ~ Offset + age.x + GMrat + Condition + Region + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
  
  summary(lmerTest::lmer(ratioRes.residuals ~ Offset*age.x + GMrat + Region + Condition + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
  
  
  ## gaba glu imbalance VS exponent 
  lunaize(ggplot(data = fooofMRS 
                 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(y = ratioRes.residuals, x = Exponent, by =luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm",alpha=.8,size=1)) + scale_color_manual(values=c("gold3", "blue4"))+ 
    ylab("Glu GABA Imbalance") + xlab("Exponent")+ theme(text = element_text(size = 30))
  
  summary(lmerTest::lmer(ratioRes.residuals ~ Exponent + age.x + GMrat +Condition+ Region + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
  
  summary(lmerTest::lmer(ratioRes.residuals ~ Exponent*age.x + GMrat + Region + Condition + (1|subjID), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC")))
  
  p.adjust((0.0144 ), method = "bonferroni", n = 4)
  
  
  ## Johnson Neyman Plots & Condition == "eyesClosed"
  
  fit <- lmer(Offset ~ age*Glu.Cr.adj + Region + Condition + (1|luna),
              data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(visitno == 1))
  
  johnson_neyman(model = fit, pred = Glu.Cr.adj, modx = age, alpha = .05)
  
  
  ## Mediation
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
  
  
  mediationMatrix <- fooofMRS %>% dplyr::select(ratioRes.residuals, Offset, Condition, Region, luna, visitno, age.x) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
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

FOOOF_SNR_MRS <- function() {
  fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_20230111.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  SNR <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/frontalRegions_40minus20_SNR.csv')
  
  
  fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F)
  
  MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  idx <- which(!is.na(MRSregionsRes$GABA.Cr.adj) & !is.na(MRSregionsRes$Glu.Cr.adj))
  gabaglu.lm <- lm(Glu.Cr.adj ~ GABA.Cr.adj + Region, data = MRSregionsRes[idx,])
  MRSregionsRes$ratioRes.residuals <- NA
  MRSregionsRes[idx,]$ratioRes.residuals <- abs(gabaglu.lm$residuals)
  
  ## MERGE MRS AND FOOOF
  fooofMRS <- merge(MRSregionsRes, fooofPeaks_age, by = c("luna", "Region", "visitno")) 
  
  
  ## GABA glu imbalance measure VS Power   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 16, 'Adol', 'Adult')))
                 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = Offset, y = ratioRes.residuals, by =luna, color = Region)) + 
            geom_line(aes(group=interaction(luna,Region))) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Region, linetype=Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.9,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4"))
  
  
  lm.model <- lmer(Power ~ ratioRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Power ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Power ~ ratioRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  
  ## GLU VS AGE
  lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))
                 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = gluRes.residuals, by =luna, color = Region)) + 
            geom_line(aes(group=interaction(luna,Region))) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Region, linetype=Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.9,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4"))
  
  gam.model <-  gam(Glu.Cr.adj ~ s(age, k = 3) + GMrat + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("Glu.Cr.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  ## GABA VS AGE
  lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))
                 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = gabaRes.residuals, by =luna, color = Region)) + 
            geom_line(aes(group=interaction(luna,Region))) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Region, linetype=Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4"))
  
  gam.model <-  gam(GABA.Cr.adj ~ s(age, k = 3)+ GMrat + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("GABA.Cr.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  ## RATIO VS AGE
  lunaize(ggplot(data = MRSregionsRes %>% mutate(ageGroup = as.factor(ifelse(age <= 16, 'Adol', 'Adult')))
                 %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(x = age, y = ratioRes.residuals, by =luna, color = Region)) + 
            geom_line(aes(group=interaction(luna,Region))) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = Region, linetype=Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4"))
  
  
  gam.model <-  gam(Ratio.adj ~ s(age, k = 3) +GMrat + Region, data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("Ratio.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)+GMrat")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summary(model$gam)
  
  
  ## GLU VS Exponent
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent, x = gluRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.8) + 
            geom_point(alpha=.8) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(Exponent ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  
  
  lm.model <- lmer(Exponent ~ gluRes.residuals * Condition +Region+ age.x + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  ## GLU VS Offset + scale_color_gradient2(midpoint=20, low="gold2", mid="lightblue4", high="blue4") )
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Offset, x = gluRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.8) + 
            geom_point(alpha=.8) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(Offset ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  
  
  lm.model <- lmer(Offset ~ gluRes.residuals * Condition +Region + age.x +  (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## GABA VS Exponent
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 16, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent, x = gabaRes.residuals, by = luna, color = age.x)) + geom_line(aes(group=luna)) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = 1), method="lm") + facet_wrap(.~Condition) + scale_color_gradient2(midpoint=22, low="gold1", mid="blue2", high="blue4"))
  
  lm.model <- lmer(Exponent ~ gabaRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  lm.model <- lmer(Exponent ~ gabaRes.residuals * Condition + Region+ age.x + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  ## GABA VS Offset
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 16, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Offset, x = gabaRes.residuals, by = luna, color = age.x)) + geom_line(aes(group=luna)) + geom_point(alpha=.3) + 
            geom_smooth(aes(group = 1), method="lm") + facet_wrap(.~Condition) + scale_color_gradient2(midpoint=22, low="gold1", mid="blue2", high="blue4"))
  
  lm.model <- lmer(Offset ~ gabaRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  
  
  lm.model <- lmer(Offset ~ gabaRes.residuals * Condition +Region + age.x +  (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  ## Ratio VS Exponent
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent, x = ratioRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.8) + 
            geom_point(alpha=.8) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  lm.model <- lmer(Exponent ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  
  lm.model <- lmer(Exponent ~ ratioRes.residuals * Condition +Region+ age.x +  (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  
  ## Ratio VS Offset 
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Offset, x = ratioRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.8) + 
            geom_point(alpha=.8) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  
  
  lm.model <- lmer(Offset ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  
  lm.model <- lmer(Offset ~ ratioRes.residuals *Condition + Region + age.x +(1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  summ(lm.model)
  car::Anova(lm.model)
  
  
  
  p.adjust((0.0006905), method = "bonferroni", n = 2)
  
  
  
  ## Johnson Neyman Plots & Condition == "eyesClosed"
  
  fit <- lmer(Offset ~ age*Glu.Cr.adj + Region + Condition + (1|luna),
              data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(visitno == 1))
  
  johnson_neyman(model = fit, pred = Glu.Cr.adj, modx = age, alpha = .05)
  
  
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

peaksMRS <- function() {
  
  fooofPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalPeaks_allBands.csv')
  fooofPeaks_age <- merge(fooofPeaks, agefile, by = 'Subject')
  fooofPeaks_age <- fooofPeaks_age %>% separate(Subject,c("luna","vdate"), remove=F) 
  fooofPeaks_age$luna <- as.factor(fooofPeaks_age$luna)
  
  
  MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  idx <- which(!is.na(MRSregionsRes$GABA.Cr.adj) & !is.na(MRSregionsRes$Glu.Cr.adj))
  gabaglu.lm <- lm(Glu.Cr.adj ~ GABA.Cr.adj + Region, data = MRSregionsRes[idx,])
  MRSregionsRes$ratioRes.residuals <- NA
  MRSregionsRes[idx,]$ratioRes.residuals <- abs(gabaglu.lm$residuals)
  
  ## MERGE MRS AND FOOOF
  fooofMRS <- merge(MRSregionsRes, fooofPeaks_age, by = c("luna", "Region", "visitno")) 
 
  
  ## GABA glu imbalance measure VS Power   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Power, x = ratioRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Power ~ ratioRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Power ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Power ~ ratioRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
   ## GLU VS Power   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"), 
                 aes(y = Power, x = gluRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Power ~ gluRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Power ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Power ~ gluRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  ## GLU VS Bandwidth   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Bandwidth, x = gluRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Bandwidth ~ gluRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Bandwidth ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Bandwidth ~ gluRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ## GLU VS center frequency   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Center.Frequency, x = gluRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Center.Frequency ~ gluRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Center.Frequency ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Center.Frequency ~ gluRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ## GABA VS Power   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Power, x = gabaRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Power ~ gabaRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Power ~ gabaRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Power ~ gabaRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  ## GABA VS Bandwidth   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Bandwidth, x = gabaRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Bandwidth ~ gabaRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Bandwidth ~ gabaRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Bandwidth ~ gabaRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ## GABA VS center frequency   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Center.Frequency, x = gabaRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Center.Frequency ~ gabaRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Center.Frequency ~ gabaRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Center.Frequency ~ gabaRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ## Ratio VS Power   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Power, x = ratioRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Power ~ ratioRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Power ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Power ~ ratioRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ## Ratio VS Bandwidth   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Bandwidth, x = ratioRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Bandwidth ~ ratioRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Bandwidth ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Bandwidth ~ ratioRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Gamma"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ## Ratio VS center frequency   +  scale_color_manual(values=c("gold3", "blue4")))
  lunaize(ggplot(data = fooofMRS %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Center.Frequency, x = ratioRes.residuals, by = luna, color = Band))+ 
            geom_line(aes(group=interaction(luna,Region,Condition,Band)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = Band, color = Band), method="lm", alpha = 0.5))
  
  
  lm.model <- lmer(Center.Frequency ~ ratioRes.residuals + age.x + Condition + Region + Band + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Center.Frequency ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Center.Frequency ~ ratioRes.residuals * age.x + Condition +Region + (1|luna), data = fooofMRS %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(Band == "Beta"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  
  ## Mediation
  
  mediationMatrix <- fooofMRS %>% dplyr::select(ratioRes.residuals, Offset, Condition, Region, luna, visitno, age.x, Band) %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on offset (c)
  model.0 <- lmer(Power ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
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
  model.Y <- lmer(Power ~ gluRes.residuals + age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = 'gluRes.residuals', boot = FALSE, sims = 1000)
  (summary(results))
  
  
}

fooofvsMRSvsBehavior <- function() {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/03.BehaviorAcrossAge.R", envir = knitr::knit_global(), chdir = TRUE)
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv')
  Behavior <- Behavior_Sublevel_Maria() %>% separate(Subject,c("luna","vdate"), remove=F)
  
  
  fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/RegionalFOOOFMeasures_20230111.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  
  fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F)
  
  
  MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  
  idx <- which(!is.na(MRSregionsRes$GABA.Cr.adj) & !is.na(MRSregionsRes$Glu.Cr.adj))
  gabaglu.lm <- lm(Glu.Cr.adj ~ GABA.Cr.adj + Region, data = MRSregionsRes[idx,])
  MRSregionsRes$ratioRes.residuals <- NA
  MRSregionsRes$ratioRes.residuals_abs <- NA
  
  MRSregionsRes[idx,]$ratioRes.residuals <- (gabaglu.lm$residuals)
  MRSregionsRes[idx,]$ratioRes.residuals_abs <- abs(gabaglu.lm$residuals)
  
  MRSregionsRes$GluMinusGABA <- MRSregionsRes$gluRes.residuals - MRSregionsRes$gabaRes.residuals
  
  
  ## MERGE MRS AND FOOOF
  fooofMRS <- merge(MRSregionsRes, fooof_regional_age, by = c("luna", "Region", "visitno", "sex")) %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult')))
  
  fooofMRSbehavior <- merge(fooofMRS, Behavior, by = c("luna", "visitno"))
  
  write.csv(fooofMRSbehavior, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/fooofMRSbehavior_20230313.csv')
  
  ## accuracy VS Exponent
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(x = absBestError, y = Exponent, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  absBestError+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Exponent ~  absBestError* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.inverse <- lmer(Exponent ~  absBestError* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.inverse)
  
  AIC(lm.model,lm.inverse)
  
  fit <- lmer(Exponent ~ age.x*absBestError + Region + Condition + (1|luna),
              data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  
  johnson_neyman(model = fit, pred = absBestError, modx = age.x, alpha = .05)
  
  
  ## accuracy var VS Exponent
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  absBestError_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Exponent~ absBestError_sd * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Exponent ~  absBestError_sd* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(Exponent ~  absBestError_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  
  ## accuracy var VS offset
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Offset, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy Var") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Offset ~  absBestError_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Offset~  absBestError_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Offset ~  absBestError_sd* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(Offset ~  absBestError_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  ## Johnson Neyman Plots 
  
  fit <- lmer(absBestError_sd ~ age.x*Offset + Region + Condition + (1|luna),
              data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  
  johnson_neyman(model = fit, pred = Offset, modx = age.x, alpha = .05)
  
  
  ## accuracy  VS offset
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Offset, x = absBestError, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Offset~  absBestError+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Offset~  absBestError* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Offset ~  absBestError* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(Offset ~  absBestError* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  
  ## latency  VS offset
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 20, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Offset, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Latency") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Offset~  mgsLatency+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC") )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Offset ~ mgsLatency * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Offset ~  mgsLatency* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(Offset ~  mgsLatency* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  ## Johnson Neyman Plots 
  
  fit <- lmer(mgsLatency ~ age.x*Offset + Region + Condition + (1|luna),
              data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  
  johnson_neyman(model = fit, pred = Offset, modx = age.x, alpha = .05)
  
  
  
  
  ## latency var  VS offset
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Offset, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group =1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Latency Var") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Offset~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Offset ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Offset ~  mgsLatency_sd* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(Offset ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  
  ## latency  VS exponent
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 22, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Latency") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  mgsLatency+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Exponent~  mgsLatency* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Exponent ~  mgsLatency* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(Exponent ~  mgsLatency* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  ## Johnson Neyman Plots 
  
  fit <- lmer( Exponent~ inverseAge*mgsLatency + Region + Condition + (1|luna),
              data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  
  johnson_neyman(model = fit, pred = mgsLatency, modx = inverseAge, alpha = .05)
  
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(y = Offset, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + xlab("Latency") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))
  
  lm.model <- lmer(Exponent ~  mgsLatency + Condition + Region + (1|luna), data = fooofMRSbehavior%>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30")))  %>% filter(Region == "LDLPFC" | Region == "RDLPFC")%>% filter(ageGroup == "23-30"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  
  ## latency var  VS exponent
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group =1), method="lm", alpha = 0.8)) + xlab("Latency Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Exponent~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Exponent~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(Exponent ~  mgsLatency_sd* age.x + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(Exponent ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  ## acc VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = absBestError, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy") + ylab("ratioRes.residuals")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError+ inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError* inverseAge  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(ratioRes.residuals ~  absBestError* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(ratioRes.residuals ~  absBestError* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  ## acc var  VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy Var") + ylab("ratioRes.residuals")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError_sd+ inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError_sd* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(ratioRes.residuals ~  absBestError_sd* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(ratioRes.residuals ~  absBestError_sd* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  
  ## Lat VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Latency") + ylab("ratioRes.residuals")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( GluGABAimbalance~  mgsLatency+ inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( GluGABAimbalance~  mgsLatency* inverseAge  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(GluGABAimbalance ~  mgsLatency* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(GluGABAimbalance ~  mgsLatency* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  ## Johnson Neyman Plots 
  
  fit <- lmer( GluGABAimbalance~ inverseAge*mgsLatency + Region + Condition + (1|luna),
               data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  
  johnson_neyman(model = fit, pred = mgsLatency, modx = inverseAge, alpha = .05)
  
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(y = GluGABAimbalance, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
    geom_point(alpha=.5) + geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + xlab("Latency") + ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))
  
  ## Lat Var VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) +xlab("Latency Var") + ylab("ratioRes.residuals")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(ratioRes.residuals ~  mgsLatency_sd+ inverseAge  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  mgsLatency_sd* inverseAge  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(ratioRes.residuals ~  mgsLatency_sd* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(ratioRes.residuals ~  mgsLatency_sd* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  
  ## acc VS ratio
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = absBestError, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError+ age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError* age.x  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(ratioRes.residuals ~  absBestError* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(ratioRes.residuals ~  absBestError* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  ## acc var  VS ratio
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError_sd+ age.xage.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError_sd* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(ratioRes.residuals ~  absBestError_sd* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(ratioRes.residuals ~  absBestError_sd* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  
  ## Lat VS ratio
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( ratioRes.residuals~  mgsLatency+ age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  mgsLatency* age.x  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(ratioRes.residuals ~  mgsLatency* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(ratioRes.residuals ~  mgsLatency* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), aes(y = ratioRes.residuals, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
    geom_point(alpha=.5) + geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))
  
  ## Lat Var VS ratio
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = ratioRes.residuals, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) +xlab("Latency Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(ratioRes.residuals ~  mgsLatency_sd+ age.x  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  mgsLatency_sd* age.x  + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  lm.model <- lmer(ratioRes.residuals ~  mgsLatency_sd* age.x + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  lm.inverse <- lmer(ratioRes.residuals ~  mgsLatency_sd* inverseAge + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  AIC(lm.model,lm.inverse)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## regress age out of the measures 
  
  idx <- which(!is.na(fooofMRSbehavior$Exponent) & !is.na(fooofMRSbehavior$age.x))
  Exponent.lm <- lm(Exponent ~ inverseAge + Region, data = fooofMRSbehavior[idx,])  
  fooofMRSbehavior$Exponent_res <- NA
  fooofMRSbehavior[idx,]$Exponent_res <- (Exponent.lm$residuals)
  
  idx <- which(!is.na(fooofMRSbehavior$Offset) & !is.na(fooofMRSbehavior$age.x))
  Offset.lm <- lm(Offset ~ inverseAge + Region, data = fooofMRSbehavior[idx,])  
  fooofMRSbehavior$Offset_res <- NA
  fooofMRSbehavior[idx,]$Offset_res <- (Offset.lm$residuals)


  ## accuracy  VS exponent
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent_res, x = absBestError, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Accuracy") + ylab("Exponent age Res")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lmer(Exponent ~  absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Exponent ~ absBestError * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
 
  
  
  ## lat  VS exponent
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) %>% filter(Region == "LDLPFC" | Region == "RDLPFC"), 
                 aes(y = Exponent_res, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Latency") + ylab("Exponent age Res")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lmer(Exponent ~  mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Exponent ~ mgsLatency * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  ## Mediation
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
  
  idx <- which(!is.na(fooofMRSbehavior$Exponent) & !is.na(fooofMRSbehavior$Condition))
 exponent.lm <- lm(Exponent ~ Condition, data = fooofMRSbehavior[idx,])
 fooofMRSbehavior$exponent.residuals <- NA

 fooofMRSbehavior[idx,]$exponent.residuals <- (exponent.lm$residuals)
  
  mediationMatrix <- fooofMRSbehavior %>% mutate(ageGroup = as.factor(ifelse(age.x <= 20, 'Adol', 'Adult'))) %>% dplyr::select(GluGABAimbalance, gluRes.residuals, Exponent, Condition, Region, luna, visitno, inverseAge, ageGroup, GMrat) %>% filter(Region == "LDLPFC" | Region == "RDLPFC") %>% filter(ageGroup == "Adult")
  
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  

  #the effect of exponent on imbalance 
  model.0 <- lme4::lmer(GluGABAimbalance ~ inverseAge + Region + GMrat+(1|luna), data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  summ(model.0)
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  #the effect of exponent on latency 
  model.M <- lme4::lmer(gluRes.residuals ~ inverseAge + Region + (1|luna), data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  summ(model.M)
  
  #the effect of imbalance on latency
  model.Y <- lme4::lmer(GluGABAimbalance ~ gluRes.residuals+ inverseAge +GMrat+ Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = "inverseAge", mediator = "gluRes.residuals", boot = FALSE, sims = 1000)
  (summary(results))
  
 
}