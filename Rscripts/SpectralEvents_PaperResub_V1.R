
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
library(lme4)
library(mgcViz)
library(tidymv)

behavior <- function() {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/03.BehaviorAcrossAge.R", envir = knitr::knit_global(), chdir = TRUE)
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv')
  Behavior <- Behavior_Sublevel_Maria()
  trialLevelBehavior <- Behavior_TrialLevel_Maria()
  
  behaviorAge <-Behavior %>% filter(visitno == 1)
  
  ### Significant Periods of Growth 
  
  # BEST SACCADE
  gam.model <- gam(absBestError ~ s(age), data = behaviorAge)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  lunaize(gam_growthrate_plot(behaviorAge, gam.model, gam.growthrate, agevar = 'age', yvar = 'absBestError', draw_points = T, xplotname = "Age", yplotname = "Mean MGS Accuracy (degs)"))
  
  # BEST SACCADE VARIABILITY
  gam.model <- gam(absBestError_sd ~ s(age), data = behaviorAge)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  lunaize(gam_growthrate_plot(behaviorAge, gam.model, gam.growthrate, agevar = 'age', yvar = 'absBestError_sd', draw_points = F))
  
  # MGS LATENCY
  gam.model <- gam(mgsLatency ~ s(age), data = behaviorAge)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  lunaize(gam_growthrate_plot(behaviorAge, gam.model, gam.growthrate, agevar = 'age', yvar = 'mgsLatency', draw_points = F))
  
  # MGS LATENCY VARIABILITY
  gam.model <- gam(mgsLatency_sd ~ s(age), data = behaviorAge)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  lunaize(gam_growthrate_plot(behaviorAge, gam.model, gam.growthrate, agevar = 'age', yvar = 'mgsLatency_sd', draw_points = F))
  
  
  
  ## VGS
  (lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, 
                                                    y = vgsLatency)) + geom_point() + stat_smooth(method = "gam" )))
  
  # VGS LATENCY
  gam.model <- gam(vgsLatency ~ s(age), data = behaviorAge)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  lunaize(gam_growthrate_plot(behaviorAge, gam.model, gam.growthrate, agevar = 'age', yvar = 'vgsLatency', draw_points = F))
  
  # VGS LATENCY VARIABILITY
  gam.model <- gam(vgsLatency_sd ~ s(age), data = behaviorAge)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  lunaize(gam_growthrate_plot(behaviorAge, gam.model, gam.growthrate, agevar = 'age', yvar = 'vgsLatency_sd', draw_points = F))
  
}

spectralEvents_WholeBrain <- function () {
  
  #### Spectral Events
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/02.FrequencyAcrossAge.R", envir = knitr::knit_global(), chdir = TRUE)
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv')
  
  channelLocations <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')
  
  ## Load in Delay
  GammaDelay <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
  
  
  GammaDelay_Age <- merge(GammaDelay, agefile, by = "Subject")
  GammaDelay_Age$log_Gamma_Power <- log1p(GammaDelay_Age$Gamma_Trial_Power)
  GammaDelay_Age$inverseAge <- 1/GammaDelay_Age$age
  GammaDelay_Age_Channel <- merge(GammaDelay_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  GammaDelay_Age_Channel$Task <- 'Delay'
  
  
  # outlier detection 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  GammaDelay_Age_Channel_new <- GammaDelay_Age_Channel
  
  cols = names(GammaDelay_Age_Channel[c(4:6,11)])
  for (col in cols) {
    
    GammaDelay_Age_Channel_group <- GammaDelay_Age_Channel %>% group_by(Subject) %>% group_by(Channel)
    indx <- outliers(GammaDelay_Age_Channel_group[[col]])
    
    GammaDelay_Age_Channel_new[[col]] <- Map(replace, GammaDelay_Age_Channel_new[[col]], indx, NA)
    GammaDelay_Age_Channel_new[[col]] <- as.numeric(GammaDelay_Age_Channel_new[[col]])
    
  }  
  
  GammaDelay_Age_Channel_new <- GammaDelay_Age_Channel_new
  
  GammaDelay_Age_Channel_new_SubjectLevel_mean <- aggregate(.~Subject+Channel, data = GammaDelay_Age_Channel_new[c(1:6,11)], mean, na.rm = T)
  GammaDelay_Age_Channel_new_SubjectLevel_sd <- aggregate(.~Subject+Channel, data = GammaDelay_Age_Channel_new[c(1:6,11)], sd, na.rm = T)
  
  GammaDelay_Sublevel <- merge(GammaDelay_Age_Channel_new_SubjectLevel_mean, GammaDelay_Age_Channel_new_SubjectLevel_sd, by = c("Subject", "Channel"), suffixes = c("", "_sd"))
  names(GammaDelay_Sublevel) <- gsub('Gamma_','',names(GammaDelay_Sublevel))
  
  
  GammaDelay_Sublevel <- merge(GammaDelay_Sublevel, agefile, by = "Subject")
  GammaDelay_Sublevel$inverseAge <- 1/GammaDelay_Sublevel$age
  GammaDelay_Sublevel <- merge(GammaDelay_Sublevel, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  GammaDelay_Sublevel$Task <- 'Delay'
  
  
  
  # Load in Fix
  GammaFix <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_FIX.csv') 
  
  
  GammaFix_Age <- merge(GammaFix, agefile, by = "Subject")
  GammaFix_Age$log_Gamma_Power <- log1p(GammaFix_Age$Gamma_Trial_Power)
  GammaFix_Age$inverseAge <- 1/GammaFix_Age$age
  GammaFix_Age_Channel <- merge(GammaFix_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  GammaFix_Age_Channel$Task <- 'Fix'
  
  
  # outlier detection 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  GammaFix_Age_Channel_new <- GammaFix_Age_Channel
  
  cols = names(GammaFix_Age_Channel_new[c(4:6,11)])
  for (col in cols) {
    
    GammaFix_Age_Channel_group <- GammaFix_Age_Channel_new %>% group_by(Subject) %>% group_by(Channel)
    indx <- outliers(GammaFix_Age_Channel_group[[col]])
    
    GammaFix_Age_Channel_group[[col]] <- Map(replace, GammaFix_Age_Channel_group[[col]], indx, NA)
    GammaFix_Age_Channel_group[[col]] <- as.numeric(GammaFix_Age_Channel_group[[col]])
    
  }  
  
  GammaFix_Age_Channel_group <- GammaFix_Age_Channel_group
  
  
  
  GammaFix_Age_Channel_new_SubjectLevel_mean <- aggregate(.~Subject+Channel, data = GammaFix_Age_Channel_group[c(1:6,11)], mean, na.rm = T)
  GammaFix_Age_Channel_new_SubjectLevel_sd <- aggregate(.~Subject+Channel, data = GammaFix_Age_Channel_group[c(1:6,11)], sd, na.rm = T)
  
  GammaFix_Sublevel <- merge(GammaFix_Age_Channel_new_SubjectLevel_mean, GammaFix_Age_Channel_new_SubjectLevel_sd, by = c("Subject", "Channel"), suffixes = c("", "_sd"))
  names(GammaFix_Sublevel) <- gsub('Gamma_','',names(GammaFix_Sublevel))
  
  GammaFix_Sublevel <- merge(GammaFix_Sublevel, agefile, by = "Subject")
  GammaFix_Sublevel$inverseAge <- 1/GammaFix_Sublevel$age
  GammaFix_Sublevel <- merge(GammaFix_Sublevel, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  GammaFix_Sublevel$Task <- 'Fix'
  
  
  ## Combine delay and fix
  
  allGamma <- rbind(GammaFix_Sublevel, GammaDelay_Sublevel)
  
}

spectralEvents_individualChannels <- function () {
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv')
  agefile$inverseAge <- 1/agefile$age
  
  gammaAllChannels <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
  gammaAllChannels$Epoch <- 'Delay'
  
  gammaAllChannelsFix <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_FIX.csv') 
  gammaAllChannelsFix$Epoch <- 'Fix'
  
  gammaDelayFix <- rbind(gammaAllChannels, gammaAllChannelsFix)
  
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs <- chanLocs[c(-30),] #remove POz
  chanLocs$urchan <- (1:63)
  
  
  gammaDelayFix$urchan <- gammaDelayFix$Channel
  
  gammaDelayFix_locs <- merge(gammaDelayFix, chanLocs, by = "urchan")
  
  
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  gammaDelayFix_locs_new <- gammaDelayFix_locs
  
  
  outlier_cols <- names(gammaDelayFix_locs[c(5:7)])
  na_outliers <- function(x) ifelse(!outliers(x), x, NA)
  gammaDelayFix_locs_new <- gammaDelayFix_locs %>% group_by(Subject, Channel, Epoch) %>% mutate(across({outlier_cols}, na_outliers)) %>% ungroup 

  gammaDelayFix_locs_new <- merge(agefile, gammaDelayFix_locs_new, by = c("Subject"), all = T)
  
  
  aggregateTrials <- aggregate(cbind(Gamma_Trial_Power, Gamma_Event_Number, Gamma_Event_Duration, age, inverseAge, visitno) ~ Subject + labels + Epoch, data = gammaDelayFix_locs_new%>% filter(Trial < 97), mean)
  
  
  aggregateTrials_sd <- aggregate(cbind(Gamma_Trial_Power, Gamma_Event_Number, Gamma_Event_Duration) ~ Subject + labels + Epoch, data = gammaDelayFix_locs_new%>% filter(Trial < 97), sd)
  
  regionLevel <- merge(aggregateTrials, aggregateTrials_sd, by = c("Subject", "labels", "Epoch"), suffixes = c("", "_sd")) %>% filter(visitno == 1)
  
  # outlier on region level 
  regionLevel_new <- regionLevel

  outlier_cols <- names(regionLevel[c(4:6, 10:12)])
  na_outliers <- function(x) ifelse(!outliers(x), x, NA)
  regionLevel_new <- regionLevel %>% group_by(labels, Epoch) %>% mutate(across({outlier_cols}, na_outliers)) %>% ungroup 
  
  
  write.csv(regionLevel_new, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/allChannels_averagedOverTrials.csv')
  
  regionLevel_new<- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC_averagedAcrossChannelsandTrials.csv')
  
  # method = "gam", formula = y ~ splines::bs(x, 3)
  lunaize(ggplot(data = gammaDelayFix_locs_new , 
                 aes(x = age, y = Gamma_Event_Number, group = labels)) + geom_smooth(aes(color = labels), method = lm, formula = 'y ~ I(1/x)', alpha = 0))

  
  
  
  lunaize(ggplot(channelLevel_locs %>% mutate(isDLPFC = labels %in% c('F3','F5','F7','F4','F6','F8')), aes(x = age, y = Gamma_Trial_Power, color = isDLPFC)) + geom_smooth(aes(group = labels), method=lm, formula = 'y ~ I(1/x)',alpha=0.1,size=1)) + scale_color_manual(values=c("lightgray", "blue4"))
  
  inverse.model <- lm(Gamma_Event_Duration_sd ~ inverseAge, data = channelLevel_locs %>% filter(labels == 'T7'))
  anova(inverse.model)
  p.adjust(0.0004972, method = 'bonferroni', n = 9)
  
  
  lm.model <- lm(Gamma_Trial_Power ~ age, data = channelLevel_locs %>% filter(labels == 'P10'))
  gam.model <- gam(Gamma_Trial_Power ~ s(age, k =4), data = channelLevel_locs %>% filter(labels == 'P10'))
  
  AIC(lm.model, gam.model, inverse.model)
  
  
  
  ## add in behavior
  channelLevelBehavior <- merge(Behavior, channelLevel_locs, by = c('Subject', 'visitno', 'age', 'inverseAge'))
  
  lunaize(ggplot(data = channelLevelBehavior , 
                 aes(y = mgsLatency, x = Gamma_Event_Duration, group = labels)) + geom_smooth(aes(color = labels),  method = "lm", alpha = 0.1))
  
  lm.model <- lm(mgsLatency_sd ~ Gamma_Event_Duration_sd*age, data = channelLevelBehavior %>% filter(labels == 'F5'))
anova(lm.model)
  
  
  gam.model <- gam(absBestError ~ s(Gamma_Trial_Power, k=4) + s(age, k=4), data = channelLevelBehavior %>% filter(labels == 'F8'))
  summary(gam.model)
  
  p.adjust(0.008652, method = 'bonferroni', n = 9)
  
  
  
  
  
  AIC(lm.model, gam.model)
  
  
}


spectralEvents_DLPFCs <- function () {
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv')
  agefile$inverseAge <- 1/agefile$age
  
  gammaAllChannels <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
  gammaAllChannels$Epoch <- 'Delay'
  
  gammaAllChannelsFix <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_FIX.csv') 
  gammaAllChannelsFix$Epoch <- 'Fix'
  
  gammaDelayFix <- rbind(gammaAllChannels, gammaAllChannelsFix)
  
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs <- chanLocs[c(-30),] #remove POz
  chanLocs$urchan <- (1:63)
  
  
  gammaDelayFix$urchan <- gammaDelayFix$Channel
  gammaDelayFix$logGammaPower <- log1p(gammaDelayFix$Gamma_Trial_Power)
  
  gammaDelayFix_locs <- merge(gammaDelayFix, chanLocs, by = "urchan")
  
  RDLPFC <- filter(gammaDelayFix_locs, gammaDelayFix_locs$labels == 'F4' | gammaDelayFix_locs$labels == 'F6'| gammaDelayFix_locs$labels == 'F8')
  RDLPFC$Region <- 'R DLPFC'
  
  LDLPFC <- filter(gammaDelayFix_locs, gammaDelayFix_locs$labels == 'F3' | gammaDelayFix_locs$labels == 'F5'| gammaDelayFix_locs$labels == 'F7')
  LDLPFC$Region <- 'L DLPFC'
  
  bothDLPFCs <- rbind(RDLPFC, LDLPFC)
  
  DLPFCs_age <- merge(agefile, bothDLPFCs, by = c("Subject"), all = T)
  
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  DLPFCs_age_new <- DLPFCs_age
  
  cols = names(DLPFCs_age[c(10:12,14)])
  for (col in cols) {
    
    DLPFCs_age_grouped <- DLPFCs_age %>% group_by(Subject)%>% group_by(Channel) %>% group_by(Epoch)
    
    indx <- outliers(DLPFCs_age_grouped[[col]])
    
    DLPFCs_age_new[[col]] <- Map(replace, DLPFCs_age_new[[col]], indx, NA)
    DLPFCs_age_new[[col]] <- as.numeric(DLPFCs_age_new[[col]])
    
  } 
  
    write.csv(DLPFCs_age_new, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC_channels_TrialLevel.csv')

  aggregateTrials <- aggregate(cbind(Gamma_Trial_Power, logGammaPower, Gamma_Event_Number, Gamma_Event_Duration, age, inverseAge, visitno) ~ Subject + Region + Epoch, data = DLPFCs_age_new%>% filter(Trial < 97), mean)
  
  
  aggregateTrials_sd <- aggregate(cbind(Gamma_Trial_Power, logGammaPower, Gamma_Event_Number, Gamma_Event_Duration) ~ Subject + Region + Epoch, data = DLPFCs_age_new%>% filter(Trial < 97), sd)
  
  regionLevel <- merge(aggregateTrials, aggregateTrials_sd, by = c("Subject", "Region", "Epoch"), suffixes = c("", "_sd")) %>% filter(visitno == 1)

  # outlier on region level 
  regionLevel_new <- regionLevel
  
  cols = names(regionLevel[c(4:7, 11:14)])
  for (col in cols) {
    
    regionLevel_new <- regionLevel %>% group_by(Subject, Epoch) %>% mutate(col = ifelse(!outliers(col), col, NA)) %>% ungroup
  
    
  } 
  write.csv(regionLevel_new, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC_averagedAcrossChannelsandTrials.csv')
  
  regionLevel_new<- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC_averagedAcrossChannelsandTrials.csv')
  
  # creating delay-fix
  
  regionLevel_new <- regionLevel_new%>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult')))

  
  lunaize(ggplot(data = regionLevel_new %>% filter(age <=25)  , 
                 aes(x = age, y = Gamma_Trial_Power, color = Epoch, linetype=Region))  +
            geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Trial Power")+ theme(text = element_text(size = 30))
  
  #delay 
  regionLevel_new$Region <- as.factor(regionLevel_new$Region)
  
  gam.model <-  gam(Gamma_Trial_Power  ~ s(age, k=3) + Region  + Epoch  , data = regionLevel_new, random = ~(1|Subject))
  summary(gam.model)
  
  # interaction with region
  regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group
  model_formula <- as.formula("Gamma_Trial_Power ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch")
  # Note we keep fx = T for reliable p-values.
  model <- mgcv::gamm(model_formula,
                      random = list(Subject=~1),
                      data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
  summary(model$gam)
  
  ## interaction with epoch
  regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
  model_formula <- as.formula("Gamma_Trial_Power ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- mgcv::gamm(model_formula,
                      random = list(Subject=~1),
                      data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
  summary(model$gam)
  
  
  
  
  lunaize(ggplot(data = regionLevel_new %>% filter(age <= 25) , 
                 aes(x = age, y = Gamma_Trial_Power_sd, color = Epoch, linetype=Region))  +
            geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Trial Power Variability")+ theme(text = element_text(size = 30))
  
  #delay 
  regionLevel_new$Region <- as.factor(regionLevel_new$Region)
  
  gam.model <-  gam(Gamma_Trial_Power_sd  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new, random = ~(1|Subject))
  summary(gam.model)
  
  # interaction with region
  regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group
  model_formula <- as.formula("Gamma_Trial_Power_sd ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch")
  # Note we keep fx = T for reliable p-values.
  model <- mgcv::gamm(model_formula,
                      random = list(Subject=~1),
                      data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
  summary(model$gam)
  
  ## interaction with epoch
  regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
  model_formula <- as.formula("Gamma_Trial_Power_sd ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- mgcv::gamm(model_formula,
                      random = list(Subject=~1),
                      data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
  summary(model$gam)
  
  
  # method y ~I(1/x)
  lunaize(ggplot(data = regionLevel_new %>% filter(age <=25) , 
                 aes(x = age, y = Gamma_Event_Number, color = Epoch, linetype=Region))  +
            geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Number")+ theme(text = element_text(size = 30))
  
  #delay 
  regionLevel_new$Region <- as.factor(regionLevel_new$Region)
  
gam.model <-  gam(Gamma_Event_Number  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new, random = ~(1|Subject))
 summary(gam.model)
 
 # interaction with region
 regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group
 model_formula <- as.formula("Gamma_Event_Number ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
               random = list(Subject=~1),
               data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 ## interaction with epoch
 regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
 model_formula <- as.formula("Gamma_Event_Number ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
                     random = list(Subject=~1),
                     data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 
 # method = "gam", formula = y ~ splines::bs(x, 3)
 lunaize(ggplot(data = regionLevel_new %>% filter(age<=25) , 
                aes(x = age, y = Gamma_Event_Number_sd, color = Epoch, linetype=Region))  +
           geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Number Variability")+ theme(text = element_text(size = 30))
 
 #delay 
 regionLevel_new$Region <- as.factor(regionLevel_new$Region)
 
 gam.model <-  gam(Gamma_Event_Number_sd  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new, random = ~(1|Subject))
 summary(gam.model)
 
 # interaction with region
 regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group
 model_formula <- as.formula("Gamma_Event_Number_sd ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
                     random = list(Subject=~1),
                     data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 ## interaction with epoch
 regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
 model_formula <- as.formula("Gamma_Event_Number_sd ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
                     random = list(Subject=~1),
                     data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 
 
 lunaize(ggplot(data = regionLevel_new%>% filter(age<=25), 
                aes(x = age, y = Gamma_Event_Duration, color = Epoch, linetype=Region))  +
           geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Duration")+ theme(text = element_text(size = 30))
 
 #delay 
 regionLevel_new$Region <- as.factor(regionLevel_new$Region)
 
 gam.model <-  gam(Gamma_Event_Duration  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new, random = ~(1|Subject))
 summary(gam.model)
 
 # interaction with region
 regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group
 model_formula <- as.formula("Gamma_Event_Duration ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
                     random = list(Subject=~1),
                     data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 ## interaction with epoch
 regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
 model_formula <- as.formula("Gamma_Event_Duration ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
                     random = list(Subject=~1),
                     data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 
 
 
 
 lunaize(ggplot(data = regionLevel_new %>% filter(age<= 25) , 
                aes(x = age, y = Gamma_Event_Duration_sd, color = Epoch, linetype=Region))  +
           geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Duration Variability")+ theme(text = element_text(size = 30))
 
 #delay 
 regionLevel_new$Region <- as.factor(regionLevel_new$Region)
 
 gam.model <-  gam(Gamma_Event_Duration_sd  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new, random = ~(1|Subject))
 summary(gam.model)
 
 # interaction with region
 regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group
 model_formula <- as.formula("Gamma_Event_Duration_sd ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
                     random = list(Subject=~1),
                     data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 ## interaction with epoch
 regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
 model_formula <- as.formula("Gamma_Event_Duration_sd ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region")
 # Note we keep fx = T for reliable p-values.
 model <- mgcv::gamm(model_formula,
                     random = list(Subject=~1),
                     data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
 summary(model$gam)
 
 
 #sig periods of growth
 gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
 plist <- gam_growthrate_plot(regionLevel_new, gam.model, gam.growthrate, agevar = 'age', yvar = 'Gamma_Trial_Power', draw_points = T, xplotname = "Age", yplotname = "Gamma_Event_Duration_sd")
 plist$tile <- plist$tile +  scale_fill_gradient2(low = "red", mid = "white", high = "red")
 gam_growthrate_plot_combine(plist$ageplot, plist$tile)
 
 p.adjust(5.6e-05 , method = 'bonferroni', n = 12)
 
 
 
 delay <- regionLevel_new %>% subset(., Epoch == 'Delay') %>% dplyr::select(Subject, inverseAge, age, Epoch, Region, Gamma_Trial_Power,Gamma_Trial_Power_sd, Gamma_Event_Number,Gamma_Event_Number_sd, Gamma_Event_Duration,Gamma_Event_Duration_sd, visitno)
 fix <- regionLevel_new %>% subset(., Epoch == 'Fix') %>% dplyr::select(Subject, inverseAge, age, Epoch, Region, Gamma_Trial_Power,Gamma_Trial_Power_sd, Gamma_Event_Number,Gamma_Event_Number_sd, Gamma_Event_Duration,Gamma_Event_Duration_sd, visitno)
 delayFix <- merge(delay, fix, by = c("Subject", "age", "inverseAge", "Region", "visitno"), suffixes = c("_delay","_fix"))
 
 
  
  ## add in behavior
  channelLevelBehavior <- merge(Behavior, regionLevel_new, by = c('Subject', 'visitno', 'age', 'inverseAge'))
  
  EEGtrialLevelBehavior <- merge(trialLevelBehavior, DLPFCs_age_new, by = c('Subject', 'visitno', 'age', 'Trial'), all = T)
  
  write.csv(channelLevelBehavior, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC_averagedAcrossChannelsandTrials_Behavior.csv')
  write.csv(EEGtrialLevelBehavior, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC_channels_TrialLevel_Behavior.csv')
  
  gamma_long <- rbind(
    channelLevelBehavior %>% dplyr::select(Subject, inverseAge, age, Epoch, Region, Gamma_Trial_Power,Gamma_Trial_Power_sd, Gamma_Event_Number,Gamma_Event_Number_sd, Gamma_Event_Duration,Gamma_Event_Duration_sd,  Latency=vgsLatency, LatencySD=vgsLatency_sd, absBestError, absBestError_sd) %>% mutate(type='VGS') %>% subset(., Epoch == 'Fix'),
    channelLevelBehavior %>% dplyr::select(Subject, inverseAge, age, Epoch, Region, Gamma_Trial_Power, Gamma_Trial_Power_sd,Gamma_Event_Number, Gamma_Event_Number_sd, Gamma_Event_Duration,  Gamma_Event_Duration_sd, Latency=mgsLatency, LatencySD=mgsLatency_sd,absBestError, absBestError_sd) %>% mutate(type='MGS')%>% subset(., Epoch == 'Delay')
  )
  gamma_long$Subject <- as.factor(gamma_long$Subject) 
  gamma_long<-  gamma_long%>% separate(Subject,c("luna","vdate"), remove=F)
  
  
  gamma_long_delayMinusFix <- merge(Behavior, delayFix, by = c('Subject', 'age', 'inverseAge'))
  
  
 #color = ageGroup,
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = Latency, x = Gamma_Trial_Power, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power") +ylab("MGS Latency")+ theme(text = element_text(size = 30))
  
  library(jtools)
  
  lm.model <- lm(Latency ~ Gamma_Trial_Power + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Latency ~ Gamma_Trial_Power * age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = Latency, x = Gamma_Event_Number, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number") +ylab("MGS Latency")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(Latency ~ Gamma_Event_Number + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Latency ~ Gamma_Event_Number * age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = Latency, x = Gamma_Event_Duration, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Duration") +ylab("MGS Latency")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(Latency ~ Gamma_Event_Duration + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Latency ~ Gamma_Event_Duration * age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = LatencySD, x = Gamma_Trial_Power_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power Variability") +ylab("MGS Latency Variability")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(LatencySD ~ Gamma_Trial_Power_sd + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(LatencySD ~ Gamma_Trial_Power_sd * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = LatencySD, x = Gamma_Event_Number_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number Variability") +ylab("MGS Latency Variability")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(LatencySD ~ Gamma_Event_Number_sd + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(LatencySD ~ Gamma_Event_Number_sd * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = LatencySD, x = Gamma_Event_Duration_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number Variability") +ylab("MGS Latency Variability")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(LatencySD ~ Gamma_Event_Duration_sd + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(LatencySD ~ Gamma_Event_Duration_sd * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError, x = Gamma_Trial_Power, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power") +ylab("MGS Accuracy ")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(absBestError ~ Gamma_Trial_Power + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(absBestError ~ Gamma_Trial_Power * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError_sd, x = Gamma_Trial_Power_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power Variability") +ylab("MGS Accuracy Variability")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(absBestError_sd ~ Gamma_Trial_Power_sd + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(absBestError_sd ~ Gamma_Trial_Power_sd * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError, x = Gamma_Event_Number, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number") +ylab("MGS Accuracy ")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(absBestError ~ Gamma_Event_Number + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(absBestError ~ Gamma_Event_Number * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError_sd, x = Gamma_Event_Number_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number Variability") +ylab("MGS Accuracy Variability")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(absBestError_sd ~ Gamma_Event_Number_sd + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(absBestError_sd ~ Gamma_Event_Number_sd * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError, x = Gamma_Event_Duration, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Duration") +ylab("MGS Accuracy ")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(absBestError ~ Gamma_Event_Duration + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(absBestError ~ Gamma_Event_Duration * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError_sd, x = Gamma_Event_Duration_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Duration Variability") +ylab("MGS Accuracy Variability")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lm(absBestError_sd ~ Gamma_Event_Duration_sd + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(absBestError_sd ~ Gamma_Event_Duration_sd * age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  p.adjust(0.00532, method = 'bonferroni', n = 6)
  
  
  delayminusfix <- function () {

    delayFix$Gamma_Trial_Power_delayMinusFix <- delayFix$Gamma_Trial_Power_delay - delayFix$Gamma_Trial_Power_fix
    delayFix$Gamma_Trial_Power_sd_delayMinusFix <- delayFix$Gamma_Trial_Power_sd_delay - delayFix$Gamma_Trial_Power_sd_fix
    delayFix$Gamma_Event_Number_delayMinusFix <- delayFix$Gamma_Event_Number_delay - delayFix$Gamma_Event_Number_fix
    delayFix$Gamma_Event_Number_sd_delayMinusFix <- delayFix$Gamma_Event_Number_sd_delay - delayFix$Gamma_Event_Number_sd_fix
    delayFix$Gamma_Event_Duration_delayMinusFix <- delayFix$Gamma_Event_Duration_delay - delayFix$Gamma_Event_Duration_fix
    delayFix$Gamma_Event_Duration_sd_delayMinusFix <- delayFix$Gamma_Event_Duration_sd_delay - delayFix$Gamma_Event_Duration_sd_fix
    channelLevelBehavior <- merge(Behavior, delayFix, by = c('Subject', 'visitno', 'age', 'inverseAge'))
    
    lunaize(ggplot(data = channelLevelBehavior, 
                   aes(x = absBestError_sd, y = Gamma_Event_Number_sd_delayMinusFix, color = Region))  +
              geom_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("MGS Latency") +ylab("Event Number Delay-Fix")+ theme(text = element_text(size = 30))
    
    #delay 
    regionLevel_new$Region <- as.factor(regionLevel_new$Region)
    
    gam.model <-  gam(Gamma_Event_Number_delayMinusFix  ~ absBestError + s(age, k=3) + Region   , data = channelLevelBehavior, random = ~(1|Subject))
    summary(gam.model)
    
    # interaction with region
    regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group
    model_formula <- as.formula("Gamma_Event_Number_delayMinusFix ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch")
    # Note we keep fx = T for reliable p-values.
    model <- mgcv::gamm(model_formula,
                        random = list(Subject=~1),
                        data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
    summary(model$gam)
    
    ## interaction with epoch
    regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
    model_formula <- as.formula("Gamma_Event_Number_delayMinusFix ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region")
    # Note we keep fx = T for reliable p-values.
    model <- mgcv::gamm(model_formula,
                        random = list(Subject=~1),
                        data = regionLevel_new %>% filter(Region == "L DLPFC" | Region == "R DLPFC"))
    summary(model$gam)
    
    
    
    
    
    
    
  }
  
}

powerVsNumber <- function() {
  library(ggbeeswarm)
  
  regionLevel_new_ageGroups <- regionLevel_new  %>% mutate(ageGroup = as.factor(ifelse(age<15, '10-14',ifelse(age >= 15 & age<18, '15-17',ifelse(age  >= 18 & age<22, '18-21',ifelse(age >= 22 & age<26, '22-25', "over 25")))))) %>% aggregate(.~Subject+Epoch, ., mean)
  
  
  groupAvg <- regionLevel_new_ageGroups %>% group_by(ageGroup, Epoch) %>% 
    summarise(meanPower = mean(Gamma_Trial_Power, na.rm=T), 
              sdPower = sd(Gamma_Trial_Power, na.rm=T), 
              meanEvnum = mean(Gamma_Event_Number, na.rm=T), 
              sdEvnum = sd(Gamma_Event_Number, na.rm=T), 
              n = n()/2) %>% # divide by two since 2 hemis
    mutate(sePower = sdPower / sqrt(n), seEvnum = sdEvnum / sqrt(n))

  
  
  lunaize(ggplot(data = groupAvg, 
                 aes(x = ageGroup, y = meanPower, color = Epoch, fill=Epoch)) + 
            geom_bar(stat = "identity", position="dodge")) + 
            geom_errorbar(aes(ymin=meanPower-sePower, ymax=meanPower+sePower), position = position_dodge(0.9), width = 0.25, color="black")
  
  lunaize(ggplot(data = groupAvg, 
                 aes(x = ageGroup, y = meanEvnum, color = Epoch, fill=Epoch)) + 
            geom_bar(stat = "identity", position="dodge")) + 
    geom_errorbar(aes(ymin=meanEvnum-seEvnum, ymax=meanEvnum+seEvnum), position = position_dodge(0.9), width = 0.25, color="black") + coord_cartesian(ylim=c(2,2.8))
  
  regionLevel_new_ageGroups_wide <- regionLevel_new_ageGroups %>% tidyr::pivot_wider(id_cols = c('Subject','Region','ageGroup','age'), names_from='Epoch', values_from=c('Gamma_Trial_Power', 'Gamma_Event_Number', 'Gamma_Event_Duration')) %>%
    mutate(Gamma_Power_Diff = Gamma_Trial_Power_Delay - Gamma_Trial_Power_Fix,
           Gamma_Evnum_Diff = Gamma_Event_Number_Delay - Gamma_Event_Number_Fix) %>%
    mutate(zscorePower_Delay = zscore(Gamma_Trial_Power_Delay), 
           zscorePower_Fix = zscore(Gamma_Trial_Power_Fix), 
           zscoreNumber_Delay = zscore(Gamma_Event_Number_Delay),
           zscoreNumber_Fix = zscore(Gamma_Event_Number_Fix))
  
  
  groupAvgDiff <- regionLevel_new_ageGroups_wide %>% group_by(ageGroup) %>% 
    summarise(meanPowerDiff = mean(Gamma_Power_Diff, na.rm=T), 
              sdPowerDiff = sd(Gamma_Power_Diff, na.rm=T), 
              meanEvnumDiff = mean(Gamma_Evnum_Diff, na.rm=T), 
              sdEvnumDiff = sd(Gamma_Evnum_Diff, na.rm=T), 
              nDiff = n()/2) %>% # divide by two since 2 hemis
    mutate(sePowerDiff = sdPowerDiff / sqrt(nDiff), seEvnumDiff = sdEvnumDiff / sqrt(nDiff))
  
  #+ coord_cartesian(ylim=c(-.05, .3))
  lunaize(ggplot(data = groupAvgDiff %>% filter(ageGroup <5), 
                 aes(x = ageGroup, y = meanEvnumDiff)) + 
            geom_bar(stat = "identity", position="dodge", fill = "lightblue")) + 
    geom_errorbar(aes(ymin=meanEvnumDiff-seEvnumDiff, ymax=meanEvnumDiff+seEvnumDiff), position = position_dodge(0.9), width = 0.25, color="black") 
  
  write.csv(groupAvgDiff, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/powervsNumber_.csv')
  
 
  lunaize(ggplot(data =   rbind(groupAvgDiff %>% dplyr::select(ageGroup, meanPowerDiff, sePowerDiff) %>% mutate(type='Power') %>% rename(value=meanPowerDiff, se=sePowerDiff),
                                groupAvgDiff %>% dplyr::select(ageGroup, meanEvnumDiff , seEvnumDiff) %>% mutate(type='EventNumber') %>% rename(value=meanEvnumDiff, se=seEvnumDiff)), 
                 aes(x = ageGroup, y = value, color=type, fill=type)) + 
            geom_bar(stat = "identity", position="dodge")) + 
    geom_errorbar(aes(ymin=value-se, ymax=value+se), position = position_dodge(0.9), width = 0.25, color="black")
  
  
  ## Check ratios 
  avgRegions <- aggregate(.~Subject,regionLevel_new_ageGroups_wide[c(1,4:16)], mean)
  
  avgRegions$EventPowerRatio_delay <- avgRegions$Gamma_Event_Number_Delay/avgRegions$Gamma_Trial_Power_Delay
  avgRegions$EventPowerRatio_fix <- avgRegions$Gamma_Event_Number_Fix/avgRegions$Gamma_Trial_Power_Fix
  avgRegions$zscoredDiff_Delay <- avgRegions$zscorePower_Delay - avgRegions$zscoreNumber_Delay
  avgRegions$zscoredDiff_Fix <- avgRegions$zscorePower_Fix - avgRegions$zscoreNumber_Fix

  avgRegions_long <- avgRegions %>% pivot_longer(cols = c("zscoredDiff_Fix", "zscoredDiff_Delay"), names_to= "Epoch", values_to="zscoredDiff")
  
  
    lunaize(ggplot(data = avgRegions_long, 
                 aes(x = age, y = zscoredDiff, color=Epoch)) + stat_smooth(method="lm", formula=y~I(1/x)))
  
    gam.model <-  gam(zscoredDiff  ~ s(age, k=3) +s(age, k=3, by = as.factor(Epoch)) + Epoch  , data = avgRegions_long, random = ~(1|Subject))
    summary(gam.model)
    
    
    
  #modeling
  my.lm <- lm(Gamma_Event_Number ~ inverseAge*Epoch*Region, data = regionLevel_new)

    dt <- expand.grid(age = c(seq(14,26,by=1)), Epoch = c('Delay','Fix')) %>% mutate(Region = 'L DLPFC', inverseAge = 1/age)

  dt$xEvnum <- predict(my.lm, dt)
  
  myPower.lm <- lm(Gamma_Trial_Power ~ inverseAge*Epoch*Region, data = regionLevel_new)
  dt$xPower <- predict(myPower.lm, dt)
  
  
  lunaize(ggplot(data = dt  , 
                 aes(x = xEvnum, y = xPower, color = age, group=Epoch)) + geom_point(aes(shape=Epoch, color=age), size = 5)) + scale_size_manual(4) + geom_path(aes(group=age))
  
  
  #Dual axes
  combineDLPFCs <- regionLevel_new %>% dplyr::select(Subject, age,Epoch, Gamma_Trial_Power, Gamma_Event_Number) %>% aggregate(.~Subject+Epoch, data=., mean, na.rm=T)
  

 lunaize( ggplot(combineDLPFCs, aes(x = age, y = Gamma_Event_Number, linetype=Epoch)) +
    stat_smooth(method='lm', formula= y~I(1/x), aes(y = Gamma_Event_Number, color = "Event"), alpha = 0.1) +
    stat_smooth(method='lm', formula= y~I(1/x),aes(y = Gamma_Trial_Power/4, color = "Power"), alpha = 0.1) +
    scale_y_continuous(sec.axis = sec_axis(~.*4, name="Power")) +
    labs(x = "Age", y = "Event Number", color = "") +
    scale_color_manual(values = c("gold3", "blue4")))
  
  ## delay-fix 
 
 regionLevel_new_ageGroups <- delayFix %>% filter(age >=12 & age <32) %>% mutate(ageGroup = as.factor(ifelse(age >= 12 & age<15, '12-14',ifelse(age >= 15 & age<18, '15-17',ifelse(age  >= 18 & age<22, '18-21',ifelse(age >= 22 & age<26, '22-25', "over 25"))))))
 
 
 ageGroupTrialPowerDF <- regionLevel_new_ageGroups %>% dplyr::select(ageGroup,Gamma_Trial_Power_delayMinusFix) %>% aggregate(.~ageGroup, data=., mean, na.rm=T)%>% rename(SpectralMeasure = Gamma_Trial_Power_delayMinusFix)
 ageGroupTrialPowerDF$Measure = "Trial Power"
 
 ageGroupEventNumberDF <- regionLevel_new_ageGroups %>% dplyr::select(ageGroup,Gamma_Event_Number_delayMinusFix) %>% aggregate(.~ageGroup, data=., mean, na.rm=T) %>% rename(SpectralMeasure = Gamma_Event_Number_delayMinusFix)
 ageGroupEventNumberDF$Measure = "Event Number"
 
 allDF <- rbind(ageGroupTrialPowerDF, ageGroupEventNumberDF)
 
 
 #barplot
 lunaize(ggplot(data = allDF  , 
                aes(x = ageGroup, y = SpectralMeasure, fill=Measure)) + geom_bar(stat = "identity", position="dodge")) 
 
 
 
 + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25, color="black")
 
}

spectralVsTraditional <- function() {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  traditionalGammaPower_allChannels <- traditionalEEG_IndividualChannels()
  traditionalGammaPower_allChannels$Epoch <- as.factor(traditionalGammaPower_allChannels$Epoch)
  traditionalGammaPower_allChannels$labels <- traditionalGammaPower_allChannels$Channel
  
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs <- chanLocs[c(-30),] #remove POz
  chanLocs$urchan <- (1:63)
  
  traditionalGammaPower_allChannelsLocs <- merge(traditionalGammaPower_allChannels, chanLocs, by = "labels")
  
  lunaize(ggplot(traditionalGammaPower_allChannelsLocs %>% filter(Epoch == "Delay"), aes(x = -Y, y = X, fill = Power, z = Power, label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.03, low="white", mid="red", high="firebrick4")  + theme(text = element_text(size = 30)))
  
  
  RDLPFC <- filter(traditionalGammaPower_allChannels, traditionalGammaPower_allChannels$Channel == 'F4' | traditionalGammaPower_allChannels$Channel == 'F6'| traditionalGammaPower_allChannels$Channel == 'F8')
  RDLPFC$Region <- 'R DLPFC'
  
  LDLPFC <- filter(traditionalGammaPower_allChannels, traditionalGammaPower_allChannels$Channel == 'F3' | traditionalGammaPower_allChannels$Channel == 'F5'| traditionalGammaPower_allChannels$Channel == 'F7')
  LDLPFC$Region <- 'L DLPFC'
  
  bothDLPFCs <- rbind(RDLPFC, LDLPFC)
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  DLPFCs_age_new <- bothDLPFCs
    
    DLPFCs_age_new <- bothDLPFCs %>% group_by(Subject, Epoch) %>% mutate(Power = ifelse(!outliers(Power), Power, NA)) %>% ungroup
    
    
  AvgDLPFC <- DLPFCs_age_new %>% group_by(Subject, Epoch, age, visitno) %>% filter(Epoch == 'Delay'|Epoch == 'Fix') %>% summarise(meanPower_trad = mean(Power, na.rm=T), sdPower_trad = sd(Power, na.rm=T)) %>% tidyr::pivot_wider(id_cols = c('Subject','age'), names_from='Epoch', values_from=c('meanPower_trad', 'sdPower_trad'))
 
  
   lunaize(ggplot(data = AvgDLPFC, 
                 aes(x = age, y = meanPower_trad, color=Epoch)) +  geom_point(alpha=.3) +
             geom_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 4, fx = T),alpha=.2,size=.75)) 
   
   
   lunaize(ggplot(data = AvgDLPFC, 
                  aes(x = age, y = meanPower_trad, color=Epoch)) + stat_smooth(method="lm", formula=y~I(1/x)))
   
  
    gam.model <-  gam(meanPower_trad  ~ s(age, k=3) + Epoch, data = AvgDLPFC, random = ~(1|Subject))
   summary(gam.model)
  
 spectralTrad <- merge(AvgDLPFC,avgRegions_long, by = "Subject")
 
 spectralTrad$eventPowerPercentage_Delay <- (spectralTrad$Gamma_Event_Number_Delay*spectralTrad$Gamma_Trial_Power_Delay)/spectralTrad$meanPower_trad_Delay
  
 spectralTrad$eventPowerPercentage_Fix <- (spectralTrad$Gamma_Event_Number_Fix*spectralTrad$Gamma_Trial_Power_Fix)/spectralTrad$meanPower_trad_Fix
 
 spectralTrad$eventTradPowerRatio_Delay <- spectralTrad$Gamma_Event_Number_Delay/spectralTrad$meanPower_trad_Delay
 spectralTrad$eventTradPowerRatio_Fix <- spectralTrad$Gamma_Event_Number_Fix/spectralTrad$meanPower_trad_Fix
 
 
 
 spectralTrad_long <- spectralTrad %>% pivot_longer(cols = c("eventTradPowerRatio_Delay", "eventTradPowerRatio_Fix"), names_to= "Epoch2", values_to="eventPowerRatio")
  
   lunaize(ggplot(data = spectralTrad_long, 
                 aes(x = age.x, y = eventPowerRatio, color=Epoch2)) + stat_smooth(method="lm", formula=y~I(1/x)))
  
  gam.model <-  gam(eventPowerRatio  ~ s(age.x, k=3) +s(age.x, k=3, by = as.factor(Epoch2)) + Epoch2  , data = spectralTrad_long, random = ~(1|Subject))
  summary(gam.model)
  
  
  ## Mediation
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
  
  regionLevel_new<- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC_averagedAcrossChannelsandTrials.csv')
  
 
  
  regionLevel_new_wide <- regionLevel_new %>% tidyr::pivot_wider(id_cols = c('Subject','Region','age'), names_from='Epoch', values_from=c('Gamma_Trial_Power', 'Gamma_Event_Number', 'Gamma_Event_Duration')) 
  
   avgRegions <- aggregate(.~Subject,regionLevel_new_wide[c(1,3:9)], mean)
     
  spectralTrad_med <- merge(AvgDLPFC,avgRegions, by = "Subject")
  
  other_vars <- c("Gamma_Event_Number_Fix","meanPower_trad_Fix")
  spectralTrad_med$invage <- 1/spectralTrad_med$age.x
  mediationAnalysis_invage(spectralTrad_med, other_vars)
  mediationAnalysis(spectralTrad_med, other_vars)
  
  other_vars <- c("Gamma_Event_Number_Delay","meanPower_trad_Delay")
  spectralTrad_med$invage <- 1/spectralTrad_med$age.x
  mediationAnalysis_invage(spectralTrad_med, other_vars)
  mediationAnalysis(spectralTrad_med, other_vars)
  
  
  
  
  mediationMatrix <-spectralTrad_med 
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ age.x", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ age.x", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ age.x +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = other_vars[1], boot = FALSE, sims = 1000)
  
  
  
  
  
  
  
  
  #check for mediation inluding epoch in model
  AvgDLPFC_epochs <- DLPFCs_age_new %>% group_by(Subject, Epoch, age, visitno) %>% 
    filter(Epoch == 'Delay'|Epoch == 'Fix') %>%
    summarise(meanPower_trad = mean(Power, na.rm=T), 
              sdPower_trad = sd(Power, na.rm=T))
  
  regionLevel_AvgDLPFC <- regionLevel_new %>% group_by(Subject, Epoch) %>% 
    summarise(meanPower = mean(Gamma_Trial_Power, na.rm=T), 
              sdPower = sd(Gamma_Trial_Power, na.rm=T), 
              meanEvnum = mean(Gamma_Event_Number, na.rm=T), 
              sdEvnum = sd(Gamma_Event_Number, na.rm=T),
              meanDuration = mean(Gamma_Event_Duration, na.rm=T))
  
  spectralTrad_med_withEpoch <- merge(AvgDLPFC_epochs,regionLevel_AvgDLPFC, by = c("Subject", "Epoch"))
  
  other_vars <- c("meanEvnum","meanPower_trad","Epoch")
  
  spectralTrad_med_withEpoch$invage <- 1/spectralTrad_med_withEpoch$age

  mediationAnalysis_Epoch_invage(spectralTrad_med_withEpoch, other_vars)
  mediationAnalysis_Epoch(spectralTrad_med_withEpoch, other_vars)
  
}


correlationBetweenSpectralMeasures <- function() {
    library(corrplot)

  #grab only adolescence 
  
  regionLevel_new_avgDLPFC <- aggregate(.~Subject, data=regionLevel_new, mean )
  
  regionLevel_adol <- regionLevel_new %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(ageGroup=="Adol") %>% subset(., Epoch == 'Fix') %>% dplyr::select(Subject,Gamma_Trial_Power,Gamma_Trial_Power_sd, Gamma_Event_Number,Gamma_Event_Number_sd, Gamma_Event_Duration,Gamma_Event_Duration_sd) %>% aggregate(.~Subject, data=., mean )
    
  
  AdolCorrelations <- cor(regionLevel_adol[2:7], method="pearson", use = "complete.obs")
  
  
  corrplot(AdolCorrelations, method = "color", order = "alphabet", tl.col = "black",cl.pos = 'b', tl.pos = 'n' , cl.cex = 1,addCoef.col = 'black', number.cex=3,cl.ratio=.3, col= colorRampPalette(c("navy","white","firebrick"))(20))
  
  
  regionLevel_adult <- regionLevel_new %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) %>% filter(ageGroup=="Adult") %>% subset(., Epoch == 'Fix')%>% dplyr::select(Subject,Gamma_Trial_Power,Gamma_Trial_Power_sd, Gamma_Event_Number,Gamma_Event_Number_sd, Gamma_Event_Duration,Gamma_Event_Duration_sd) %>% aggregate(.~Subject, data=., mean )
  
  AdultCorrelations <- cor(regionLevel_adult[2:7], method="pearson", use = "complete.obs")
  

  corrplot(AdultCorrelations, method = "color", order = "alphabet", tl.col = "black",cl.pos = 'b', tl.pos = 'n' , addCoef.col = 'black', number.cex=3,cl.ratio=.3, col= colorRampPalette(c("navy","white","firebrick"))(20))  
  
  

 }


individualChannels_SpectralvsTraditional <- function() {
  
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  traditionalGammaPower_allChannels <- traditionalEEG_IndividualChannels()
  traditionalGammaPower_allChannels$Epoch <- as.factor(traditionalGammaPower_allChannels$Epoch)
  traditionalGammaPower_allChannels$labels <- traditionalGammaPower_allChannels$Channel
  
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs <- chanLocs[c(-30),] #remove POz
  chanLocs$urchan <- (1:63)
  
  traditionalGammaPower_allChannelsLocs <- merge(traditionalGammaPower_allChannels, chanLocs, by = "labels")
  
  grabOneChannel_trad <- traditionalGammaPower_allChannelsLocs %>% filter(Channel == "FT9")
  
  grabOneChannel_trad_wide <- grabOneChannel_trad %>% tidyr::pivot_wider(id_cols = c('Subject','labels','age'), names_from='Epoch', values_from=c('Power')) 
  
  # spectral events
  allChannels<- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/allChannels_averagedOverTrials.csv')
  
  grabOneChannel_spectral <- allChannels %>% filter(labels == "FT9")
  
  
  
  lunaize(ggplot(data = grabOneChannel_spectral %>% filter(age <=25)  , 
                 aes(x = age, y = Gamma_Event_Number, color = Epoch))  +
            geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Number")+ theme(text = element_text(size = 30))
  
  gam.model <-  gamm(Gamma_Event_Number  ~ s(age, k=3)  + Epoch  , data = grabOneChannel_spectral%>% filter(age <=25)  , random=list(Subject=~1))
  summary(gam.model$gam)
  
  ## interaction with epoch
  grabOneChannel_spectral$oEpoch <- ordered(grabOneChannel_spectral$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
  model_formula <- as.formula("Gamma_Event_Number ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- mgcv::gamm(model_formula,
                      random = list(Subject=~1),
                      data = grabOneChannel_spectral%>% filter(age <=25))
  summary(model$gam)
  
  
  
  
  
  lunaize(ggplot(data = grabOneChannel_spectral %>% filter(age <=25)  , 
                 aes(x = age, y = Gamma_Trial_Power, color = Epoch))  +
            geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Trial Power")+ theme(text = element_text(size = 30))
  
  gam.model <-  gamm(Gamma_Trial_Power  ~ s(age, k=3)  + Epoch  , data = grabOneChannel_spectral%>% filter(age <=25)  , random=list(Subject=~1))
  summary(gam.modelgam)
  
  ## interaction with epoch
  grabOneChannel_spectral$oEpoch <- ordered(grabOneChannel_spectral$Epoch, levels = c("Fix","Delay")) # Fix will be the reference group
  model_formula <- as.formula("Gamma_Trial_Power ~ oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- mgcv::gamm(model_formula,
                      random = list(Subject=~1),
                      data = grabOneChannel_spectral%>% filter(age <=25))
  summary(model$gam)
  
  
  ## Mediation on individual channels
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/06.Mediation.R", envir = knitr::knit_global(), chdir = TRUE)
  
  
  grabOneChannel_spectral_wide <- grabOneChannel_spectral %>% tidyr::pivot_wider(id_cols = c('Subject','labels','age'), names_from='Epoch', values_from=c('Gamma_Trial_Power', 'Gamma_Event_Number', 'Gamma_Event_Duration')) 
  
  
  spectralTrad_med <- merge(grabOneChannel_spectral_wide,grabOneChannel_trad_wide, by = "Subject")
  
  other_vars <- c("Gamma_Trial_Power_Fix","Fix")
  spectralTrad_med$invage <- 1/spectralTrad_med$age.x
  mediationAnalysis_invage(spectralTrad_med, other_vars)
  mediationAnalysis(spectralTrad_med, other_vars)
  
  other_vars <- c("Gamma_Event_Number_Delay","Delay")
  spectralTrad_med$invage <- 1/spectralTrad_med$age.x
  mediationAnalysis_invage(spectralTrad_med, other_vars)
  mediationAnalysis(spectralTrad_med, other_vars)
  
  
  
  
  
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ age.x", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ age.x", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ age.x +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = other_vars[1], boot = FALSE, sims = 1000)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}









