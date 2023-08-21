
# 1.0 Libraries ----

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
library(jtools)


# 2.0 Define Functions ----

outliers <- function(x) {
  
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
  
  
}


# 3.0 Load in Behavior ----
  
  sys.source("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/spectral_event_scripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  sys.source("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/spectral_event_scripts/03.BehaviorAcrossAge.R", envir = knitr::knit_global(), chdir = TRUE)
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv')
  Behavior <- Behavior_Sublevel_Maria()
  trialLevelBehavior <- Behavior_TrialLevel_Maria()
  
  behaviorAge <-Behavior %>% filter(visitno == 1)
  
  
  
# 4.0 Behavior vs Age ----
  
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
  

# 5.0 Load DLPFC Spectral Event Data ----
  
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
  
  ## 5.1 Define DLPFC Channels ----
  
  RDLPFC <- filter(gammaDelayFix_locs, gammaDelayFix_locs$labels == 'F4' | gammaDelayFix_locs$labels == 'F6'| gammaDelayFix_locs$labels == 'F8')
  RDLPFC$Region <- 'R DLPFC'
  
  LDLPFC <- filter(gammaDelayFix_locs, gammaDelayFix_locs$labels == 'F3' | gammaDelayFix_locs$labels == 'F5'| gammaDelayFix_locs$labels == 'F7')
  LDLPFC$Region <- 'L DLPFC'
  
  bothDLPFCs <- rbind(RDLPFC, LDLPFC)
  
  DLPFCs_age <- merge(agefile, bothDLPFCs, by = c("Subject"), all = T)
  

  ## 5.2 Outlier Detection ----
  
  DLPFCs_age_new <- DLPFCs_age %>% group_by(Subject, Channel, Epoch) %>% mutate(Gamma_Trial_Power = ifelse(!outliers(Gamma_Trial_Power), Gamma_Trial_Power, NA)) %>% ungroup
  
  DLPFCs_age_new <- DLPFCs_age_new %>% group_by(Subject, Channel, Epoch) %>% mutate(Gamma_Event_Number = ifelse(!outliers(Gamma_Event_Number), Gamma_Event_Number, NA)) %>% ungroup

  DLPFCs_age_new <- DLPFCs_age_new %>% group_by(Subject, Channel, Epoch) %>% mutate(Gamma_Event_Duration = ifelse(!outliers(Gamma_Event_Duration), Gamma_Event_Duration, NA)) %>% ungroup
  
  DLPFCs_age_new <- DLPFCs_age_new %>% group_by(Subject, Channel, Epoch) %>% mutate(logGammaPower = ifelse(!outliers(logGammaPower), logGammaPower, NA)) %>% ungroup
  
  ## 5.3 Average Trials ----

  aggregateTrials <- aggregate(cbind(Gamma_Trial_Power, logGammaPower, Gamma_Event_Number, Gamma_Event_Duration, age, inverseAge, visitno) ~ Subject + Region + Epoch, data = DLPFCs_age_new%>% filter(Trial < 97), mean)
  
  
  aggregateTrials_sd <- aggregate(cbind(Gamma_Trial_Power, logGammaPower, Gamma_Event_Number, Gamma_Event_Duration) ~ Subject + Region + Epoch, data = DLPFCs_age_new%>% filter(Trial < 97), sd)
  
  regionLevel <- merge(aggregateTrials, aggregateTrials_sd, by = c("Subject", "Region", "Epoch"), suffixes = c("", "_sd")) %>% filter(visitno == 1)

  ## 5.4 Subject level outlier detection ----
  
  regionLevel_new <- regionLevel

  outlier_cols <- names(regionLevel[c(4:7, 11:14)])
  na_outliers <- function(x) ifelse(!outliers(x), x, NA)
  regionLevel_new <- regionLevel %>% group_by(Region, Epoch) %>% mutate(across({outlier_cols}, na_outliers)) %>% ungroup 
  regionLevel_new$Region <- as.factor(regionLevel_new$Region)


# 6.0 Spectral Event vs Age ----

  ## 6.1 Trial Power ----
  
  lunaize(ggplot(data = regionLevel_new %>% filter(age <= 25) , 
                 aes(x = age, y = Gamma_Trial_Power, linetype = Epoch, color=Region))  +
            geom_smooth(aes(group = interaction(Region, Epoch)), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Trial Power")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  gam.model <-  gamm(Gamma_Trial_Power  ~ s(age, k=3) + Region  + Epoch  , data = regionLevel_new,  random=list(Subject=~1))
  summary(gam.model$gam)
  
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
  
  
  ## 6.2 Trial Power Variability ----
  
  lunaize(ggplot(data = regionLevel_new %>% filter(age <= 25) , 
                 aes(x = age, y = Gamma_Trial_Power_sd, color = Epoch, linetype=Region)) +
            geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Trial Power Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
 
  gam.model <-  gamm(Gamma_Trial_Power_sd  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new,  random=list(Subject=~1))
  summary(gam.model$gam)
  
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
  
  ## 6.3 Event Number ----
  
  # method y ~I(1/x)
  lunaize(ggplot(data = regionLevel_new %>% filter(age <=25) , 
                 aes(x = age, y = Gamma_Event_Number, color = Epoch, linetype=Region))  +
            geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Number")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  

gam.model <-  gamm(Gamma_Event_Number  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new,  random=list(Subject=~1))
 summary(gam.model$gam)
 
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
 
 ## 6.4 Event Number Variability ----
 
 # method = "gam", formula = y ~ splines::bs(x, 3)
 lunaize(ggplot(data = regionLevel_new %>% filter(age<=25) , 
                aes(x = age, y = Gamma_Event_Number_sd, color = Epoch, linetype=Region))  +
           geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Number Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
 
 gam.model <-  gamm(Gamma_Event_Number_sd  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new, random=list(Subject=~1))
 summary(gam.model$gam)
 
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
 
 ## 6.5 Event Duration ----
 
 lunaize(ggplot(data = regionLevel_new%>% filter(age<=25), 
                aes(x = age, y = Gamma_Event_Duration, color = Epoch, linetype=Region))  +
           geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Duration")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
 
 gam.model <-  gamm(Gamma_Event_Duration  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new, random=list(Subject=~1))
 summary(gam.model$gam)
 
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
 
 
## 6.6 Event Duration Variability ----
 
 lunaize(ggplot(data = regionLevel_new %>% filter(age<= 25) , 
                aes(x = age, y = Gamma_Event_Duration_sd, color = Epoch, linetype=Region))  +
           geom_smooth(aes(group = Epoch),method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.2,size=.75)) + xlab("Age") +ylab("Gamma Event Duration Variability")+ theme(text = element_text(size = 30))
 
 
 gam.model <-  gamm(Gamma_Event_Duration_sd  ~ s(age, k=3) + Region + Epoch  , data = regionLevel_new,random=list(Subject=~1))
 summary(gam.model$gam)
 
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

 
 # 7.0 Merge EEG and Behavior ----
 
  ## add in behavior
  channelLevelBehavior <- merge(Behavior, regionLevel_new, by = c('Subject', 'visitno', 'age', 'inverseAge'))

  
  gamma_long <- rbind(
    channelLevelBehavior %>% dplyr::select(Subject, inverseAge, age, Epoch, Region, Gamma_Trial_Power,Gamma_Trial_Power_sd, Gamma_Event_Number,Gamma_Event_Number_sd, Gamma_Event_Duration,Gamma_Event_Duration_sd,  Latency=vgsLatency, LatencySD=vgsLatency_sd, absBestError, absBestError_sd) %>% mutate(type='VGS') %>% subset(., Epoch == 'Fix'),
    channelLevelBehavior %>% dplyr::select(Subject, inverseAge, age, Epoch, Region, Gamma_Trial_Power, Gamma_Trial_Power_sd,Gamma_Event_Number, Gamma_Event_Number_sd, Gamma_Event_Duration,  Gamma_Event_Duration_sd, Latency=mgsLatency, LatencySD=mgsLatency_sd,absBestError, absBestError_sd) %>% mutate(type='MGS')%>% subset(., Epoch == 'Delay')
  )
  gamma_long$Subject <- as.factor(gamma_long$Subject) 
  gamma_long<-  gamma_long%>% separate(Subject,c("luna","vdate"), remove=F)
  

  ## 7.1 Trial Power vs Latency ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = Latency, x = Gamma_Trial_Power, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power") +ylab("MGS Latency")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

  lm.model <- lm( Gamma_Trial_Power ~ Latency + age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm( Gamma_Trial_Power ~  Latency* age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## 7.2 Event Number vs Latency ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = Latency, x = Gamma_Event_Number, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number") +ylab("MGS Latency")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Event_Number ~  Latency+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm( Gamma_Event_Number~  Latency* age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.3 Event Duration vs Latency ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = Latency, x = Gamma_Event_Duration, color = Region)) + geom_smooth(method = "lm", alpha = 0.1)) + xlab("Event Duration") +ylab("MGS Latency")+ theme(text = element_text(size = 30)) + theme(legend.position='none')
  
  
  lm.model <- lm( Gamma_Event_Duration~  Latency+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm( Gamma_Event_Duration~  Latency* age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.4 Trial Power Var vs Latency Var ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = LatencySD, x = Gamma_Trial_Power_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power Variability") +ylab("MGS Latency Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Trial_Power_sd ~  LatencySD+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Trial_Power_sd ~  LatencySD* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.5 Event Number Var vs Latency Var ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = LatencySD, x = Gamma_Event_Number_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number Variability") +ylab("MGS Latency Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Event_Number_sd ~  LatencySD+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Event_Number_sd ~  LatencySD* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.6 Event Duration Var vs Latency Var ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = LatencySD, x = Gamma_Event_Duration_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number Variability") +ylab("MGS Latency Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Event_Duration_sd ~  LatencySD+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Event_Duration_sd ~  LatencySD* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.7 Trial Power vs Accuracy ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError, x = Gamma_Trial_Power, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power") +ylab("MGS Accuracy ")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Trial_Power ~  absBestError+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Trial_Power ~  absBestError* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.8 Trial Power Var vs Accuracy Var ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError_sd, x = Gamma_Trial_Power_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Trial Power Variability") +ylab("MGS Accuracy Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Trial_Power_sd ~  absBestError_sd+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Trial_Power_sd ~  absBestError_sd* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.9 Event Number vs Accuracy ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError, x = Gamma_Event_Number, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number") +ylab("MGS Accuracy ")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Event_Number ~  absBestError+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Event_Number ~  absBestError* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.10 Event Number Var vs Accuracy Var ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError_sd, x = Gamma_Event_Number_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Number Variability") +ylab("MGS Accuracy Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Event_Number_sd ~  absBestError_sd+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Event_Number_sd ~  absBestError_sd* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.11 Event Duration vs Accuracy ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError, x = Gamma_Event_Duration, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Duration") +ylab("MGS Accuracy ")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Event_Duration ~  absBestError+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Event_Duration ~  absBestError* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## 7.12 Event Duration Var vs Accuracy Var ----
  
  lunaize(ggplot(data = gamma_long  %>% filter(Epoch == 'Delay') %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))) ,
                 aes(y = absBestError_sd, x = Gamma_Event_Duration_sd, color = Region)) + geom_smooth(method = "lm", alpha = 0.1))+ xlab("Event Duration Variability") +ylab("MGS Accuracy Variability")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lm(Gamma_Event_Duration_sd ~  absBestError_sd+ age + Region , data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lm(Gamma_Event_Duration_sd~  absBestError_sd* age + Region, data = gamma_long %>% filter(Epoch == 'Delay'))
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  p.adjust(0.00532, method = 'bonferroni', n = 6)
  
  
# 8.0 Traditional EEG ----
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  traditionalGammaPower_allChannels <- traditionalEEG_IndividualChannels()
  traditionalGammaPower_allChannels$Epoch <- as.factor(traditionalGammaPower_allChannels$Epoch)
  traditionalGammaPower_allChannels$labels <- traditionalGammaPower_allChannels$Channel
  
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs <- chanLocs[c(-30),] #remove POz
  chanLocs$urchan <- (1:63)
  
  ## 8.1 Plot across Channels ----
  
  traditionalGammaPower_allChannelsLocs <- merge(traditionalGammaPower_allChannels, chanLocs, by = "labels")
  
  lunaize(ggplot(traditionalGammaPower_allChannelsLocs %>% filter(Epoch == "Delay"), aes(x = -Y, y = X, fill = Power, z = Power, label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.03, low="white", mid="red", high="firebrick4")  + theme(text = element_text(size = 30)))
  
  ## 8.2 Average DLPFC Electrodes ----
  
  RDLPFC <- filter(traditionalGammaPower_allChannels, traditionalGammaPower_allChannels$Channel == 'F4' | traditionalGammaPower_allChannels$Channel == 'F6'| traditionalGammaPower_allChannels$Channel == 'F8')
  RDLPFC$Region <- 'R DLPFC'
  
  LDLPFC <- filter(traditionalGammaPower_allChannels, traditionalGammaPower_allChannels$Channel == 'F3' | traditionalGammaPower_allChannels$Channel == 'F5'| traditionalGammaPower_allChannels$Channel == 'F7')
  LDLPFC$Region <- 'L DLPFC'
  
  bothDLPFCs <- rbind(RDLPFC, LDLPFC)
  
  ## 8.3 Outlier Detection ----
  
  DLPFCs_age_new <- bothDLPFCs
    
  DLPFCs_age_new <- bothDLPFCs %>% group_by(Subject, Epoch) %>% mutate(Power = ifelse(!outliers(Power), Power, NA)) %>% ungroup
    
  AvgDLPFC <- DLPFCs_age_new %>% group_by(Subject, Epoch, age, visitno) %>% filter(Epoch == 'Delay'|Epoch == 'Fix') %>% summarise(meanPower_trad = mean(Power, na.rm=T), sdPower_trad = sd(Power, na.rm=T)) 
  
  AvgDLPFC$inverseAge <- 1/AvgDLPFC$age
  
  ### 8.3.1 Merge Traditional with Spectral Events ----
  
  avgRegions <- aggregate(.~Subject+Epoch,regionLevel_new[c(1,3:9, 11:14)], mean)
  
  spectralTrad <- merge(AvgDLPFC,avgRegions, by = c("Subject", "Epoch"))
 
 
  ## 8.4 Total Power vs Age ----
  
  # gam
   lunaize(ggplot(data = spectralTrad, 
                 aes(x = age.x, y = meanPower_trad, color = Epoch)) +
             geom_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 4, fx = T),alpha=.2,size=.75)) 
   

  gam.model <-  gamm(meanPower_trad  ~ s(age.x, k=3) + Epoch, data = spectralTrad,  random=list(Subject=~1))
  summary(gam.model$gam)
   
   
   # linear
   lunaize(ggplot(data = spectralTrad, 
                  aes(x = age.x, y = meanPower_trad, color = Epoch)) + stat_smooth(method="lm", formula=y~I(1/x), alpha = 0.2))+ ylab("Total Power") +xlab("Age")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
   
  
   lm.model <- lm(meanPower_trad ~ age.x + Epoch, data = spectralTrad)
   car::Anova(lm.model)
   summ(lm.model)
   
   
  
   ## 8.5 Mediation ----
   
   AvgDLPFCwide <- DLPFCs_age_new %>% group_by(Subject, Epoch, age, visitno) %>% filter(Epoch == 'Delay'|Epoch == 'Fix') %>% summarise(meanPower_trad = mean(Power, na.rm=T), sdPower_trad = sd(Power, na.rm=T)) %>% tidyr::pivot_wider(id_cols = c('Subject','age'), names_from='Epoch', values_from=c('meanPower_trad', 'sdPower_trad'))
   
   regionLevel_new_wide <- regionLevel_new %>% tidyr::pivot_wider(id_cols = c('Subject','Region','age'), names_from='Epoch', values_from=c('Gamma_Trial_Power', 'Gamma_Event_Number', 'Gamma_Event_Duration')) 
   
   avgRegions <- aggregate(.~Subject,regionLevel_new_wide[c(1,3:9)], mean)
   
   spectralTrad_med <- merge(AvgDLPFCwide,avgRegions, by = "Subject", all.x = T, all.y = T)

  
  ### 8.5.1 Event Number Fix and Total Power Fix ----
  
   #the effect of age on traditional power (c)
   model.0 <- lm(meanPower_trad_Fix ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.0))
   print(summary(model.0))
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   
   #the effect of age on spectral measure  (a)
   model.M <- lm(Gamma_Event_Number_Fix ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.M))
   print(summary(model.M))
   
   #the effect of spectral measure on trad power while controlling for age (b) 
   model.Y <- lm(meanPower_trad_Fix ~ age.x + Gamma_Event_Number_Fix, data = spectralTrad_med)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   
   
   results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = 'Gamma_Event_Number_Fix', boot = FALSE, sims = 1000)
   summary(results)
   

  ### 8.5.2 Event Number Delay and Total Power Delay ----
  
   #the effect of age on traditional power (c)
   model.0 <- lm(meanPower_trad_Delay ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.0))
   print(summary(model.0))
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   
   #the effect of age on spectral measure  (a)
   model.M <- lm(Gamma_Event_Number_Delay ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.M))
   print(summary(model.M))
   
   #the effect of spectral measure on trad power while controlling for age (b) 
   model.Y <- lm(meanPower_trad_Delay ~ age.x + Gamma_Event_Number_Delay, data = spectralTrad_med)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   
   
   results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = 'Gamma_Event_Number_Delay', boot = FALSE, sims = 1000)
   
   summary(results)
   
  
  ### 8.5.3 Trial Power Fix and Total Power Fix ----
  
   #the effect of age on traditional power (c)
   model.0 <- lm(meanPower_trad_Fix ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.0))
   print(summary(model.0))
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   
   #the effect of age on spectral measure  (a)
   model.M <- lm(Gamma_Trial_Power_Fix ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.M))
   print(summary(model.M))
   
   #the effect of spectral measure on trad power while controlling for age (b) 
   model.Y <- lm(meanPower_trad_Fix ~ age.x + Gamma_Trial_Power_Fix, data = spectralTrad_med)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   
   
   results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = 'Gamma_Trial_Power_Fix', boot = FALSE, sims = 1000)
   summary(results)
   
   
   ### 8.5.4 Trial Power Delay and Total Power Delay ----
  
   #the effect of age on traditional power (c)
   model.0 <- lm(meanPower_trad_Delay ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.0))
   print(summary(model.0))
   #anova(model.00, model.0) # significant improvement in fit when you include age 
   
   
   #the effect of age on spectral measure  (a)
   model.M <- lm(Gamma_Trial_Power_Delay ~ age.x, data = spectralTrad_med)
   print(car::Anova(model.M))
   print(summary(model.M))
   
   #the effect of spectral measure on trad power while controlling for age (b) 
   model.Y <- lm(meanPower_trad_Delay ~ age.x + Gamma_Trial_Power_Delay, data = spectralTrad_med)
   print(car::Anova(model.Y))
   print(summary(model.Y))
   
   
   results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = 'Gamma_Trial_Power_Delay', boot = FALSE, sims = 1000)
  
  summary(results)
 
  
  # 9.0 Loops ----
  
  ## 9.1 Spectral Event vs Age Loop ----
  
  eegVars <- c('Gamma_Trial_Power','Gamma_Trial_Power_sd', 'Gamma_Event_Number', 'Gamma_Event_Number_sd', 'Gamma_Event_Duration', 'Gamma_Event_Duration_sd')
  
  output <- c()  
  for (eegVar in eegVars) {
    
    model <- paste0(eegVar, ' ~ ', 's(age, k = 3) + Region + Epoch')
    model.out <- (mgcv::gamm(as.formula(model), random=list(Subject=~1), data = regionLevel_new))
    
    F <- summary(model.out$gam)$s.table[1,3]
    p <- summary(model.out$gam)$s.table[1,4]
    RegionF <- summary(model.out$gam)$pTerms.table[1,2]
    Regionp <- summary(model.out$gam)$pTerms.table[1,3]
    epochF <- summary(model.out$gam)$pTerms.table[2,2]
    epochp <- summary(model.out$gam)$pTerms.table[2,3]
    
    pcor <-  p.adjust((p), method = "bonferroni", n = 6)
    Regionpcor <- p.adjust((Regionp), method = "bonferroni", n = 6)
    epochpcor <- p.adjust((epochp), method = "bonferroni", n = 6)
    
    
    output <- rbind(output, data.frame(eegVar, F, p, pcor, RegionF, Regionp, Regionpcor, epochF, epochp, epochpcor))
  }
  output  
  
  ### 9.1.1 Spectral Event Age Interactions ----
  
  eegVars <- c('Gamma_Trial_Power','Gamma_Trial_Power_sd', 'Gamma_Event_Number', 'Gamma_Event_Number_sd', 'Gamma_Event_Duration', 'Gamma_Event_Duration_sd')
  regionLevel_new$Epoch <- as.factor(regionLevel_new$Epoch)
  regionLevel_new$Subject <- as.factor(regionLevel_new$Subject)
  
  #### 9.1.1.1 Interaction with Region ----
  regionLevel_new$oRegion <- ordered(regionLevel_new$Region, levels = c("L DLPFC","R DLPFC")) # LDLPFC will be the reference group

  
  output <- c()  
  for (eegVar in eegVars) {
    
    fullSubjs <- regionLevel_new %>% filter(!is.na(!!as.symbol(eegVar))) %>% group_by(Subject) %>% tally() %>% filter(n >= 3)
    
    model <- paste0(eegVar, ' ~ ', 'oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Epoch')
    model.out <- mgcv::gamm(as.formula(model), random=list(Subject=~1), data = regionLevel_new %>% filter(Subject %in% fullSubjs$Subject))
    
    F <- summary(model.out$gam)$s.table[2,3]
    p <- summary(model.out$gam)$s.table[2,4]
    
    pcor <-  p.adjust((p), method = "bonferroni", n = 6)
  
    output <- rbind(output, data.frame(eegVar, F, p, pcor))
  }
   output  
   
   
   #### 9.1.1.1 Interaction with Epoch ----
   regionLevel_new$oEpoch <- ordered(regionLevel_new$Epoch, levels = c("Fix","Delay")) # fix will be the reference group
   
   
   output <- c()  
   for (eegVar in eegVars) {
     
     fullSubjs <- regionLevel_new %>% filter(!is.na(!!as.symbol(eegVar))) %>% group_by(Subject) %>% tally() %>% filter(n >= 3)
     
     model <- paste0(eegVar, ' ~ ', 'oEpoch + s(age, k = 4, fx = T) + s(age, by = oEpoch, k = 4, fx = T) + Region')
     model.out <- mgcv::gamm(as.formula(model), random=list(Subject=~1), data = regionLevel_new %>% filter(Subject %in% fullSubjs$Subject))
     
     F <- summary(model.out$gam)$s.table[2,3]
     p <- summary(model.out$gam)$s.table[2,4]
     
     pcor <-  p.adjust((p), method = "bonferroni", n = 6)
     
     output <- rbind(output, data.frame(eegVar, F, p, pcor))
   }
   output  
  
  ## 9.2 Spectral Event vs Behavior loop: Delay Period ----

  ### 9.2.1 Average Main Effects ---- 
  
  eegVars <- c('Gamma_Trial_Power', 'Gamma_Event_Number', 'Gamma_Event_Duration')  
  behVars <- c('absBestError','Latency')
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '+ age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Delay')))  
      b <- model.out$coefficients[2,1]
      t <- model.out$coefficients[2,3]
      p <- model.out$coefficients[2,4]
      
      ageb <- model.out$coefficients[3,1]
      aget <- model.out$coefficients[3,3]
      agep <- model.out$coefficients[3,4]
      
      regionb <- model.out$coefficients[4,1]
      regiont <- model.out$coefficients[4,3]
      regionp <- model.out$coefficients[4,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
      agepcor <- p.adjust((agep), method = "bonferroni", n = 4)
      regionpcor <- p.adjust((regionp), method = "bonferroni", n = 4)
      
      
      
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor, ageb, aget, agep, agepcor, regionb, regionp, regiont, regionpcor))
      
    }
  }  
  output  
  
  ### 9.2.2 Average Age Interactions ----
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '*age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Delay')))  
      b <- model.out$coefficients[5,1]
      t <- model.out$coefficients[5,3]
      p <- model.out$coefficients[5,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
     
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor))
      
    }
  }  
  output  

 ### 9.2.3 Variability Main Effects ----

  eegVars <- c('Gamma_Trial_Power_sd', 'Gamma_Event_Number_sd', 'Gamma_Event_Duration_sd')  
  behVars <- c('absBestError_sd','LatencySD')
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '+ age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Delay')))  
      b <- model.out$coefficients[2,1]
      t <- model.out$coefficients[2,3]
      p <- model.out$coefficients[2,4]
      
      ageb <- model.out$coefficients[3,1]
      aget <- model.out$coefficients[3,3]
      agep <- model.out$coefficients[3,4]
      
      regionb <- model.out$coefficients[4,1]
      regiont <- model.out$coefficients[4,3]
      regionp <- model.out$coefficients[4,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
      agepcor <- p.adjust((agep), method = "bonferroni", n = 4)
      regionpcor <- p.adjust((regionp), method = "bonferroni", n = 4)
      
      
      
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor, ageb, aget, agep, agepcor, regionb, regionp, regiont, regionpcor))
      
    }
  }  
  output  

  ### 9.2.4 Variability Age Interactions ----
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '*age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Delay')))  
      b <- model.out$coefficients[5,1]
      t <- model.out$coefficients[5,3]
      p <- model.out$coefficients[5,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
      
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor))
      
    }
  }  
  output  


  ## 9.3 Spectral Event vs Behavior loop: Fixation Period----
  
  ### 9.3.1 Average Main Effects ---- 
  
  eegVars <- c('Gamma_Trial_Power', 'Gamma_Event_Number', 'Gamma_Event_Duration')  
  behVars <- c('absBestError','Latency')
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '+ age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Fix')))  
      b <- model.out$coefficients[2,1]
      t <- model.out$coefficients[2,3]
      p <- model.out$coefficients[2,4]
      
      ageb <- model.out$coefficients[3,1]
      aget <- model.out$coefficients[3,3]
      agep <- model.out$coefficients[3,4]
      
      regionb <- model.out$coefficients[4,1]
      regiont <- model.out$coefficients[4,3]
      regionp <- model.out$coefficients[4,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
      agepcor <- p.adjust((agep), method = "bonferroni", n = 4)
      regionpcor <- p.adjust((regionp), method = "bonferroni", n = 4)
      
      
      
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor, ageb, aget, agep, agepcor, regionb, regionp, regiont, regionpcor))
      
    }
  }  
  output  
  
  ### 9.3.2 Average Age Interactions ----
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '*age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Delay')))  
      b <- model.out$coefficients[5,1]
      t <- model.out$coefficients[5,3]
      p <- model.out$coefficients[5,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
      
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor))
      
    }
  }  
  output  
  
  ### 9.3.3 Variability Main Effects ----
  
  eegVars <- c('Gamma_Trial_Power_sd', 'Gamma_Event_Number_sd', 'Gamma_Event_Duration_sd')  
  behVars <- c('absBestError_sd','LatencySD')
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '+ age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Fix')))  
      b <- model.out$coefficients[2,1]
      t <- model.out$coefficients[2,3]
      p <- model.out$coefficients[2,4]
      
      ageb <- model.out$coefficients[3,1]
      aget <- model.out$coefficients[3,3]
      agep <- model.out$coefficients[3,4]
      
      regionb <- model.out$coefficients[4,1]
      regiont <- model.out$coefficients[4,3]
      regionp <- model.out$coefficients[4,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
      agepcor <- p.adjust((agep), method = "bonferroni", n = 4)
      regionpcor <- p.adjust((regionp), method = "bonferroni", n = 4)
      
      
      
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor, ageb, aget, agep, agepcor, regionb, regionp, regiont, regionpcor))
      
    }
  }  
  output  
  
  ### 9.3.4 Variability Age Interactions ----
  
  output <- c()  
  for (eegVar in eegVars) {
    for (behVar in behVars) {
      model <- as.formula(paste0(eegVar, ' ~ ', behVar, '*age + Region'))
      model.out <- summary(lm(model, data = gamma_long %>% filter(Epoch == 'Fix')))  
      b <- model.out$coefficients[5,1]
      t <- model.out$coefficients[5,3]
      p <- model.out$coefficients[5,4]
      
      pcor <- p.adjust((p), method = "bonferroni", n = 4)
      
      output <- rbind(output, data.frame(eegVar, behVar, b, t, p, pcor))
      
    }
  }  
  output  
  
  
  