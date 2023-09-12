## 
## remove outliders from power timeseries, mean and sd, remove outliers from summary
##
## 202306xxSM - init
## 20230623WF - see spectral_events_wide.R for what goes into merge7T
##              watch out for aggregate's default behavior: complete.case (no row with any NA col)
##
# Define Functions ----
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  
# Load in Data Frames ----
## Previous spectral event subjects ----

gammaAllChannels_old <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
gammaAllChannels_old$Epoch <- 'Delay'


## New spectral event subjects ----

gammaAllChannels_new <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_additionalSubjects_20230530.csv') 
gammaAllChannels_new$Epoch <- 'Delay'

## Combine new and old subjects ----

gammaAllChannels_allSubs <- rbind(gammaAllChannels_old,gammaAllChannels_new)

## Channel locations ----
chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
chanLocs <- chanLocs[c(-30),] #remove POz
chanLocs$urchan <- (1:63)

# Merge spectral events with Channel locations
gammaAllChannels_allSubs$urchan <- gammaAllChannels_allSubs$Channel
gammaAllChannels_allSubs$logGammaPower <- log1p(gammaAllChannels_allSubs$Gamma_Trial_Power)

gammaAllChannels_allSubs_locs <- merge(gammaAllChannels_allSubs, chanLocs, by = "urchan")


# Define DLPFC Channels ----

RDLPFC <- filter(gammaAllChannels_allSubs_locs, gammaAllChannels_allSubs_locs$labels == 'F4' | gammaAllChannels_allSubs_locs$labels == 'F6'| gammaAllChannels_allSubs_locs$labels == 'F8')
RDLPFC$Region <- 'R DLPFC'

LDLPFC <- filter(gammaAllChannels_allSubs_locs, gammaAllChannels_allSubs_locs$labels == 'F3' | gammaAllChannels_allSubs_locs$labels == 'F5'| gammaAllChannels_allSubs_locs$labels == 'F7')
LDLPFC$Region <- 'L DLPFC'

bothDLPFCs <- rbind(RDLPFC, LDLPFC)


# Outlier Detection ----

bothDLPFCs_new <- bothDLPFCs %>% filter(Trial < 97) %>% group_by(Subject, Channel, Epoch) %>% mutate(Gamma_Trial_Power = ifelse(!outliers(Gamma_Trial_Power), Gamma_Trial_Power, NA)) %>% ungroup

bothDLPFCs_new <- bothDLPFCs_new %>% filter(Trial < 97) %>% group_by(Subject, Channel, Epoch) %>% mutate(Gamma_Event_Number = ifelse(!outliers(Gamma_Event_Number), Gamma_Event_Number, NA)) %>% ungroup

bothDLPFCs_new <- bothDLPFCs_new %>% filter(Trial < 97) %>% group_by(Subject, Channel, Epoch) %>% mutate(Gamma_Event_Duration = ifelse(!outliers(Gamma_Event_Duration), Gamma_Event_Duration, NA)) %>% ungroup

bothDLPFCs_new <- bothDLPFCs_new %>% filter(Trial < 97) %>% group_by(Subject, Channel, Epoch) %>% mutate(logGammaPower = ifelse(!outliers(logGammaPower), logGammaPower, NA)) %>% ungroup

# Average Trials ----
# 20230623WF - na.action=NULL so each column gets all the data available (otherwise complete.case removes any row with any NA)
agg_fml <- cbind(Gamma_Trial_Power, logGammaPower, Gamma_Event_Number, Gamma_Event_Duration) ~ Subject + Region + Epoch
dlpfc_96 <- bothDLPFCs_new %>% filter(Trial < 97)
aggregateTrials    <- aggregate(dlpfc_96, agg_fml, function(x) mean(x,na.rm=T), na.action=NULL)
aggregateTrials_sd <- aggregate(dlpfc_96, agg_fml, function(x) sd(x,na.rm=T)  , na.action=NULL)

regionLevel <- merge(aggregateTrials, aggregateTrials_sd, by = c("Subject", "Region", "Epoch"), suffixes = c("", "_sd")) 

# Subject level outlier detection ----

regionLevel_new <- regionLevel

outlier_cols <- names(regionLevel[c(4:11)])
na_outliers <- function(x) ifelse(!outliers(x), x, NA)
regionLevel_new <- regionLevel %>% group_by(Region, Epoch) %>% mutate(across({outlier_cols}, na_outliers)) %>% ungroup 
regionLevel_new$Region <- as.factor(regionLevel_new$Region)

write.csv(regionLevel_new, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_DLPFCs_spectralEvents_20230622.csv')
