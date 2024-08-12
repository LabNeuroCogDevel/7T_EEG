
# Define Functions ----
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
  }
  
# Delay Period ----
## Load in Data Frames ----
### Previous spectral event subjects ----

gammaAllChannels_old <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Spectral_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
gammaAllChannels_old$Epoch <- 'Delay'


### New spectral event subjects ----

gammaAllChannels_new <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Spectral_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_additionalSubjects_20230530.csv') 
gammaAllChannels_new$Epoch <- 'Delay'

### Combine new and old subjects ----

gammaAllChannels_allSubs <- rbind(gammaAllChannels_old,gammaAllChannels_new)

### Channel locations ----
chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/ChannelLocs.csv')
chanLocs <- chanLocs[c(-30),] #remove POz
chanLocs$urchan <- (1:63)

# Merge spectral events with Channel locations
gammaAllChannels_allSubs$urchan <- gammaAllChannels_allSubs$Channel
gammaAllChannels_allSubs$logGammaPower <- log1p(gammaAllChannels_allSubs$Gamma_Trial_Power)

gammaAllChannels_allSubs_locs <- merge(gammaAllChannels_allSubs, chanLocs, by = "urchan")

gammaAllChannels_allSubs_locs_new <- gammaAllChannels_allSubs_locs

# outlier detection to rid bad trials
outlier_cols <- names(gammaAllChannels_allSubs_locs[c(5:7,9)])
na_outliers <- function(x) ifelse(!outliers(x), x, NA)
gammaAllChannels_allSubs_locs_new <- gammaAllChannels_allSubs_locs %>% group_by(Subject, Channel, Epoch) %>% mutate(across({outlier_cols}, na_outliers)) %>% ungroup 


### Average Trials ----
agg_fml <- cbind(Gamma_Trial_Power, logGammaPower, Gamma_Event_Number, Gamma_Event_Duration) ~ Subject + Epoch + urchan
dlpfc_96 <- gammaAllChannels_allSubs_locs_new %>% filter(Trial < 97)
aggregateTrials    <- aggregate(dlpfc_96, agg_fml, function(x) mean(x,na.rm=T), na.action=NULL)
aggregateTrials_sd <- aggregate(dlpfc_96, agg_fml, function(x) sd(x,na.rm=T)  , na.action=NULL)

SubLevelDelay <- merge(aggregateTrials, aggregateTrials_sd, by = c("Subject", "Epoch", "urchan"), suffixes = c("", "_sd")) 

write.csv(SubLevelDelay, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Spectral_Analysis/Spectral_events_analysis/GammaSpectralEvents_DelayOnly_channelLevel_merge7tEEG_20231110.csv')


# Fixation Period ----
## Load in Data Frames ----
### Previous spectral event subjects ----

gammaAllChannels_old <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Spectral_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_FIX.csv') 
gammaAllChannels_old$Epoch <- 'Fix'

### New spectral event subjects (NEED TO DO THIS) ----

gammaAllChannels_new <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Spectral_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_additionalSubjects_20230530.csv') 
gammaAllChannels_new$Epoch <- 'Fix'

### Combine new and old subjects ----
gammaAllChannels_allSubs <- gammaAllChannels_old

#gammaAllChannels_allSubs <- rbind(gammaAllChannels_old,gammaAllChannels_new) ADD BACK IN ONCE YOU RUN SPECTRAL EVENTS ON NEW SUBS FOR FIX 

### Channel locations ----
chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/ChannelLocs.csv')
chanLocs <- chanLocs[c(-30),] #remove POz
chanLocs$urchan <- (1:63)

# Merge spectral events with Channel locations
gammaAllChannels_allSubs$urchan <- gammaAllChannels_allSubs$Channel
gammaAllChannels_allSubs$logGammaPower <- log1p(gammaAllChannels_allSubs$Gamma_Trial_Power)

gammaAllChannels_allSubs_locs <- merge(gammaAllChannels_allSubs, chanLocs, by = "urchan")

gammaAllChannels_allSubs_locs_new <- gammaAllChannels_allSubs_locs

# outlier detection to rid bad trials
outlier_cols <- names(gammaAllChannels_allSubs_locs[c(5:7,9)])
na_outliers <- function(x) ifelse(!outliers(x), x, NA)
gammaAllChannels_allSubs_locs_new <- gammaAllChannels_allSubs_locs %>% group_by(Subject, Channel, Epoch) %>% mutate(across({outlier_cols}, na_outliers)) %>% ungroup 


### Average Trials ----
agg_fml <- cbind(Gamma_Trial_Power, logGammaPower, Gamma_Event_Number, Gamma_Event_Duration) ~ Subject + Epoch + urchan
dlpfc_96 <- gammaAllChannels_allSubs_locs_new %>% filter(Trial < 97)
aggregateTrials    <- aggregate(dlpfc_96, agg_fml, function(x) mean(x,na.rm=T), na.action=NULL)
aggregateTrials_sd <- aggregate(dlpfc_96, agg_fml, function(x) sd(x,na.rm=T)  , na.action=NULL)

SubLevelFix <- merge(aggregateTrials, aggregateTrials_sd, by = c("Subject", "Epoch", "urchan"), suffixes = c("", "_sd")) 

write.csv(SubLevel, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Spectral_Analysis/Spectral_events_analysis/GammaSpectralEvents_FixOnly_channelLevel_merge7tEEG_20231110.csv')


# Combine Delay and Fix ----
allSubsallEpochs <- rbind(SubLevelFix, SubLevelDelay)
write.csv(allSubsallEpochs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Spectral_Analysis/Spectral_events_analysis/GammaSpectralEvents_DelayandFix_channelLevel_merge7tEEG_20231110.csv')


