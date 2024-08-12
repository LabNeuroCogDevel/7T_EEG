
# looking at behavior eeg and age interactions 
# only looking at the eeg measures that had significant age related changes cause who cares otherwise 
# only going to look at it at the subject level cause I just want everything on the same level 

behvaiorAgeInteractions_Method2 <- function () {

source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R")

alldf <- DelayOnly_Sublevel()
alldata_delayOnly <- alldf$alldata_delayOnly

avgBehavior <- Behavior_Sublevel_Maria()

allData <- merge(alldata_delayOnly, avgBehavior, by = c("Subject", "age"))

# behvaior by age 
lunaize(ggplot(data = allData, aes(x = age, y = mgsLatency)) + stat_smooth(method = "gam") + geom_point())
behaviorAge <- lm(mgsLatency ~ age, data = allData)
summary(behaviorAge)

# eeg by age 
lunaize(ggplot(data = allData, aes(x = age, y = Alpha.log1p_Trial_Power)) + stat_smooth(method = "gam") + geom_point())
behaviorAge <- lm(Theta.Event_Number ~ age, data = allData)
summary(behaviorAge)

# behavior eeg age interaction 
model <- lm(absPositionError ~ Gamma.log1p_Trial_Power_Variability*age, data = allData)
summary(model)


## individual regions interactions 
individualChannelDF <- DelayOnly_IndividualChannels_TrialLevel()
individualChannel_Regions <- IndividualChannels_GroupedintoRegions()
allFrontal <- individualChannel_Regions$allFrontal
avgFrontal <- aggregate(.~Subject, data = allFrontal, mean)
allFrontal_sd <- individualChannel_Regions$allFrontal_sd

frontal_behavior <- merge(avgFrontal, allFrontal_sd, by = "Subject", suffixes = c("", "_sd")) %>% merge(., avgBehavior, by = c("Subject", "age"))


# outlier detection 

outliers <- function(x) {
  
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
  
  
}

frontal_behavior_new <- frontal_behavior

cols = names(frontal_behavior[c(6:17,21:32)])
for (col in cols) {
  
  indx <- outliers(frontal_behavior[[col]])
  
  frontal_behavior_new[[col]] <- Map(replace, frontal_behavior_new[[col]], indx, NA)
  frontal_behavior_new[[col]] <- as.numeric(frontal_behavior_new[[col]])
  
}  


# eeg by age 
lunaize(ggplot(data = frontal_behavior_new, aes(x = age, y = Theta_Event_Duration)) + stat_smooth(method = "gam") + geom_point())
behaviorAge <- lm(Theta_Event_Duration ~ age, data = frontal_behavior_new)
summary(behaviorAge)

# behavior eeg age interaction 
model <- lm(absPositionError ~ Theta_Trial_Power*age, data = frontal_behavior_new)
summary(model)

#DLPFC
allDLPFC <- individualChannel_Regions$allDLPFC
avgDLPFC <- aggregate(.~Subject, data = allDLPFC, mean)
allDLPFC_sd <- individualChannel_Regions$allDLPFC_sd

DLPFC_behavior <- merge(avgDLPFC, allDLPFC_sd, by = "Subject", suffixes = c("", "_sd")) %>% merge(., avgBehavior, by = c("Subject", "age"))


# outlier detection 

outliers <- function(x) {
  
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
  
  
}

DLPFC_behavior_new <- DLPFC_behavior

cols = names(DLPFC_behavior[c(6:17,21:32)])
for (col in cols) {
  
  indx <- outliers(DLPFC_behavior[[col]])
  
  DLPFC_behavior_new[[col]] <- Map(replace, DLPFC_behavior_new[[col]], indx, NA)
  DLPFC_behavior_new[[col]] <- as.numeric(DLPFC_behavior_new[[col]])
  
}  


# eeg by age 
lunaize(ggplot(data = DLPFC_behavior_new, aes(x = age, y = Theta_Event_Duration)) + stat_smooth(method = "gam") + geom_point())
behaviorAge <- lm(Theta_Event_Duration ~ age, data = DLPFC_behavior_new)
summary(behaviorAge)

# behavior eeg age interaction 
model <- lm(absPositionError ~ Theta_Event_Duration*age, data = DLPFC_behavior_new)
summary(model)


# Parietal

allParietal <- individualChannel_Regions$allParietal
avgParietal <- aggregate(.~Subject, data = allParietal, mean)
allParietal_sd <- individualChannel_Regions$allParietal_sd

Parietal_behavior <- merge(avgParietal, allParietal_sd, by = "Subject", suffixes = c("", "_sd")) %>% merge(., avgBehavior, by = c("Subject", "age"))

Parietal_behavior_new <- Parietal_behavior

cols = names(Parietal_behavior[c(6:17,21:32)])
for (col in cols) {
  
  indx <- outliers(Parietal_behavior[[col]])
  
  Parietal_behavior_new[[col]] <- Map(replace, Parietal_behavior_new[[col]], indx, NA)
  Parietal_behavior_new[[col]] <- as.numeric(Parietal_behavior_new[[col]])
  
}  


# eeg by age 
lunaize(ggplot(data = Parietal_behavior_new, aes(x = age, y = log_Theta_Power_sd)) + stat_smooth(method = "gam") + geom_point())
behaviorAge <- lm(log_Theta_Power_sd ~ age, data = Parietal_behavior_new)
summary(behaviorAge)

# behavior eeg age interaction 
model <- gam(absPositionError ~ Theta_Trial_Power*age, data = Parietal_behavior_new)
summary(model)

Parietal_behavior_new$Group <- as.factor(Parietal_behavior_new$Group)

#graph by age groups
lunaize(ggplot(data = Parietal_behavior_new, aes(x = Alpha_Event_Duration, y = absPositionError, group = Group, color = Group)) + stat_smooth(method = "lm", alpha = 0.05, se = T) + ggtitle("") + xlab("Parietal Alpha Duration (s)") + ylab("Position Error (degs)")) +theme(plot.title = element_text(hjust = 0.5))



# growth model
gam.model <- gam(Alpha_Event_Duration ~ s(age), data = Parietal_behavior_new)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
powerPlot <- gam_growthrate_plot(Parietal_behavior_new, gam.model, gam.growthrate, agevar = 'age', yvar = 'Alpha_Event_Duration', draw_points = T)


# Occipital

allOccipital <- individualChannel_Regions$allOccipital
avgOccipital <- aggregate(.~Subject, data = allOccipital, mean)
allOccipital_sd <- individualChannel_Regions$allOccipital_sd

Occipital_behavior <- merge(avgOccipital, allOccipital_sd, by = "Subject", suffixes = c("", "_sd")) %>% merge(., avgBehavior, by = c("Subject", "age"))

Occipital_behavior_new <- Occipital_behavior

cols = names(Occipital_behavior[c(6:17,21:32)])
for (col in cols) {
  
  indx <- outliers(Occipital_behavior[[col]])
  
  Occipital_behavior_new[[col]] <- Map(replace, Occipital_behavior_new[[col]], indx, NA)
  Occipital_behavior_new[[col]] <- as.numeric(Occipital_behavior_new[[col]])
  
}  


# eeg by age 
lunaize(ggplot(data = Occipital_behavior_new, aes(x = age, y = log_Alpha_Power)) + stat_smooth(method = "gam") + geom_point())
behaviorAge <- lm(log_Alpha_Power ~ age, data = Occipital_behavior_new)
summary(behaviorAge)

# behavior eeg age interaction 
model <- lm(absPositionError ~ log_Alpha_Power*age, data = Occipital_behavior_new)
summary(model)

Occipital_behavior_new$Group <- as.factor(Occipital_behavior_new$Group)

#graph by age groups
lunaize(ggplot(data = Occipital_behavior_new, aes(x = log_Alpha_Power, y = absPositionError, group = Group, color = Group)) + stat_smooth(method = "lm", alpha = 0.05, se = T) + ggtitle("") + xlab("Occiptial Alpha Power (log)") + ylab("Position Error (degs)")) +theme(plot.title = element_text(hjust = 0.5))


# growth model
gam.model <- gam(log_Alpha_Power ~ s(age), data = Occipital_behavior_new)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
powerPlot <- gam_growthrate_plot(Occipital_behavior_new, gam.model, gam.growthrate, agevar = 'age', yvar = 'log_Alpha_Power', draw_points = T)



# put all regions into one dataframe
allRegions <- merge(frontal_behavior_new, Parietal_behavior_new[1:40], by = c("Subject", "age", "Group", "inverseAge.x"), suffixes = c("_frontal", "_parietal"))

colnames(Occipital_behavior_new)[c(6:21, 25:40)] <- paste(colnames(Occipital_behavior_new)[c(6:21, 25:40)], "occipital", sep = "_")

allRegions <- merge(allRegions, Occipital_behavior_new[1:40], by = c("Subject", "age", "Group", "inverseAge.x"))

save(allRegions, file = "ShaneAlldata.Rda")





# Whole Brain

allWholeBrain <- individualChannel_Regions$allWholeBrain
avgWholeBrain <- aggregate(.~Subject, data = allWholeBrain, mean)
allWholeBrain_sd <- individualChannel_Regions$allWholeBrain_sd

WholeBrain_behavior <- merge(avgWholeBrain, allWholeBrain_sd, by = "Subject", suffixes = c("", "_sd")) %>% merge(., avgBehavior, by = c("Subject", "age"))

WholeBrain_behavior_new <- WholeBrain_behavior

cols = names(WholeBrain_behavior[c(6:17,21:32)])
for (col in cols) {
  
  indx <- outliers(WholeBrain_behavior[[col]])
  
  WholeBrain_behavior_new[[col]] <- Map(replace, WholeBrain_behavior_new[[col]], indx, NA)
  WholeBrain_behavior_new[[col]] <- as.numeric(WholeBrain_behavior_new[[col]])
  
}  


# eeg by age 
lunaize(ggplot(data = WholeBrain_behavior_new, aes(x = age, y = log_Theta_Power_sd)) + stat_smooth(method = "gam") + geom_point())
behaviorAge <- lm(log_Theta_Power_sd ~ age, data = WholeBrain_behavior_new)
summary(behaviorAge)

# behavior eeg age interaction 
model <- lm(absPositionError ~ log_Theta_Power_sd*age, data = WholeBrain_behavior_new)
summary(model)

# growth model
gam.model <- gam(log_Alpha_Power ~ s(age), data = WholeBrain_behavior_new)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
powerPlot <- gam_growthrate_plot(WholeBrain_behavior_new, gam.model, gam.growthrate, agevar = 'age', yvar = 'log_Alpha_Power', draw_points = T)

}


behaviorAgeInteractions_Method1 <- function () {
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R")
  
  regiontoSpectraldata <- RegionstoSpectralEventData_SubjectLevel()
  avgBehavior <- Behavior_Sublevel_Maria()
  
  allData <- merge(regiontoSpectraldata, avgBehavior, by = c("Subject", "age"))
  
  
  # behvaior by age 
  lunaize(ggplot(data = allData, aes(x = age, y = absPositionError)) + stat_smooth(method = "gam") + geom_point())
  behaviorAge <- lm(absPositionError ~ age, data = allData)
  summary(behaviorAge)
  
  # eeg by age 
  lunaize(ggplot(data = allData, aes(x = age, y = Alpha_Event_Number_Occipital)) + stat_smooth(method = "gam") + geom_point())
  behaviorAge <- lm(log_Alpha_Power_Occipital ~ age, data = allData)
  summary(behaviorAge)
  
  # behavior eeg age interaction 
  model <- lm(absPositionError ~ log_Alpha_Power_Parietal*age, data = allData)
  summary(model)
  
  
  
}
