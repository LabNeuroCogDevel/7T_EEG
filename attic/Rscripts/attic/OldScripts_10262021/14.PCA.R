
# RUN PCA ON THE DELAY AND REST DATA 

sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)

#all data, delay and fix, at the subject level 
alldf <- Only_Take_One_Delay_Bin_and_Fix()
alldata <- alldf$alldata

#only delay variables, trial level
alldf_trialLevel <- Only_Take_One_Delay_Bin_TrialLevel()
alldata_trialLevel <- alldf_trialLevel$alldata_TrialLevel
alldata_trialLevel$Task <- "Delay"
names(alldata_trialLevel) <- sub("Gamma_","",names(alldata_trialLevel))
names(alldata_trialLevel) <- sub("Beta_","",names(alldata_trialLevel))
names(alldata_trialLevel) <- sub("Theta_","",names(alldata_trialLevel))
names(alldata_trialLevel) <- sub("Alpha_","",names(alldata_trialLevel))

#only the delay variables at the subject level 
delayOnlyVars <- grep("Delay|Subject|age|idvalues|Group|visitno|inverse", names(alldata), value = TRUE)
delayOnlyVars <- delayOnlyVars[!grepl("Trial.x|Trial.y|Peak", delayOnlyVars)]
delayOnly <- alldata[,delayOnlyVars]
delayOnly <- subset(delayOnly, select = -c(Gamma.Trial_Power_Variability_Delay, Beta.Trial_Power_Variability_Delay, Alpha.Trial_Power_Variability_Delay,Theta.Trial_Power_Variability_Delay))


delayOnly$Task <- "Delay"
names(delayOnly) <- sub("_Delay","",names(delayOnly))

#resting state at the trial level 
alldf_RS <- Resting_State_Data_TrialLevel()
alldata_RS <- alldf_RS$alldata_RS_TrialLevel
alldata_RS$Task <- "Rest"

#resting state at the subject level
alldf_SubjectLevel_RS <- Resting_State_Data_SubjectLevel()
alldata_RS_SubLevel <- alldf_SubjectLevel_RS$rest_SubLevel
alldata_RS_SubLevel$Task <- "Rest"
alldata_RS_SubLevel <- subset(alldata_RS_SubLevel, select = -Trial)

#merging delay and rest, at the subject level 

delayRest <- merge(delayOnly, alldata_RS_SubLevel[1:30], by = "Subject", suffixes = c("_Delay", "_Rest"))
delayRest_long <- rbind(delayOnly, alldata_RS_SubLevel)

numericDF <- delayRest[-c(1,30:35,44:49,64)]
pcaDF <- prcomp(numericDF)
plot(pcaDF)
summary(pcaDF)
