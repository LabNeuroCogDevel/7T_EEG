
sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)


# put all four bands on one graph 
alldf_subLevel <- DelayOnly_Sublevel()
delayOnly_SubLevel <- alldf_subLevel$alldata_delayOnly


alphavars <- grep("Alpha",names(delayOnly_SubLevel),value=TRUE)
allAlpha <- delayOnly_SubLevel[,c("Subject", "age",alphavars)]
colnames(allAlpha) <-  sub("Alpha.", "", colnames(allAlpha))
allAlpha$Band <- "Alpha"


thetavars <- grep("Theta",names(delayOnly_SubLevel),value=TRUE)
allTheta <- delayOnly_SubLevel[,c("Subject", "age",thetavars)]
colnames(allTheta) <-  sub("Theta.", "", colnames(allTheta))
allTheta$Band <- "Theta"

betavars <- grep("Beta",names(delayOnly_SubLevel),value=TRUE)
allBeta <- delayOnly_SubLevel[,c("Subject", "age",betavars)]
colnames(allBeta) <-  sub("Beta.", "", colnames(allBeta))
allBeta$Band <- "Beta"


gammavars <- grep("Gamma",names(delayOnly_SubLevel),value=TRUE)
allGamma <- delayOnly_SubLevel[,c("Subject", "age",gammavars)]
colnames(allGamma) <-  sub("Gamma.", "", colnames(allGamma))
allGamma$Band <- "Gamma"

alldata <- rbind(allAlpha,allTheta) %>% rbind(., allBeta) %>% rbind(.,allGamma)

lunaize(ggplot(data = alldata, aes(x = age, y = Event_Duration_Variability, group = Band, color = Band)) + stat_smooth(method = "gam") + ggtitle("Event Duration Variability") + xlab("Age") + ylab("Event Duration Var (s)")) +theme(plot.title = element_text(hjust = 0.5))



# alpha

regiontoSpectraldata <- RegionstoSpectralEventData_SubjectLevel()
alphavars <- grep("Alpha",names(regiontoSpectraldata),value=TRUE)
allAlpha_regions <- regiontoSpectraldata[,c("Subject", "age",alphavars)]

frontalvars <- grep("Frontal",names(allAlpha_regions),value=TRUE)
alphafrontal <- allAlpha_regions[,c("Subject", "age",frontalvars)]
alphafrontal$Band <- "Alpha"
alphafrontal$Region <- "Frontal"
colnames(alphafrontal) <-  sub("Alpha.", "", colnames(alphafrontal))
colnames(alphafrontal) <-  sub("_Frontal", "", colnames(alphafrontal))
colnames(alphafrontal) <-  sub("_sd", "_Variability", colnames(alphafrontal))
colnames(alphafrontal) <-  sub("log", "log1p_Trial", colnames(alphafrontal))


parietalvars <- grep("Parietal",names(allAlpha_regions),value=TRUE)
alphaparietal <- allAlpha_regions[,c("Subject", "age",parietalvars)]
alphaparietal$Band <- "Alpha"
alphaparietal$Region <- "Parietal"
colnames(alphaparietal) <-  sub("Alpha.", "", colnames(alphaparietal))
colnames(alphaparietal) <-  sub("_Parietal", "", colnames(alphaparietal))
colnames(alphaparietal) <-  sub("_sd", "_Variability", colnames(alphaparietal))
colnames(alphaparietal) <-  sub("log", "log1p_Trial", colnames(alphaparietal))


occipitalvars <- grep("Occipital",names(allAlpha_regions),value=TRUE)
alphaoccipital <- allAlpha_regions[,c("Subject", "age",occipitalvars)]
alphaoccipital$Band <- "Alpha"
alphaoccipital$Region <- "Occipital"
colnames(alphaoccipital) <-  sub("Alpha.", "", colnames(alphaoccipital))
colnames(alphaoccipital) <-  sub("_Occipital", "", colnames(alphaoccipital))
colnames(alphaoccipital) <-  sub("_sd", "_Variability", colnames(alphaoccipital))
colnames(alphaoccipital) <-  sub("log", "log1p_Trial", colnames(alphaoccipital))


allAlpha$Region <- "Whole Brain"

allAlphadata <- rbind(alphafrontal, alphaparietal) %>% rbind(.,alphaoccipital) %>% rbind(.,allAlpha)


lunaize(ggplot(data = allAlphadata, aes(x = age, y = Event_Duration_Variability, group = Region, color = Region)) + stat_smooth(method = "gam", alpha = 0.05, se = T) + ggtitle("Alpha Event Duration Variability") + xlab("Age") + ylab("Event Duration Var")) +theme(plot.title = element_text(hjust = 0.5))

behaviorAge <- lm(Event_Duration_Variability ~ age, data = allAlphadata[allAlphadata$Region == 'Frontal',])
summary(behaviorAge)

# merge with behavior 
avgBehavior <- Behavior_Sublevel_Maria()

allData <- merge(allAlphadata, avgBehavior, by = c("Subject", "age"))

# age interactions  
model <- lm(absPositionError ~ Event_Duration_Variability*age, data = allData[allData$Region == 'Whole Brain',])
summary(model)

# Gamma 

Gammavars <- grep("Gamma",names(regiontoSpectraldata),value=TRUE)
allGamma_regions <- regiontoSpectraldata[,c("Subject", "age",Gammavars)]

frontalvars <- grep("Frontal",names(allGamma_regions),value=TRUE)
Gammafrontal <- allGamma_regions[,c("Subject", "age",frontalvars)]
Gammafrontal$Band <- "Gamma"
Gammafrontal$Region <- "Frontal"
colnames(Gammafrontal) <-  sub("Gamma.", "", colnames(Gammafrontal))
colnames(Gammafrontal) <-  sub("_Frontal", "", colnames(Gammafrontal))
colnames(Gammafrontal) <-  sub("_sd", "_Variability", colnames(Gammafrontal))
colnames(Gammafrontal) <-  sub("log", "log1p_Trial", colnames(Gammafrontal))


parietalvars <- grep("Parietal",names(allGamma_regions),value=TRUE)
Gammaparietal <- allGamma_regions[,c("Subject", "age",parietalvars)]
Gammaparietal$Band <- "Gamma"
Gammaparietal$Region <- "Parietal"
colnames(Gammaparietal) <-  sub("Gamma.", "", colnames(Gammaparietal))
colnames(Gammaparietal) <-  sub("_Parietal", "", colnames(Gammaparietal))
colnames(Gammaparietal) <-  sub("_sd", "_Variability", colnames(Gammaparietal))
colnames(Gammaparietal) <-  sub("log", "log1p_Trial", colnames(Gammaparietal))


occipitalvars <- grep("Occipital",names(allGamma_regions),value=TRUE)
Gammaoccipital <- allGamma_regions[,c("Subject", "age",occipitalvars)]
Gammaoccipital$Band <- "Gamma"
Gammaoccipital$Region <- "Occipital"
colnames(Gammaoccipital) <-  sub("Gamma.", "", colnames(Gammaoccipital))
colnames(Gammaoccipital) <-  sub("_Occipital", "", colnames(Gammaoccipital))
colnames(Gammaoccipital) <-  sub("_sd", "_Variability", colnames(Gammaoccipital))
colnames(Gammaoccipital) <-  sub("log", "log1p_Trial", colnames(Gammaoccipital))


allGamma$Region <- "Whole Brain"

allGammadata <- rbind(Gammafrontal, Gammaparietal) %>% rbind(.,Gammaoccipital) %>% rbind(.,allGamma)


lunaize(ggplot(data = allGammadata, aes(x = age, y = Event_Duration_Variability, group = Region, color = Region)) + stat_smooth(method = "lm", alpha = 0.05, se = T) + ggtitle("Gamma Event Duration Variability") + xlab("Age") + ylab("Event Duration Var")) +theme(plot.title = element_text(hjust = 0.5))

behaviorAge <- lm(Event_Duration_Variability ~ age, data = allGammadata[allGammadata$Region == 'Whole Brain',])
summary(behaviorAge)


# merge with behavior 
avgBehavior <- Behavior_Sublevel_Maria()

allData <- merge(allGammadata, avgBehavior, by = c("Subject", "age"))

# age interactions  
model <- lm(absPositionError ~ Event_Duration_Variability*age, data = allData[allData$Region == 'Whole Brain',])
summary(model)

# plotting each electrode against age to see if there is noise (aka why dont we see the age related changes in gamma in method 2)
individualChannelDF <- DelayOnly_IndividualChannels_TrialLevel()
gammaDF <- individualChannelDF$GammaDelay_Age_Channel_new
Gamma_avgChannel_Events <- GammaDelay_Age_Channel %>% group_by(Subject, Channel, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
Gamma_avgChannel_Duration <- GammaDelay_Age_Channel %>% group_by(Subject, Channel, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
Gamma_avgChannel_Events$Channel <- as.factor(Gamma_avgChannel_Events$Channel)
Gamma_avgChannel_Duration$Channel <- as.factor(Gamma_avgChannel_Duration$Channel)


lunaize(ggplot(data = Gamma_avgChannel_Events, aes(x = age, y = Gamma_Event_Number, color = Channel)) + stat_smooth(method = "lm", alpha = 0.05, se = F) + ggtitle("Gamma Event Number") + xlab("Age") + ylab("Event Number")) +theme(plot.title = element_text(hjust = 0.5))

lunaize(ggplot(data = Gamma_avgChannel_Duration, aes(x = age, y = Gamma_Event_Duration, color = Channel)) + stat_smooth(method = "lm", alpha = 0.05, se = F) + ggtitle("Gamma Duration Number") + xlab("Age") + ylab("Event Duration")) +theme(plot.title = element_text(hjust = 0.5))



AlphaDF <- individualChannelDF$Alpha_Age_Channel_new
Alpha_avgChannel_Events <- AlphaDF %>% group_by(Subject, Channel, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
Alpha_avgChannel_Duration <- AlphaDF %>% group_by(Subject, Channel, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
Alpha_avgChannel_Events$Channel <- as.factor(Alpha_avgChannel_Events$Channel)
Alpha_avgChannel_Duration$Channel <- as.factor(Alpha_avgChannel_Duration$Channel)


lunaize(ggplot(data = Alpha_avgChannel_Events, aes(x = age, y = Alpha_Event_Number, color = Channel)) + stat_smooth(method = "lm", alpha = 0.05, se = F) + ggtitle("Alpha Event Number") + xlab("Age") + ylab("Event Number")) +theme(plot.title = element_text(hjust = 0.5))

lunaize(ggplot(data = Alpha_avgChannel_Duration, aes(x = age, y = Alpha_Event_Duration, color = Channel)) + stat_smooth(method = "lm", alpha = 0.05, se = F) + ggtitle("Alpha Duration Number") + xlab("Age") + ylab("Event Duration")) +theme(plot.title = element_text(hjust = 0.5))
