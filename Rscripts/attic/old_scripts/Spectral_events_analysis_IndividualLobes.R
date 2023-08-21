if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2,lme4, multcomp, tidyr, tidyverse,wesanderson, reshape2, corrplot, ggExtra, cowplot, ggsci, mgcv, ggpointdensity, devtools, ggpubr, itsadug, useful, rmarkdown, tinytex, gratia, gamm4)
devtools::install_github('LabneuroCogDevel/LNCDR')
library(LNCDR)


#gamma Frontal lobe
FLgammaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/FrontalLobe_Gamma_Spectral_Analysis_Table_all.csv')
FLgammaPower$GammalogPower <- log(FLgammaPower$FrontalLobe_MeanAvgPower)
FLgammaPower$GammaAge <- (FLgammaPower$age)
FLgammaPower$GammaAvgEventNumber <- FLgammaPower$FrontalLobe_MeanAvgEventNumber
FLgammaPower$GammaAvgEventDuration <- FLgammaPower$FrontalLobe_MeanAvgEventDuration
FLgammaPower$GammaPeakPower <- log(FLgammaPower$FrontalLobe_MeanpeakPower)
FLgammaPower$AgeInverse <- 1/(FLgammaPower$age)


ggplot(data=FLgammaPower[FLgammaPower$visitno < 2,], aes(x=GammaAge, y=GammalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLgammaPower[FLgammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLgammaPower[FLgammaPower$visitno < 2 & FLgammaPower$GammaAvgEventDuration < 0.029, ], aes(x=GammaAge, y=GammaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLgammaPower[FLgammaPower$visitno < 2 ,], aes(x=GammaAge, y=GammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=FLgammaPower[FLgammaPower$visitno < 2 ,], GammalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=FLgammaPower[FLgammaPower$visitno < 2 ,], GammaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=FLgammaPower[FLgammaPower$visitno < 2 & FLgammaPower$GammaAvgEventDuration < 0.029,], GammaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=FLgammaPower[FLgammaPower$visitno < 2 ,], GammaPeakPower ~ AgeInverse )
anova(lm.model)


#gamma parietal
PLgammaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ParietalLobe_Gamma_Spectral_Analysis_Table_all.csv')
PLgammaPower$GammalogPower <- log(PLgammaPower$ParietalLobe_MeanAvgPower)
PLgammaPower$GammaAge <- (PLgammaPower$age)
PLgammaPower$GammaAvgEventNumber <- PLgammaPower$ParietalLobe_MeanAvgEventNumber
PLgammaPower$GammaAvgEventDuration <- PLgammaPower$ParietalLobe_MeanAvgEventDuration
PLgammaPower$GammaPeakPower <- log(PLgammaPower$ParietalLobe_MeanpeakPower)
PLgammaPower$AgeInverse <- 1/(PLgammaPower$age)


ggplot(data=PLgammaPower[PLgammaPower$visitno < 2,], aes(x=GammaAge, y=GammalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLgammaPower[PLgammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLgammaPower[PLgammaPower$visitno < 2 & PLgammaPower$GammaAvgEventDuration < 0.029, ], aes(x=GammaAge, y=GammaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLgammaPower[PLgammaPower$visitno < 2 ,], aes(x=GammaAge, y=GammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=PLgammaPower[PLgammaPower$visitno < 2 ,], GammalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=PLgammaPower[PLgammaPower$visitno < 2 ,], GammaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=PLgammaPower[PLgammaPower$visitno < 2 & FLgammaPower$GammaAvgEventDuration < 0.029,], GammaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=PLgammaPower[PLgammaPower$visitno < 2 ,], GammaPeakPower ~ AgeInverse )
anova(lm.model)


#gamma occipital
OLgammaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/OccipitalLobe_Gamma_Spectral_Analysis_Table_all.csv')
OLgammaPower$GammalogPower <- log(OLgammaPower$OccipitalLobe_MeanAvgPower)
OLgammaPower$GammaAge <- (OLgammaPower$age)
OLgammaPower$GammaAvgEventNumber <- OLgammaPower$OccipitalLobe_MeanAvgEventNumber
OLgammaPower$GammaAvgEventDuration <- OLgammaPower$OccipitalLobe_MeanAvgEventDuration
OLgammaPower$GammaPeakPower <- log(OLgammaPower$OccipitalLobe_MeanpeakPower)
OLgammaPower$AgeInverse <- 1/(OLgammaPower$age)


ggplot(data=OLgammaPower[OLgammaPower$visitno < 2,], aes(x=GammaAge, y=GammalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLgammaPower[OLgammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLgammaPower[OLgammaPower$visitno < 2 & OLgammaPower$GammaAvgEventDuration < 0.029, ], aes(x=GammaAge, y=GammaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLgammaPower[OLgammaPower$visitno < 2 ,], aes(x=GammaAge, y=GammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=OLgammaPower[OLgammaPower$visitno < 2 ,], GammalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=OLgammaPower[OLgammaPower$visitno < 2 ,], GammaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=OLgammaPower[OLgammaPower$visitno < 2 & OLgammaPower$GammaAvgEventDuration < 0.029,], GammaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=OLgammaPower[OLgammaPower$visitno < 2 ,], GammaPeakPower ~ AgeInverse )
anova(lm.model)


#gamma central
CgammaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Central_Gamma_Spectral_Analysis_Table_all.csv')
CgammaPower$GammalogPower <- log(CgammaPower$Central_MeanAvgPower)
CgammaPower$GammaAge <- (CgammaPower$age)
CgammaPower$GammaAvgEventNumber <- CgammaPower$Central_MeanAvgEventNumber
CgammaPower$GammaAvgEventDuration <- CgammaPower$Central_MeanAvgEventDuration
CgammaPower$GammaPeakPower <- log(CgammaPower$Central_MeanpeakPower)
CgammaPower$AgeInverse <- 1/(CgammaPower$age)


ggplot(data=CgammaPower[CgammaPower$visitno < 2,], aes(x=GammaAge, y=GammalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CgammaPower[CgammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CgammaPower[CgammaPower$visitno < 2 & CgammaPower$GammaAvgEventDuration < 0.029, ], aes(x=GammaAge, y=GammaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CgammaPower[CgammaPower$visitno < 2 ,], aes(x=GammaAge, y=GammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=CgammaPower[CgammaPower$visitno < 2 ,], GammalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=CgammaPower[CgammaPower$visitno < 2 ,], GammaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=CgammaPower[CgammaPower$visitno < 2 & CgammaPower$GammaAvgEventDuration < 0.029,], GammaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=CgammaPower[CgammaPower$visitno < 2 ,], GammaPeakPower ~ AgeInverse )
anova(lm.model)






#Theta Frontal lobe
FLThetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/FrontalLobe_Theta_Spectral_Analysis_Table_all.csv')
FLThetaPower$ThetalogPower <- log(FLThetaPower$FrontalLobe_MeanAvgPower)
FLThetaPower$ThetaAge <- (FLThetaPower$age)
FLThetaPower$ThetaAvgEventNumber <- FLThetaPower$FrontalLobe_MeanAvgEventNumber
FLThetaPower$thetaAvgEventDuration <- FLThetaPower$FrontalLobe_MeanAvgEventDuration
FLThetaPower$ThetaPeakPower <- log(FLThetaPower$FrontalLobe_MeanpeakPower)
FLThetaPower$AgeInverse <- 1/(FLThetaPower$age)


ggplot(data=FLThetaPower[FLThetaPower$visitno < 2,], aes(x=ThetaAge, y=ThetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLThetaPower[FLThetaPower$visitno < 2 & FLThetaPower$ThetaAvgEventNumber <.8,], aes(x=ThetaAge, y=ThetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLThetaPower[FLThetaPower$visitno < 2 , ], aes(x=ThetaAge, y=thetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLThetaPower[FLThetaPower$visitno < 2 ,], aes(x=ThetaAge, y=ThetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=FLThetaPower[FLThetaPower$visitno < 2 ,], ThetalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=FLThetaPower[FLThetaPower$visitno < 2 & FLThetaPower$ThetaAvgEventNumber <.8,], ThetaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=FLThetaPower[FLThetaPower$visitno < 2 ,], thetaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=FLThetaPower[FLThetaPower$visitno < 2 ,], ThetaPeakPower ~ AgeInverse )
anova(lm.model)



#Theta Parietal lobe
PLThetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ParietalLobe_Theta_Spectral_Analysis_Table_all.csv')
PLThetaPower$ThetalogPower <- log(PLThetaPower$ParietalLobe_MeanAvgPower)
PLThetaPower$ThetaAge <- (PLThetaPower$age)
PLThetaPower$ThetaAvgEventNumber <- PLThetaPower$ParietalLobe_MeanAvgEventNumber
PLThetaPower$thetaAvgEventDuration <- PLThetaPower$ParietalLobe_MeanAvgEventDuration
PLThetaPower$ThetaPeakPower <- log(PLThetaPower$ParietalLobe_MeanpeakPower)
PLThetaPower$AgeInverse <- 1/(PLThetaPower$age)


ggplot(data=PLThetaPower[PLThetaPower$visitno <2 & PLThetaPower$ThetalogPower < 2.5,], aes(x=ThetaAge, y=ThetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLThetaPower[PLThetaPower$visitno < 2  & PLThetaPower$ThetaAvgEventNumber <.8,], aes(x=ThetaAge, y=ThetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLThetaPower[PLThetaPower$visitno < 2 , ], aes(x=ThetaAge, y=thetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLThetaPower[PLThetaPower$visitno < 2 & PLThetaPower$ThetaPeakPower <5,], aes(x=ThetaAge, y=ThetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for log power
lm.model <- lm(data=PLThetaPower[PLThetaPower$visitno < 2 & PLThetaPower$ThetalogPower < 2.5,], ThetalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=PLThetaPower[PLThetaPower$visitno < 2 & PLThetaPower$ThetaAvgEventNumber <.8,], ThetaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=PLThetaPower[PLThetaPower$visitno < 2 ,], thetaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=PLThetaPower[PLThetaPower$visitno < 2 & PLThetaPower$ThetaPeakPower <5,], ThetaPeakPower ~ AgeInverse )
anova(lm.model)



#Theta Occipital lobe
OLThetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/OccipitalLobe_Theta_Spectral_Analysis_Table_all.csv')
OLThetaPower$ThetalogPower <- log(OLThetaPower$OccipitalLobe_MeanAvgPower)
OLThetaPower$ThetaAge <- (OLThetaPower$age)
OLThetaPower$ThetaAvgEventNumber <- OLThetaPower$OccipitalLobe_MeanAvgEventNumber
OLThetaPower$thetaAvgEventDuration <- OLThetaPower$OccipitalLobe_MeanAvgEventDuration
OLThetaPower$ThetaPeakPower <- log(OLThetaPower$OccipitalLobe_MeanpeakPower)
OLThetaPower$AgeInverse <- 1/(OLThetaPower$age)


ggplot(data=OLThetaPower[OLThetaPower$visitno < 2 & OLThetaPower$ThetalogPower <2,], aes(x=ThetaAge, y=ThetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLThetaPower[OLThetaPower$visitno < 2 ,], aes(x=ThetaAge, y=ThetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLThetaPower[OLThetaPower$visitno < 2 ,], aes(x=ThetaAge, y=thetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLThetaPower[OLThetaPower$visitno < 2 & OLThetaPower$ThetaPeakPower <5,], aes(x=ThetaAge, y=ThetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=OLThetaPower[OLThetaPower$visitno < 2 & OLThetaPower$ThetalogPower <2,], ThetalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=OLThetaPower[OLThetaPower$visitno < 2 ,], ThetaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=OLThetaPower[OLThetaPower$visitno < 2 ,], thetaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=OLThetaPower[OLThetaPower$visitno < 2 & OLThetaPower$ThetaPeakPower <5,], ThetaPeakPower ~ AgeInverse )
anova(lm.model)


#Theta Central
CThetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Central_Theta_Spectral_Analysis_Table_all.csv')
CThetaPower$ThetalogPower <- log(CThetaPower$Central_MeanAvgPower)
CThetaPower$ThetaAge <- (CThetaPower$age)
CThetaPower$ThetaAvgEventNumber <- CThetaPower$Central_MeanAvgEventNumber
CThetaPower$thetaAvgEventDuration <- CThetaPower$Central_MeanAvgEventDuration
CThetaPower$ThetaPeakPower <- log(CThetaPower$Central_MeanpeakPower)
CThetaPower$AgeInverse <- 1/(CThetaPower$age)


ggplot(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$ThetalogPower > -1,], aes(x=ThetaAge, y=ThetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$ThetaAvgEventNumber > 0.1,], aes(x=ThetaAge, y=ThetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$thetaAvgEventDuration >0.05,], aes(x=ThetaAge, y=thetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$ThetaPeakPower <5 & CThetaPower$ThetaPeakPower >1,], aes(x=ThetaAge, y=ThetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$ThetalogPower > -1,], ThetalogPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$ThetaAvgEventNumber > 0.1,], ThetaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$thetaAvgEventDuration >0.05 ,], thetaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=CThetaPower[CThetaPower$visitno < 2 & CThetaPower$ThetaPeakPower <5 & CThetaPower$ThetaPeakPower >1,], ThetaPeakPower ~ AgeInverse )
anova(lm.model)







#theta gamma ratio frontal 
FLthetagammaratio <- merge(FLgammaPower[FLgammaPower$visitno < 2,], FLThetaPower , by='idvalues')
FLthetagammaratio$thetagammalogpower<- FLthetagammaratio$ThetalogPower / FLthetagammaratio$GammalogPower
FLthetagammaratio$thetagammaPeakPower <- FLthetagammaratio$ThetaPeakPower/ FLthetagammaratio$GammaPeakPower
FLthetagammaratio$thetagammaDuration <- FLthetagammaratio$thetaAvgEventDuration/ FLthetagammaratio$GammaAvgEventDuration


ggplot(data=FLthetagammaratio[FLthetagammaratio$thetagammalogpower> 0,], aes(x=age.x, y=thetagammalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLthetagammaratio, aes(x=age.x, y=thetagammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLthetagammaratio[FLthetagammaratio$thetagammaDuration> 2,], aes(x=age.x, y=thetagammaDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=FLthetagammaratio[FLthetagammaratio$thetagammalogpower> 0,], thetagammalogpower ~ AgeInverse.x )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=FLthetagammaratio, thetagammaPeakPower ~ AgeInverse.x )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=FLthetagammaratio[FLthetagammaratio$thetagammaDuration> 2,], thetagammaDuration ~ AgeInverse.x )
anova(lm.model)

#theta gamma ratio parietal
PLthetagammaratio <- merge(PLgammaPower[PLgammaPower$visitno < 2,], PLThetaPower , by='idvalues')
PLthetagammaratio$thetagammalogpower<- PLthetagammaratio$ThetalogPower / PLthetagammaratio$GammalogPower
PLthetagammaratio$thetagammaPeakPower <- PLthetagammaratio$ThetaPeakPower/ PLthetagammaratio$GammaPeakPower
PLthetagammaratio$thetagammaduration <- PLthetagammaratio$thetaAvgEventDuration/ PLthetagammaratio$GammaAvgEventDuration


ggplot(data=PLthetagammaratio[PLthetagammaratio$thetagammalogpower < .7,], aes(x=age.x, y=thetagammalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLthetagammaratio[PLthetagammaratio$thetagammaPeakPower <0.9,], aes(x=age.x, y=thetagammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLthetagammaratio, aes(x=age.x, y=thetagammaduration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=PLthetagammaratio, thetagammalogpower ~ AgeInverse.x )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=PLthetagammaratio[PLthetagammaratio$thetagammaPeakPower <0.9,], thetagammaPeakPower ~ AgeInverse.x )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=PLthetagammaratio, thetagammaduration ~ AgeInverse.x )
anova(lm.model)



#theta gamma ratio occipital
OLthetagammaratio <- merge(PLgammaPower[OLgammaPower$visitno < 2,], OLThetaPower , by='idvalues')
OLthetagammaratio$thetagammalogpower<- OLthetagammaratio$ThetalogPower / OLthetagammaratio$GammalogPower
OLthetagammaratio$thetagammaPeakPower <- OLthetagammaratio$ThetaPeakPower/ OLthetagammaratio$GammaPeakPower
OLthetagammaratio$thetagammaduration <- OLthetagammaratio$thetaAvgEventDuration/ OLthetagammaratio$GammaAvgEventDuration


ggplot(data=OLthetagammaratio[OLthetagammaratio$thetagammalogpower < .5,], aes(x=age.x, y=thetagammalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLthetagammaratio, aes(x=age.x, y=thetagammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLthetagammaratio, aes(x=age.x, y=thetagammaduration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for log power
lm.model <- lm(data=OLthetagammaratio[OLthetagammaratio$thetagammalogpower < .5,], thetagammalogpower ~ AgeInverse.x )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=OLthetagammaratio, thetagammaPeakPower ~ AgeInverse.x )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=OLthetagammaratio, thetagammaduration ~ AgeInverse.x )
anova(lm.model)


#theta gamma ratio Central
Cthetagammaratio <- merge(PLgammaPower[CgammaPower$visitno < 2,], CThetaPower , by='idvalues')
Cthetagammaratio$thetagammalogpower<- Cthetagammaratio$ThetalogPower / Cthetagammaratio$GammalogPower
Cthetagammaratio$thetagammaPeakPower <- Cthetagammaratio$ThetaPeakPower/ Cthetagammaratio$GammaPeakPower
Cthetagammaratio$thetagammaduration <- Cthetagammaratio$thetaAvgEventDuration/ Cthetagammaratio$GammaAvgEventDuration


ggplot(data=Cthetagammaratio[Cthetagammaratio$thetagammalogpower > 0,], aes(x=age.x, y=thetagammalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Cthetagammaratio[Cthetagammaratio$thetagammaPeakPower > .2,], aes(x=age.x, y=thetagammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Cthetagammaratio[Cthetagammaratio$thetagammaduration > 2,], aes(x=age.x, y=thetagammaduration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=Cthetagammaratio[Cthetagammaratio$thetagammalogpower > 0,], thetagammalogpower ~ AgeInverse.x )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=Cthetagammaratio[Cthetagammaratio$thetagammaPeakPower > .2,], thetagammaPeakPower ~ AgeInverse.x )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=Cthetagammaratio[Cthetagammaratio$thetagammaduration > 2,], thetagammaduration ~ AgeInverse.x )
anova(lm.model)










#Beta Frontal lobe
FLBetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/FrontalLobe_Beta_Spectral_Analysis_Table_all.csv')
FLBetaPower$BetalogPower <- log(FLBetaPower$FrontalLobe_MeanAvgPower)
FLBetaPower$BetaAge <- (FLBetaPower$age)
FLBetaPower$BetaAvgEventNumber <- FLBetaPower$FrontalLobe_MeanAvgEventNumber
FLBetaPower$BetaAvgEventDuration <- FLBetaPower$FrontalLobe_MeanAvgEventDuration
FLBetaPower$BetaPeakPower <- log(FLBetaPower$FrontalLobe_MeanpeakPower)
FLBetaPower$BetaAgeInverse <- 1/(FLBetaPower$age)


ggplot(data=FLBetaPower[FLBetaPower$visitno < 2,], aes(x=BetaAge, y=BetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLBetaPower[FLBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLBetaPower[FLBetaPower$visitno < 2 , ], aes(x=BetaAge, y=BetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLBetaPower[FLBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=FLBetaPower[FLBetaPower$visitno < 2 ,], BetalogPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=FLBetaPower[FLBetaPower$visitno < 2 ,], BetaAvgEventNumber ~ BetaAgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=FLBetaPower[FLBetaPower$visitno < 2 ,], BetaAvgEventDuration ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=FLBetaPower[FLBetaPower$visitno < 2 ,], BetaPeakPower ~ BetaAgeInverse )
anova(lm.model)


#Beta parietal lobe
PLBetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ParietalLobe_Beta_Spectral_Analysis_Table_all.csv')
PLBetaPower$BetalogPower <- log(PLBetaPower$ParietalLobe_MeanAvgPower)
PLBetaPower$BetaAge <- (PLBetaPower$age)
PLBetaPower$BetaAvgEventNumber <- PLBetaPower$ParietalLobe_MeanAvgEventNumber
PLBetaPower$BetaAvgEventDuration <- PLBetaPower$ParietalLobe_MeanAvgEventDuration
PLBetaPower$BetaPeakPower <- log(PLBetaPower$ParietalLobe_MeanpeakPower)
PLBetaPower$BetaAgeInverse <- 1/(PLBetaPower$age)


ggplot(data=PLBetaPower[PLBetaPower$visitno < 2 & PLBetaPower$BetalogPower < 2.5,], aes(x=BetaAge, y=BetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLBetaPower[PLBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLBetaPower[PLBetaPower$visitno < 2 , ], aes(x=BetaAge, y=BetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLBetaPower[PLBetaPower$visitno < 2 & PLBetaPower$BetaPeakPower <5,], aes(x=BetaAge, y=BetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=PLBetaPower[PLBetaPower$visitno < 2 & PLBetaPower$BetalogPower < 2.5,], BetalogPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=PLBetaPower[PLBetaPower$visitno < 2 ,], BetaAvgEventNumber ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=PLBetaPower[PLBetaPower$visitno < 2 ,], BetaAvgEventDuration ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=PLBetaPower[PLBetaPower$visitno < 2 & PLBetaPower$BetaPeakPower <5,], BetaPeakPower ~ BetaAgeInverse )
anova(lm.model)




#Beta Occipital lobe
OLBetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/OccipitalLobe_Beta_Spectral_Analysis_Table_all.csv')
OLBetaPower$BetalogPower <- log(OLBetaPower$OccipitalLobe_MeanAvgPower)
OLBetaPower$BetaAge <- (OLBetaPower$age)
OLBetaPower$BetaAvgEventNumber <- OLBetaPower$OccipitalLobe_MeanAvgEventNumber
OLBetaPower$BetaAvgEventDuration <- OLBetaPower$OccipitalLobe_MeanAvgEventDuration
OLBetaPower$BetaPeakPower <- log(OLBetaPower$OccipitalLobe_MeanpeakPower)
OLBetaPower$BetaAgeInverse <- 1/(OLBetaPower$age)


ggplot(data=OLBetaPower[OLBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLBetaPower[OLBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLBetaPower[OLBetaPower$visitno < 2 , ], aes(x=BetaAge, y=BetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLBetaPower[OLBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=OLBetaPower[OLBetaPower$visitno < 2 ,], BetalogPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=OLBetaPower[OLBetaPower$visitno < 2 ,], BetaAvgEventNumber ~ BetaAgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=OLBetaPower[OLBetaPower$visitno < 2 ,], BetaAvgEventDuration ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=OLBetaPower[OLBetaPower$visitno < 2 ,], BetaPeakPower ~ BetaAgeInverse )
anova(lm.model)



#Beta central
CBetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Central_Beta_Spectral_Analysis_Table_all.csv')
CBetaPower$BetalogPower <- log(CBetaPower$Central_MeanAvgPower)
CBetaPower$BetaAge <- (CBetaPower$age)
CBetaPower$BetaAvgEventNumber <- CBetaPower$Central_MeanAvgEventNumber
CBetaPower$BetaAvgEventDuration <- CBetaPower$Central_MeanAvgEventDuration
CBetaPower$BetaPeakPower <- log(CBetaPower$Central_MeanpeakPower)
CBetaPower$BetaAgeInverse <- 1/(CBetaPower$age)


ggplot(data=CBetaPower[CBetaPower$visitno < 2 & CBetaPower$BetalogPower < 2.4,], aes(x=BetaAge, y=BetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CBetaPower[CBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CBetaPower[CBetaPower$visitno < 2 , ], aes(x=BetaAge, y=BetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CBetaPower[CBetaPower$visitno < 2 ,], aes(x=BetaAge, y=BetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=CBetaPower[CBetaPower$visitno < 2 & CBetaPower$BetalogPower < 2.4,], BetalogPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=CBetaPower[CBetaPower$visitno < 2 ,], BetaAvgEventNumber ~ BetaAgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=CBetaPower[CBetaPower$visitno < 2 ,], BetaAvgEventDuration ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=CBetaPower[CBetaPower$visitno < 2 ,], BetaPeakPower ~ BetaAgeInverse )
anova(lm.model)




#theta beta ratio frontal 
FLthetabetaratio <- merge(FLBetaPower[FLBetaPower$visitno < 2,], FLThetaPower , by='idvalues')
FLthetabetaratio$thetaBetalogpower<- FLthetabetaratio$ThetalogPower / FLthetabetaratio$BetalogPower
FLthetabetaratio$thetaBetaPeakPower <- FLthetabetaratio$ThetaPeakPower/ FLthetabetaratio$BetaPeakPower
FLthetabetaratio$thetaBetaDuration <- FLthetabetaratio$thetaAvgEventDuration/ FLthetabetaratio$BetaAvgEventDuration


ggplot(data=FLthetabetaratio[], aes(x=age.x, y=thetaBetalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLthetabetaratio, aes(x=age.x, y=thetaBetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=FLthetabetaratio[], aes(x=age.x, y=thetaBetaDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=FLthetabetaratio[], thetaBetalogpower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=FLthetabetaratio, thetaBetaPeakPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=FLthetabetaratio[], thetaBetaDuration ~ BetaAgeInverse )
anova(lm.model)

#theta beta ratio parietal
PLthetaBetaratio <- merge(PLBetaPower[PLBetaPower$visitno < 2,], PLThetaPower , by='idvalues')
PLthetaBetaratio$thetaBetalogpower<- PLthetaBetaratio$ThetalogPower / PLthetaBetaratio$BetalogPower
PLthetaBetaratio$thetaBetaPeakPower <- PLthetaBetaratio$ThetaPeakPower/ PLthetaBetaratio$BetaPeakPower
PLthetaBetaratio$thetaBetaduration <- PLthetaBetaratio$thetaAvgEventDuration/ PLthetaBetaratio$BetaAvgEventDuration


ggplot(data=PLthetaBetaratio[], aes(x=age.x, y=thetaBetalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLthetaBetaratio[], aes(x=age.x, y=thetaBetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=PLthetaBetaratio, aes(x=age.x, y=thetaBetaduration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=PLthetaBetaratio, thetaBetalogpower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=PLthetaBetaratio[], thetaBetaPeakPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=PLthetaBetaratio, thetaBetaduration ~ BetaAgeInverse )
anova(lm.model)



#theta Beta ratio occipital
OLthetaBetaratio <- merge(OLBetaPower[OLBetaPower$visitno < 2,], OLThetaPower , by='idvalues')
OLthetaBetaratio$thetaBetalogpower<- OLthetaBetaratio$ThetalogPower / OLthetaBetaratio$BetalogPower
OLthetaBetaratio$thetaBetaPeakPower <- OLthetaBetaratio$ThetaPeakPower/ OLthetaBetaratio$BetaPeakPower
OLthetaBetaratio$thetaBetaduration <- OLthetaBetaratio$thetaAvgEventDuration/ OLthetaBetaratio$BetaAvgEventDuration


ggplot(data=OLthetaBetaratio[], aes(x=age.x, y=thetaBetalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLthetaBetaratio, aes(x=age.x, y=thetaBetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=OLthetaBetaratio, aes(x=age.x, y=thetaBetaduration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for log power
lm.model <- lm(data=OLthetaBetaratio[], thetaBetalogpower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=OLthetaBetaratio, thetaBetaPeakPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=OLthetaBetaratio, thetaBetaduration ~ BetaAgeInverse )
anova(lm.model)


#theta Beta ratio Central
CthetaBetaratio <- merge(CBetaPower[CBetaPower$visitno < 2,], CThetaPower , by='idvalues')
CthetaBetaratio$thetaBetalogpower<- CthetaBetaratio$ThetalogPower / CthetaBetaratio$BetalogPower
CthetaBetaratio$thetaBetaPeakPower <- CthetaBetaratio$ThetaPeakPower/ CthetaBetaratio$BetaPeakPower
CthetaBetaratio$thetaBetaduration <- CthetaBetaratio$thetaAvgEventDuration/ CthetaBetaratio$BetaAvgEventDuration


ggplot(data=CthetaBetaratio[], aes(x=age.x, y=thetaBetalogpower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CthetaBetaratio[], aes(x=age.x, y=thetaBetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=CthetaBetaratio[], aes(x=age.x, y=thetaBetaduration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#stats for log power
lm.model <- lm(data=CthetaBetaratio[], thetaBetalogpower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg peak power
lm.model <- lm(data=CthetaBetaratio[], thetaBetaPeakPower ~ BetaAgeInverse )
anova(lm.model)

#stats for avg event duration
lm.model <- lm(data=CthetaBetaratio[], thetaBetaduration ~ BetaAgeInverse )
anova(lm.model)
