if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2,lme4, multcomp, tidyr, tidyverse,wesanderson, reshape2, corrplot, ggExtra, cowplot, ggsci, mgcv, ggpointdensity, devtools, ggpubr, itsadug, useful, rmarkdown, tinytex, gratia, gamm4)
devtools::install_github('LabneuroCogDevel/LNCDR')
library(LNCDR)

#gamma all channels
gammaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_all_allevents.csv')
gammaPower$GammalogPower <- log(gammaPower$MeanAvgPower)
gammaPower$GammaAge <- (gammaPower$age)
gammaPower$GammaAvgEventNumber <- gammaPower$MeanAvgEventNumber
gammaPower$GammaAvgEventDuration <- gammaPower$MeanAvgEventDuration
gammaPower$GammaAvgEventFSpan <- gammaPower$MeanAvgEventFSpan
gammaPower$GammaPeakPower <- log(gammaPower$MeanpeakPower)
gammaPower$AgeInverse <- 1/(gammaPower$age)


ggplot(data=gammaPower[gammaPower$visitno < 2,], aes(x=GammaAge, y=GammalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPower[gammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPower[gammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPower[gammaPower$visitno < 2,], aes(x=GammaAge, y=GammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#compare age to inverse age
lm.lin <- lm(GammalogPower ~ GammaAge, data=gammaPower[gammaPower$visitno < 2,])
lm.inv <- lm(GammalogPower ~ AgeInverse, data=gammaPower[gammaPower$visitno < 2,])
summary(lm.lin)
summary(lm.inv)
AIC(lm.lin, lm.inv)

#stats for avg peak power
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2 ,], GammaPeakPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2 ,], GammaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2,], GammaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for log power
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2 ,], GammalogPower ~ AgeInverse )
anova(lm.model)





#gamma DLPFC
gammaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Gamma_Spectral_Analysis_Table_all_allevents.csv')
gammaPower$GammalogPower <- log(gammaPower$DLPFC_AvgPower)
gammaPower$GammaAge <- (gammaPower$age)
gammaPower$GammaAvgEventNumber <- gammaPower$DLPFC_AvgEventNumber
gammaPower$GammaAvgEventDuration <- gammaPower$DLPFC_AvgEventDuration
gammaPower$GammaPeakPower <- log(gammaPower$DLPFC_peakPower)
gammaPower$AgeInverse <- 1/(gammaPower$age)


ggplot(data=gammaPower[gammaPower$visitno < 2,], aes(x=GammaAge, y=GammalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPower[gammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPower[gammaPower$visitno < 2 & gammaPower$GammaAvgEventDuration < 0.029, ], aes(x=GammaAge, y=GammaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPower[gammaPower$visitno < 2,], aes(x=GammaAge, y=GammaAvgEventFSpan)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPower[gammaPower$visitno < 2 ,], aes(x=GammaAge, y=GammaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for avg peak power
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2 ,], GammaPeakPower ~ AgeInverse )
anova(lm.model)

#stats for avg event number
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2 ,], GammaAvgEventNumber ~ AgeInverse )
anova(lm.model)

#stats for avg eevent duration
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2 & gammaPower$GammaAvgEventDuration < 0.029,], GammaAvgEventDuration ~ AgeInverse )
anova(lm.model)

#stats for log power
lm.model <- lm(data=gammaPower[gammaPower$visitno < 2 ,], GammalogPower ~ AgeInverse )
anova(lm.model)


#compared fieldtrip and stephanie toolboxes
fieldtrip <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Fieldtrip_Table_all.csv')
fieldtrip$fieldtrippower <- log(fieldtrip$Var2)

bothDataSets <- merge(gammaPower[gammaPower$visitno < 2,], fieldtrip, by='idvalues')
ggplot(data=bothDataSets, aes(x=bothDataSets$logPower, y=bothDataSets$fieldtrippower)) + geom_point() + stat_smooth(method='lm')
ggplot(data=bothDataSets, aes(x=bothDataSets$AvgPower, y=bothDataSets$Var2)) + geom_point() + stat_smooth(method='lm')






#theta all channels 
thetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_all_allevents.csv')
thetaPower$thetalogPower <- log(thetaPower$MeanAvgPower)
thetaPower$thetaAge <- (thetaPower$age)
thetaPower$thetaAvgEventNumber <- thetaPower$MeanAvgEventNumber
thetaPower$thetaAvgEventDuration <- thetaPower$MeanAvgEventDuration
thetaPower$thetaAvgEventFSpan <- thetaPower$MeanAvgEventFSpan
thetaPower$thetaPeakPower <- thetaPower$MeanpeakPower
thetaPower$thetaAgeInverse <- 1/(thetaPower$age)


ggplot(data=thetaPower[thetaPower$visitno < 2 & thetaPower$thetalogPower < 2 & thetaPower$thetalogPower > -1.0,], aes(x=thetaAge, y=thetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=thetaPower[thetaPower$visitno < 2 ,], aes(x=thetaAge, y=thetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=thetaPower[thetaPower$visitno < 2,], aes(x=thetaAge, y=thetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=thetaPower[thetaPower$visitno < 2,], aes(x=thetaAge, y=thetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#theta stats
lm.model <- lm(data=thetaPower[thetaPower$visitno < 2 & thetaPower$thetalogPower > -1,], thetalogPower ~ thetaAgeInverse )
anova(lm.model)


lm.model <- lm(data=thetaPower[thetaPower$visitno < 2,], thetaAvgEventNumber ~ thetaAgeInverse )
anova(lm.model)

lm.model <- lm(data=thetaPower[thetaPower$visitno < 2,], thetaAvgEventDuration ~ thetaAgeInverse )
anova(lm.model)

lm.model <- lm(data=thetaPower[thetaPower$visitno < 2 & thetaPower$thetaPeakPower < 60,], thetaPeakPower ~ thetaAgeInverse )
anova(lm.model)






#theta DLPFC
thetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Theta_Spectral_Analysis_Table_all_allevents.csv')
thetaPower$thetalogPower <- log(thetaPower$DLPFC_AvgPower)
thetaPower$thetaAge <- (thetaPower$age)
thetaPower$thetaAvgEventNumber <- thetaPower$DLPFC_AvgEventNumber
thetaPower$thetaAvgEventDuration <- thetaPower$DLPFC_AvgEventDuration
thetaPower$thetaPeakPower <- thetaPower$DLPFC_peakPower
thetaPower$thetaAgeInverse <- 1/(thetaPower$age)


ggplot(data=thetaPower[thetaPower$visitno < 2 ,], aes(x=thetaAge, y=thetalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=thetaPower[thetaPower$visitno < 2 ,], aes(x=thetaAge, y=thetaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=thetaPower[thetaPower$visitno < 2,], aes(x=thetaAge, y=thetaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=thetaPower[thetaPower$visitno < 2 & thetaPower$thetaPeakPower < 60,], aes(x=thetaAge, y=thetaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#theta stats
lm.model <- lm(data=thetaPower[thetaPower$visitno < 2 ,], thetalogPower ~ thetaAgeInverse )
anova(lm.model)


lm.model <- lm(data=thetaPower[thetaPower$visitno < 2 ,], thetaAvgEventNumber ~ thetaAgeInverse )
anova(lm.model)

lm.model <- lm(data=thetaPower[thetaPower$visitno < 2 ,], thetaAvgEventDuration ~ thetaAgeInverse )
anova(lm.model)

lm.model <- lm(data=thetaPower[thetaPower$visitno < 2 & thetaPower$thetaPeakPower < 60,], thetaPeakPower ~ thetaAgeInverse )
anova(lm.model)





#theta gamma ratio- all channels
thetaGamma <- merge(gammaPower[gammaPower$visitno < 2,], thetaPower, by='idvalues')
thetaGamma$thetaGammaPowerRatio <- thetaGamma$thetalogPower/thetaGamma$GammalogPower
thetaGamma$thetaGammaInverseAge <- 1/thetaGamma$age.x
ggplot(data = thetaGamma[thetaGamma$thetaGammaPowerRatio< 0.5 & thetaGamma$thetaGammaPowerRatio > -0.2 , ], aes(x = age.x, y = thetaGammaPowerRatio )) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

thetaGamma$thetaGammaPeakPower <- thetaGamma$thetaPeakPower/thetaGamma$GammaPeakPower
ggplot(data = thetaGamma[thetaGamma$thetaGammaPeakPower<20,], aes(x = age.x, y = thetaGammaPeakPower )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


thetaGamma$thetaGammalengthRatio <- thetaGamma$thetaAvgEventDuration/thetaGamma$GammaAvgEventDuration
ggplot(data = thetaGamma[thetaGamma$thetaGammalengthRatio > 2,], aes(x = age.x, y = thetaGammalengthRatio )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


#theta gamma stats
lm.model <- lm(data=thetaGamma, thetaGammaPowerRatio ~ thetaGammaInverseAge)
anova(lm.model)

lm.model <- lm(data=thetaGamma[thetaGamma$thetaGammaPeakPower<20,], thetaGammaPeakPower ~ thetaGammaInverseAge)
anova(lm.model)


lm.model <- lm(data=thetaGamma[thetaGamma$thetaGammalengthRatio > 2,], thetaGammalengthRatio ~ thetaGammaInverseAge)
anova(lm.model)






#theta gamma ratio- DLPFC
thetaGamma <- merge(gammaPower[gammaPower$visitno < 2,], thetaPower, by='idvalues')
thetaGamma$thetaGammaPowerRatio <- thetaGamma$thetalogPower/thetaGamma$GammalogPower
thetaGamma$thetaGammaInverseAge <- 1/thetaGamma$age.x
ggplot(data = thetaGamma, aes(x = age.x, y = thetaGammaPowerRatio )) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

thetaGamma$thetaGammaPeakPower <- thetaGamma$thetaPeakPower/thetaGamma$GammaPeakPower
ggplot(data = thetaGamma[thetaGamma$thetaGammaPeakPower<20,], aes(x = age.x, y = thetaGammaPeakPower )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


thetaGamma$thetaGammalengthRatio <- thetaGamma$thetaAvgEventDuration/thetaGamma$GammaAvgEventDuration
ggplot(data = thetaGamma, aes(x = age.x, y = thetaGammalengthRatio )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


#theta gamma stats
lm.model <- lm(data=thetaGamma, thetaGammaPowerRatio ~ thetaGammaInverseAge)
anova(lm.model)

lm.model <- lm(data=thetaGamma[thetaGamma$thetaGammaPeakPower<20,], thetaGammaPeakPower ~ thetaGammaInverseAge)
anova(lm.model)


lm.model <- lm(data=thetaGamma, thetaGammalengthRatio ~ thetaGammaInverseAge)
anova(lm.model)



#beta all channels 
betaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_all_allevents.csv')
betaPower$betalogPower <- log(betaPower$MeanAvgPower)
betaPower$betaAge <- (betaPower$age)
betaPower$betaAvgEventNumber <- betaPower$MeanAvgEventNumber
betaPower$betaAvgEventDuration <- betaPower$MeanAvgEventDuration
betaPower$betaPeakPower <- betaPower$MeanpeakPower
betaPower$betaAgeInverse <- 1/(betaPower$age)


ggplot(data=betaPower[betaPower$visitno < 2,], aes(x=betaAge, y=betalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=betaPower[betaPower$visitno < 2,], aes(x=betaAge, y=betaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=betaPower[betaPower$visitno < 2,], aes(x=betaAge, y=betaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=betaPower[betaPower$visitno < 2,], aes(x=betaAge, y=betaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#beta stats
lm.model <- lm(data=betaPower[betaPower$visitno < 2 ,], betalogPower ~ betaAgeInverse )
anova(lm.model)


lm.model <- lm(data=betaPower[betaPower$visitno < 2,], betaAvgEventNumber ~ betaAgeInverse )
anova(lm.model)

lm.model <- lm(data=betaPower[betaPower$visitno < 2 ,], betaAvgEventDuration ~ betaAgeInverse )
anova(lm.model)

lm.model <- lm(data=betaPower[betaPower$visitno < 2,], betaPeakPower ~ betaAgeInverse )
anova(lm.model)


#beta DLPFC
DLPFCbetaPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Beta_Spectral_Analysis_Table_all_allevents.csv')
DLPFCbetaPower$betalogPower <- log(DLPFCbetaPower$DLPFC_AvgPower)
DLPFCbetaPower$betaAge <- (DLPFCbetaPower$age)
DLPFCbetaPower$betaAvgEventNumber <- DLPFCbetaPower$DLPFC_AvgEventNumber
DLPFCbetaPower$betaAvgEventDuration <- DLPFCbetaPower$DLPFC_AvgEventDuration
DLPFCbetaPower$betaPeakPower <- DLPFCbetaPower$DLPFC_peakPower
DLPFCbetaPower$betaAgeInverse <- 1/(DLPFCbetaPower$age)


ggplot(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2,], aes(x=betaAge, y=betalogPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2,], aes(x=betaAge, y=betaAvgEventNumber)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2,], aes(x=betaAge, y=betaAvgEventDuration)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2 & DLPFCbetaPower$betaPeakPower < 100,], aes(x=betaAge, y=betaPeakPower)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

#beta stats
lm.model <- lm(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2 ,], betalogPower ~ betaAgeInverse )
anova(lm.model)


lm.model <- lm(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2,], betaAvgEventNumber ~ betaAgeInverse )
anova(lm.model)

lm.model <- lm(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2 ,], betaAvgEventDuration ~ betaAgeInverse )
anova(lm.model)

lm.model <- lm(data=DLPFCbetaPower[DLPFCbetaPower$visitno < 2 & DLPFCbetaPower$betaPeakPower < 100,], betaPeakPower ~ betaAgeInverse )
anova(lm.model)





#theta beta ratio- all channels
thetaBeta <- merge(betaPower[betaPower$visitno < 2,], thetaPower, by='idvalues')

thetaBeta$thetaBetaPowerRatio <- thetaBeta$thetalogPower/thetaBeta$betalogPower
thetaBeta$thetaBetaInverseAge <- 1/thetaBeta$age.x
ggplot(data = thetaBeta[], aes(x = age.x, y = thetaBetaPowerRatio )) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

thetaBeta$thetaBetaPeakPower <- thetaBeta$thetaPeakPower/thetaBeta$betaPeakPower
ggplot(data = thetaBeta[], aes(x = age.x, y = thetaBetaPeakPower )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


thetaBeta$thetaBetalengthRatio <- thetaBeta$thetaAvgEventDuration/thetaBeta$betaAvgEventDuration
ggplot(data = thetaBeta[], aes(x = age.x, y = thetaBetalengthRatio )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


#theta gamma stats
lm.model <- lm(data=thetaBeta, thetaBetaPowerRatio ~ thetaBetaInverseAge)
anova(lm.model)

lm.model <- lm(data=thetaBeta[], thetaBetaPeakPower ~ thetaBetaInverseAge)
anova(lm.model)


lm.model <- lm(data=thetaBeta[], thetaBetalengthRatio ~ thetaBetaInverseAge)
anova(lm.model)




#theta beta ratio- DLPFC
thetaBeta <- merge(betaPower[betaPower$visitno < 2,], thetaPower, by='idvalues')

thetaBeta$thetaBetaPowerRatio <- thetaBeta$thetalogPower/thetaBeta$betalogPower
thetaBeta$thetaBetaInverseAge <- 1/thetaBeta$age.x
ggplot(data = thetaBeta[], aes(x = age.x, y = thetaBetaPowerRatio )) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')

thetaBeta$thetaBetaPeakPower <- thetaBeta$thetaPeakPower/thetaBeta$betaPeakPower
ggplot(data = thetaBeta[thetaBeta$thetaBetaPeakPower <4,], aes(x = age.x, y = thetaBetaPeakPower )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


thetaBeta$thetaBetalengthRatio <- thetaBeta$thetaAvgEventDuration/thetaBeta$betaAvgEventDuration
ggplot(data = thetaBeta[thetaBeta$thetaBetalengthRatio < 2.5,], aes(x = age.x, y = thetaBetalengthRatio )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


#theta gamma stats
lm.model <- lm(data=thetaBeta, thetaBetaPowerRatio ~ thetaBetaInverseAge)
anova(lm.model)

lm.model <- lm(data=thetaBeta[thetaBeta$thetaBetaPeakPower <4,], thetaBetaPeakPower ~ thetaBetaInverseAge)
anova(lm.model)


lm.model <- lm(data=thetaBeta[thetaBeta$thetaBetalengthRatio < 2.5,], thetaBetalengthRatio ~ thetaBetaInverseAge)
anova(lm.model)



# Gamma Avg Power SD 
GammaSDPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_SDPower_all.csv')
GammaSDPower$meanSD <- mean(GammaSDPower$SDPower)
GammaSDPower$SD <- sd(GammaSDPower$SDPower)
GammaSDPower$gammaInverseAge <- 1/ GammaSDPower$age

ggplot(data = GammaSDPower[GammaSDPower$visitno < 2 & GammaSDPower$SDPower < 2.5* GammaSDPower$SD,], aes(x = age, y = SDPower))+ geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data = GammaSDPower[GammaSDPower$visitno < 2 & GammaSDPower$SDPower < 2.5* GammaSDPower$SD,], SDPower ~ gammaInverseAge)
anova(lm.model)



# Theta Avg Power SD 
ThetaSDPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_SD_Power.csv')
ThetaSDPower$meanSD <- mean(ThetaSDPower$SDPower)
ThetaSDPower$SD <- sd(ThetaSDPower$SDPower)


ggplot(data = ThetaSDPower[ThetaSDPower$visitno < 2 & ThetaSDPower$SDPower < 2.5* ThetaSDPower$SD,], aes(x = age, y = SDPower))+ geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')


#Beta Avg Power SD
BetaSDPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_SD_Power.csv')
BetaSDPower$meanSD <- mean(BetaSDPower$SDPower)
BetaSDPower$SD <- sd(BetaSDPower$SDPower)
BetaSDPower$betaInverseAge <- 1/ BetaSDPower$age

ggplot(data = BetaSDPower[BetaSDPower$visitno < 2 & BetaSDPower$SDPower < 2.5* BetaSDPower$SD,], aes(x = age, y = SDPower))+ geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data = BetaSDPower[BetaSDPower$visitno < 2 & BetaSDPower$SDPower < 2.5* BetaSDPower$SD,], SDPower ~ betaInverseAge)
anova(lm.model)


