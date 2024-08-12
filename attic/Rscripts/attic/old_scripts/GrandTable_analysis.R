if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2,lme4, multcomp, tidyr, tidyverse,wesanderson, reshape2, corrplot, ggExtra, cowplot, ggsci, mgcv, ggpointdensity, devtools, ggpubr, itsadug, useful, rmarkdown, tinytex, gratia, gamm4, car, GGally, CCA)
devtools::install_github('LabneuroCogDevel/LNCDR')
library(devtools)
devtools::install_github("strengejacke/sjPlot")
library(LNCDR)
library(data.table)


GrandTable <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/GrandTable_withAges.csv')
GammaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_all_data.csv')
AlphaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_all_data.csv')
BetaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_all_data.csv')
ThetaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_all_data.csv')
agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
agefile$Subject <- agefile$idvalues 

GrandTable <- merge(GammaT, AlphaT, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
GrandTable <- merge(GrandTable, BetaT, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
GrandTable <- merge(GrandTable, ThetaT, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
GrandTable <- merge(GrandTable, agefile, by = 'Subject', all.x = TRUE, all.y = TRUE)

GrandTable$log_Gamma_Trial_Power_mean <- log(GrandTable$Gamma_Trial_Power_Mean)
GrandTable$log_Gamma_Trial_Power_SD <- log(GrandTable$Gamma_Trial_Power_SD)
GrandTable$log_Theta_Trial_Power_SD <- log(GrandTable$Theta_Trial_Power_SD)
GrandTable$log_Theta_Trial_Power_mean <- log(GrandTable$Theta_Trial_Power_Mean)
GrandTable$log_Beta_Trial_Power_mean <- log(GrandTable$Beta_Trial_Power_Mean)
GrandTable$log_Beta_Trial_Power_SD <- log(GrandTable$Beta_Trial_Power_SD)
GrandTable$log_Alpha_Trial_Power_SD <- log(GrandTable$Alpha_Trial_Power_SD)
GrandTable$log_Alpha_Trial_Power_mean <- log(GrandTable$Alpha_Trial_Power_Mean)



Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior.csv')

Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
Behavior <- merge(Behavior, agefile, by = 'Subject', all.x = TRUE, all.y = TRUE)

Everything <- merge(GrandTable, Behavior, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
Everything$absPositionError <- abs(Everything$PositionError)


## looking for effects between behavior and eeg measures controlling for age 
names(Everything)
varstobescale<-c("log_Gamma_Trial_Power_mean","Gamma_Event_Number_Mean","Gamma_Event_Duration_Mean")
#Everything[Everything$visitno.x < 2 & Everything$log_Gamma_Trial_Power_mean != '-Inf',varstobescale]<-lapply(Everything[Everything$visitno.x < 2 & Everything$log_Gamma_Trial_Power_mean != '-Inf',varstobescale],function(x),scale)
modelFit <- lmer(log_Gamma_Trial_Power_mean ~ age.x + Gamma_Event_Number_Mean + Gamma_Event_Duration_Mean  + (1|Subject), data = Everything[Everything$visitno.x < 2 & Everything$log_Gamma_Trial_Power_mean != '-Inf',])
car::Anova(modelFit)


require(ggplot2)
V1data<-Everything[Everything$visitno.x < 2 & Everything$log_Gamma_Trial_Power_mean != '-Inf',]
ggplot(V1data,aes(x=Gamma_Event_Number_Mean,y=log_Gamma_Trial_Power_mean,colour=Gamma_Event_Duration_Mean))+geom_point()
hist()

summary(modelFit)


ggplot(data=  Everything[Everything$visitno.x < 2,], aes(x = Beta_Event_Duration_Mean, y = mgsLatency))  + stat_smooth(method = 'lm')
ggdensity(data = Everything[Everything$visitno.x < 2,], x = 'mgsLatency', y = '..density..')



ggplot(data=  Everything[Everything$visitno.x < 2 ,], aes(x = Beta_Event_Duration_Mean, y = mgsLatency, color = factor(Group.x) ))  +stat_smooth(method = 'lm', se = FALSE) 


#average everyones values
Avg_absPositionError <- aggregate(absPositionError ~ Subject, Everything, mean)
BetweenVar_absPositionError <- aggregate(absPositionError ~ Subject, Everything, sd)

Avg_absPositionError <- merge(Avg_absPositionError, agefile, by = 'Subject')
Avg_Gamma_Event_Duration_Mean <- aggregate(Gamma_Event_Duration_Mean ~ Subject, Everything, mean)
Avg_Gamma_Event_Duration_Mean <- merge(Avg_Gamma_Event_Duration_Mean, agefile, by = 'Subject')
Gamma_event_duration_positionError <- merge(Avg_absPositionError, Avg_Gamma_Event_Duration_Mean, by= 'Subject' )

Avg_Gamma_Event_Duration_BetweenVar <- aggregate(Gamma_Event_Duration_Mean ~ Subject, Everything, sd)
Avg_Gamma_Event_Duration_BetweenVar <- merge(Avg_Gamma_Event_Duration_BetweenVar, agefile, by = 'Subject')
Gamma_event_duration_betweenVar_positionError <- merge(BetweenVar_absPositionError, Avg_Gamma_Event_Duration_BetweenVar, by= 'Subject' )

ggplot(data = Gamma_event_duration_betweenVar_positionError[Gamma_event_duration_betweenVar_positionError$visitno <2,], aes(x = Gamma_Event_Duration_Mean, y =absPositionError)) + geom_point(aes(color = age)) + scale_color_gradient2(low = "blue", high = "red") + stat_smooth(method='lm')



Avg_mgsLatency <- aggregate(mgsLatency ~ Subject, Everything, mean)
Avg_mgsLatency <- merge(Avg_mgsLatency, agefile, by = 'Subject')
Avg_Gamma_Event_Duration_Mean <- aggregate(Gamma_Event_Duration_Mean ~ Subject, Everything, mean)
Avg_Gamma_Event_Duration_Mean <- merge(Avg_Gamma_Event_Duration_Mean, agefile, by = 'Subject')
Gamma_event_duration_mgsLatency <- merge(Avg_mgsLatency, Avg_Gamma_Event_Duration_Mean, by= 'Subject' )


Avg_Beta_Event_Duration_Mean <- aggregate(Beta_Event_Duration_Mean ~ Subject, Everything, mean)
Avg_Beta_Event_Duration_Mean <- merge(Avg_Beta_Event_Duration_Mean, agefile, by = 'Subject')
Beta_event_duration_positionError <- merge(Avg_absPositionError, Avg_Beta_Event_Duration_Mean, by= 'Subject' )
Beta_event_duration_mgsLatency <- merge(Avg_mgsLatency, Avg_Beta_Event_Duration_Mean, by= 'Subject' )


Avg_Beta_Event_Number_Mean <- aggregate(Beta_Event_Number_Mean ~ Subject, Everything, mean)
Avg_Beta_Event_Number_Mean <- merge(Avg_Beta_Event_Number_Mean, agefile, by = 'Subject')
Beta_event_number_positionError <- merge(Avg_absPositionError, Avg_Beta_Event_Number_Mean, by= 'Subject' )


Avg_Alpha_trial_power_SD <- aggregate(log_Alpha_Trial_Power_SD ~ Subject, Everything[Everything$log_Alpha_Trial_Power_SD != '-Inf',], mean)
Avg_Alpha_trial_power_SD <- merge(Avg_Alpha_trial_power_SD, agefile, by = 'Subject')
Alpha_Trial_power_SD_mgsLatency <- merge(Avg_mgsLatency, Avg_Alpha_trial_power_SD, by= 'Subject' )
ggplot(data = Alpha_Trial_power_SD_mgsLatency[Alpha_Trial_power_SD_mgsLatency$visitno.x <2,], aes(x = log_Alpha_Trial_Power_SD, y =mgsLatency)) + geom_point(aes(color = age.x)) + scale_color_gradient2(low = "blue", high = "red") + stat_smooth(method='lm')


ggplot(data = Gamma_event_duration_positionError[Gamma_event_duration_positionError$visitno.x < 2 & Gamma_event_duration_positionError$absPositionError <10,], aes(x = Gamma_Event_Duration_Mean, y =absPositionError)) + geom_point(aes(color = age.x)) + scale_color_gradient2(low = "blue", high = "red") + stat_smooth(method='lm')
ggplot(data = Gamma_event_duration_mgsLatency[Gamma_event_duration_mgsLatency$visitno.x <2,], aes(x = Gamma_Event_Duration_Mean, y =mgsLatency)) + geom_point(aes(color = age.x)) + scale_color_gradient2(low = "blue", high = "red") + stat_smooth(method='lm')

ggplot(data = Beta_event_duration_mgsLatency[Beta_event_duration_mgsLatency$visitno.x <2,], aes(x = Beta_Event_Duration_Mean, y =mgsLatency)) + geom_point(aes(color = age.x)) + scale_color_gradient2(low = "blue", high = "red") + stat_smooth(method='lm')
ggplot(data = Beta_event_number_positionError[Beta_event_number_positionError$visitno.x <2 & Beta_event_number_positionError$absPositionError <10,], aes(x = Beta_Event_Number_Mean, y =absPositionError)) + geom_point(aes(color = age.x)) + scale_color_gradient2(low = "blue", high = "red") + stat_smooth(method='lm')
ggplot(data = Beta_event_duration_positionError[Beta_event_duration_positionError$visitno.x< 2 & Beta_event_duration_positionError$absPositionError <10,], aes(x = Beta_Event_Duration_Mean, y =absPositionError)) + geom_point(aes(color = age.x)) + scale_color_gradient2(low = "blue", high = "red") + stat_smooth(method='lm')



lm.model <- lm(data=Gamma_event_duration_mgsLatency[Gamma_event_duration_mgsLatency$visitno.x <2,], Gamma_Event_Duration_Mean ~ mgsLatency)
anova(lm.model)

lm.model <- lm(data=Beta_event_number_positionError[Beta_event_number_positionError$absPositionError <10,], Beta_Event_Number_Mean ~ absPositionError)
anova(lm.model)



## Look into alpha values 
#trial power sd
Avg_Alpha_trial_power_SD <- aggregate(log_Alpha_Trial_Power_SD ~ Subject, Everything[Everything$log_Alpha_Trial_Power_SD != '-Inf',], mean)
Avg_Alpha_trial_power_SD <- merge(Avg_Alpha_trial_power_SD, agefile, by = 'Subject')
Avg_Alpha_trial_power_SD$ageInverse <- 1/Avg_Alpha_trial_power_SD$age
ggplot(data = Avg_Alpha_trial_power_SD[Avg_Alpha_trial_power_SD$visitno <2,], aes(x = age, y = log_Alpha_Trial_Power_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Alpha_trial_power_SD[Avg_Alpha_trial_power_SD$visitno < 2,], log_Alpha_Trial_Power_SD ~ ageInverse)
anova(lm.model)

#trial power mean
Avg_Alpha_trial_power_mean <- aggregate(log_Alpha_Trial_Power_mean ~ Subject, Everything[Everything$log_Alpha_Trial_Power_mean != '-Inf',], mean)
Avg_Alpha_trial_power_mean <- merge(Avg_Alpha_trial_power_mean, agefile, by = 'Subject')
Avg_Alpha_trial_power_mean$ageInverse <- 1/Avg_Alpha_trial_power_mean$age
ggplot(data = Avg_Alpha_trial_power_mean[Avg_Alpha_trial_power_mean$visitno < 2,], aes(x = age, y = log_Alpha_Trial_Power_mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Alpha_trial_power_mean[Avg_Alpha_trial_power_mean$visitno < 2,], log_Alpha_Trial_Power_mean ~ ageInverse)
anova(lm.model)

#event number 
Avg_Alpha_Event_Number <- aggregate(Alpha_Event_Number_Mean ~ Subject, Everything, mean)
Avg_Alpha_Event_Number <- merge(Avg_Alpha_Event_Number, agefile, by = 'Subject')
Avg_Alpha_Event_Number$ageInverse <- 1/Avg_Alpha_Event_Number$age

ggplot(data = Avg_Alpha_Event_Number[Avg_Alpha_Event_Number$visitno <2,], aes(x = age, y = Alpha_Event_Number_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Alpha_Event_Number[Avg_Alpha_Event_Number$visitno <2,], Alpha_Event_Number_Mean ~ ageInverse)
anova(lm.model)


Avg_Alpha_Event_Number_SD <- aggregate(Alpha_Event_Number_SD ~ Subject, Everything, mean)
Avg_Alpha_Event_Number_SD <- merge(Avg_Alpha_Event_Number_SD, agefile, by = 'Subject')
Avg_Alpha_Event_Number_SD$ageInverse <- 1/Avg_Alpha_Event_Number_SD$age

ggplot(data = Avg_Alpha_Event_Number_SD[Avg_Alpha_Event_Number_SD$visitno <2,], aes(x = age, y = Alpha_Event_Number_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Alpha_Event_Number_SD[Avg_Alpha_Event_Number_SD$visitno <2,], Alpha_Event_Number_SD ~ ageInverse)
anova(lm.model)


#event duration 
Avg_Alpha_Event_Duration <- aggregate(Alpha_Event_Duration_Mean ~ Subject, Everything, mean)
Avg_Alpha_Event_Duration <- merge(Avg_Alpha_Event_Duration, agefile, by = 'Subject')
Avg_Alpha_Event_Duration$ageInverse <- 1/Avg_Alpha_Event_Duration$age

ggplot(data = Avg_Alpha_Event_Duration[Avg_Alpha_Event_Duration$visitno <2,], aes(x = age, y = Alpha_Event_Duration_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Alpha_Event_Duration[Avg_Alpha_Event_Duration$visitno <2,], Alpha_Event_Duration_Mean ~ ageInverse)
anova(lm.model)


Avg_Alpha_Event_Duration_SD <- aggregate(Alpha_Event_Duration_SD ~ Subject, Everything, mean)
Avg_Alpha_Event_Duration_SD <- merge(Avg_Alpha_Event_Duration_SD, agefile, by = 'Subject')
Avg_Alpha_Event_Duration_SD$ageInverse <- 1/Avg_Alpha_Event_Duration_SD$age

ggplot(data = Avg_Alpha_Event_Duration_SD[Avg_Alpha_Event_Duration_SD$visitno <2,], aes(x = age, y = Alpha_Event_Duration_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Alpha_Event_Duration_SD[Avg_Alpha_Event_Duration_SD$visitno <2,], Alpha_Event_Duration_SD ~ ageInverse)
anova(lm.model)



## Look into Gamma values 
#trial power sd
Avg_Gamma_trial_power_SD <- aggregate(log_Gamma_Trial_Power_SD ~ Subject, Everything[Everything$log_Gamma_Trial_Power_SD != '-Inf',], mean)
Avg_Gamma_trial_power_SD <- merge(Avg_Gamma_trial_power_SD, agefile, by = 'Subject')
Avg_Gamma_trial_power_SD$ageInverse <- 1/Avg_Gamma_trial_power_SD$age
ggplot(data = Avg_Gamma_trial_power_SD[Avg_Gamma_trial_power_SD$visitno <2,], aes(x = age, y = log_Gamma_Trial_Power_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Gamma_trial_power_SD[Avg_Gamma_trial_power_SD$visitno < 2,], log_Gamma_Trial_Power_SD ~ ageInverse)
anova(lm.model)

#trial power mean
Avg_Gamma_trial_power_mean <- aggregate(log_Gamma_Trial_Power_mean ~ Subject, Everything[Everything$log_Gamma_Trial_Power_mean != '-Inf',], mean)
Avg_Gamma_trial_power_mean <- merge(Avg_Gamma_trial_power_mean, agefile, by = 'Subject')
Avg_Gamma_trial_power_mean$ageInverse <- 1/Avg_Gamma_trial_power_mean$age
ggplot(data = Avg_Gamma_trial_power_mean[Avg_Gamma_trial_power_mean$visitno < 2,], aes(x = age, y = log_Gamma_Trial_Power_mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Gamma_trial_power_mean[Avg_Gamma_trial_power_mean$visitno < 2,], log_Gamma_Trial_Power_mean ~ ageInverse)
anova(lm.model)

#event number 
Avg_Gamma_Event_Number <- aggregate(Gamma_Event_Number_Mean ~ Subject, Everything, mean)
Avg_Gamma_Event_Number <- merge(Avg_Gamma_Event_Number, agefile, by = 'Subject')
Avg_Gamma_Event_Number$ageInverse <- 1/Avg_Gamma_Event_Number$age

ggplot(data = Avg_Gamma_Event_Number[Avg_Gamma_Event_Number$visitno <2,], aes(x = age, y = Gamma_Event_Number_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Gamma_Event_Number[Avg_Gamma_Event_Number$visitno <2,], Gamma_Event_Number_Mean ~ ageInverse)
anova(lm.model)


Avg_Gamma_Event_Number_SD <- aggregate(Gamma_Event_Number_SD ~ Subject, Everything, mean)
Avg_Gamma_Event_Number_SD <- merge(Avg_Gamma_Event_Number_SD, agefile, by = 'Subject')
Avg_Gamma_Event_Number_SD$ageInverse <- 1/Avg_Gamma_Event_Number_SD$age

ggplot(data = Avg_Gamma_Event_Number_SD[Avg_Gamma_Event_Number_SD$visitno <2,], aes(x = age, y = Gamma_Event_Number_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Gamma_Event_Number_SD[Avg_Gamma_Event_Number_SD$visitno <2,], Gamma_Event_Number_SD ~ ageInverse)
anova(lm.model)


#event duration 
Avg_Gamma_Event_Duration <- aggregate(Gamma_Event_Duration_Mean ~ Subject, Everything, mean)
Avg_Gamma_Event_Duration <- merge(Avg_Gamma_Event_Duration, agefile, by = 'Subject')
Avg_Gamma_Event_Duration$ageInverse <- 1/Avg_Gamma_Event_Duration$age

ggplot(data = Avg_Gamma_Event_Duration[Avg_Gamma_Event_Duration$visitno <2,], aes(x = age, y = Gamma_Event_Duration_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Gamma_Event_Duration[Avg_Gamma_Event_Duration$visitno <2,], Gamma_Event_Duration_Mean ~ ageInverse)
anova(lm.model)


Avg_Gamma_Event_Duration_SD <- aggregate(Gamma_Event_Duration_SD ~ Subject, Everything, mean)
Avg_Gamma_Event_Duration_SD <- merge(Avg_Gamma_Event_Duration_SD, agefile, by = 'Subject')
Avg_Gamma_Event_Duration_SD$ageInverse <- 1/Avg_Gamma_Event_Duration_SD$age

ggplot(data = Avg_Gamma_Event_Duration_SD[Avg_Gamma_Event_Duration_SD$visitno <2,], aes(x = age, y = Gamma_Event_Duration_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Gamma_Event_Duration_SD[Avg_Gamma_Event_Duration_SD$visitno <2,], Gamma_Event_Duration_SD ~ ageInverse)
anova(lm.model)




## Look into theta values 
#trial power sd
Avg_Theta_trial_power_SD <- aggregate(log_Theta_Trial_Power_SD ~ Subject, Everything[Everything$log_Theta_Trial_Power_SD != '-Inf',], mean)
Avg_Theta_trial_power_SD <- merge(Avg_Theta_trial_power_SD, agefile, by = 'Subject')
Avg_Theta_trial_power_SD$ageInverse <- 1/Avg_Theta_trial_power_SD$age
ggplot(data = Avg_Theta_trial_power_SD[Avg_Theta_trial_power_SD$visitno <2,], aes(x = age, y = log_Theta_Trial_Power_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Theta_trial_power_SD[Avg_Theta_trial_power_SD$visitno < 2,], log_Theta_Trial_Power_SD ~ ageInverse)
anova(lm.model)

#trial power mean
Avg_Theta_trial_power_mean <- aggregate(log_Theta_Trial_Power_mean ~ Subject, Everything[Everything$log_Theta_Trial_Power_mean != '-Inf',], mean)
Avg_Theta_trial_power_mean <- merge(Avg_Theta_trial_power_mean, agefile, by = 'Subject')
Avg_Theta_trial_power_mean$ageInverse <- 1/Avg_Theta_trial_power_mean$age
ggplot(data = Avg_Theta_trial_power_mean[Avg_Theta_trial_power_mean$visitno < 2,], aes(x = age, y = log_Theta_Trial_Power_mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Theta_trial_power_mean[Avg_Theta_trial_power_mean$visitno < 2,], log_Theta_Trial_Power_mean ~ ageInverse)
anova(lm.model)

#event number 
Avg_Theta_Event_Number <- aggregate(Theta_Event_Number_Mean ~ Subject, Everything, mean)
Avg_Theta_Event_Number <- merge(Avg_Theta_Event_Number, agefile, by = 'Subject')
Avg_Theta_Event_Number$ageInverse <- 1/Avg_Theta_Event_Number$age

ggplot(data = Avg_Theta_Event_Number[Avg_Theta_Event_Number$visitno <2,], aes(x = age, y = Theta_Event_Number_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Theta_Event_Number[Avg_Theta_Event_Number$visitno <2,], Theta_Event_Number_Mean ~ ageInverse)
anova(lm.model)


Avg_Theta_Event_Number_SD <- aggregate(Theta_Event_Number_SD ~ Subject, Everything, mean)
Avg_Theta_Event_Number_SD <- merge(Avg_Theta_Event_Number_SD, agefile, by = 'Subject')
Avg_Theta_Event_Number_SD$ageInverse <- 1/Avg_Theta_Event_Number_SD$age

ggplot(data = Avg_Theta_Event_Number_SD[Avg_Theta_Event_Number_SD$visitno <2,], aes(x = age, y = Theta_Event_Number_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Theta_Event_Number_SD[Avg_Theta_Event_Number_SD$visitno <2,], Theta_Event_Number_SD ~ ageInverse)
anova(lm.model)


#event duration 
Avg_Theta_Event_Duration <- aggregate(Theta_Event_Duration_Mean ~ Subject, Everything, mean)
Avg_Theta_Event_Duration <- merge(Avg_Theta_Event_Duration, agefile, by = 'Subject')
Avg_Theta_Event_Duration$ageInverse <- 1/Avg_Theta_Event_Duration$age

ggplot(data = Avg_Theta_Event_Duration[Avg_Theta_Event_Duration$visitno <2,], aes(x = age, y = Theta_Event_Duration_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Theta_Event_Duration[Avg_Theta_Event_Duration$visitno <2,], Theta_Event_Duration_Mean ~ ageInverse)
anova(lm.model)


Avg_Theta_Event_Duration_SD <- aggregate(Theta_Event_Duration_SD ~ Subject, Everything, mean)
Avg_Theta_Event_Duration_SD <- merge(Avg_Theta_Event_Duration_SD, agefile, by = 'Subject')
Avg_Theta_Event_Duration_SD$ageInverse <- 1/Avg_Theta_Event_Duration_SD$age

ggplot(data = Avg_Theta_Event_Duration_SD[Avg_Theta_Event_Duration_SD$visitno <2,], aes(x = age, y = Theta_Event_Duration_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Theta_Event_Duration_SD[Avg_Theta_Event_Duration_SD$visitno <2,], Theta_Event_Duration_SD ~ ageInverse)
anova(lm.model)



## Look into Beta values 
#trial power sd
Avg_Beta_trial_power_SD <- aggregate(log_Beta_Trial_Power_SD ~ Subject, Everything[Everything$log_Beta_Trial_Power_SD != '-Inf',], mean)
Avg_Beta_trial_power_SD <- merge(Avg_Beta_trial_power_SD, agefile, by = 'Subject')
Avg_Beta_trial_power_SD$ageInverse <- 1/Avg_Beta_trial_power_SD$age
ggplot(data = Avg_Beta_trial_power_SD[Avg_Beta_trial_power_SD$visitno <2,], aes(x = age, y = log_Beta_Trial_Power_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Beta_trial_power_SD[Avg_Beta_trial_power_SD$visitno < 2,], log_Beta_Trial_Power_SD ~ ageInverse)
anova(lm.model)

#trial power mean
Avg_Beta_trial_power_mean <- aggregate(log_Beta_Trial_Power_mean ~ Subject, Everything[Everything$log_Beta_Trial_Power_mean != '-Inf',], mean)
Avg_Beta_trial_power_mean <- merge(Avg_Beta_trial_power_mean, agefile, by = 'Subject')
Avg_Beta_trial_power_mean$ageInverse <- 1/Avg_Beta_trial_power_mean$age
ggplot(data = Avg_Beta_trial_power_mean[Avg_Beta_trial_power_mean$visitno < 2,], aes(x = age, y = log_Beta_Trial_Power_mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Beta_trial_power_mean[Avg_Beta_trial_power_mean$visitno < 2,], log_Beta_Trial_Power_mean ~ ageInverse)
anova(lm.model)

#event number 
Avg_Beta_Event_Number <- aggregate(Beta_Event_Number_Mean ~ Subject, Everything, mean)
Avg_Beta_Event_Number <- merge(Avg_Beta_Event_Number, agefile, by = 'Subject')
Avg_Beta_Event_Number$ageInverse <- 1/Avg_Beta_Event_Number$age

ggplot(data = Avg_Beta_Event_Number[Avg_Beta_Event_Number$visitno < 2,], aes(x = age, y = Beta_Event_Number_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Beta_Event_Number[Avg_Beta_Event_Number$visitno <2,], Beta_Event_Number_Mean  ~ ageInverse)
anova(lm.model)


Avg_Beta_Event_Number_SD <- aggregate(Beta_Event_Number_SD ~ Subject, Everything, mean)
Avg_Beta_Event_Number_SD <- merge(Avg_Beta_Event_Number_SD, agefile, by = 'Subject')
Avg_Beta_Event_Number_SD$ageInverse <- 1/Avg_Beta_Event_Number_SD$age

ggplot(data = Avg_Beta_Event_Number_SD[Avg_Beta_Event_Number_SD$visitno <2,], aes(x = age, y = Beta_Event_Number_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Beta_Event_Number_SD[Avg_Beta_Event_Number_SD$visitno <2,], Beta_Event_Number_SD ~ ageInverse)
anova(lm.model)


#event duration 
Avg_Beta_Event_Duration <- aggregate(Beta_Event_Duration_Mean ~ Subject, Everything, mean)
Avg_Beta_Event_Duration <- merge(Avg_Beta_Event_Duration, agefile, by = 'Subject')
Avg_Beta_Event_Duration$ageInverse <- 1/Avg_Beta_Event_Duration$age

ggplot(data = Avg_Beta_Event_Duration[Avg_Beta_Event_Duration$visitno <2,], aes(x = age, y = Beta_Event_Duration_Mean )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Beta_Event_Duration[Avg_Beta_Event_Duration$visitno <2,], Beta_Event_Duration_Mean ~ ageInverse)
anova(lm.model)


Avg_Beta_Event_Duration_SD <- aggregate(Beta_Event_Duration_SD ~ Subject, Everything, mean)
Avg_Beta_Event_Duration_SD <- merge(Avg_Beta_Event_Duration_SD, agefile, by = 'Subject')
Avg_Beta_Event_Duration_SD$ageInverse <- 1/Avg_Beta_Event_Duration_SD$age

ggplot(data = Avg_Beta_Event_Duration_SD[Avg_Beta_Event_Duration_SD$visitno <2,], aes(x = age, y = Beta_Event_Duration_SD )) + geom_point() + stat_smooth(method='lm',  formula='y~I(1/x)')

lm.model <- lm(data=Avg_Beta_Event_Duration_SD[Avg_Beta_Event_Duration_SD$visitno <2,], Beta_Event_Duration_SD ~ ageInverse)
anova(lm.model)





















