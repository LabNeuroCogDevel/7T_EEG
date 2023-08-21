if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2,lme4, multcomp, tidyr, tidyverse,wesanderson, reshape2, corrplot, ggExtra, cowplot, ggsci, mgcv, ggpointdensity, devtools, ggpubr, itsadug, useful, rmarkdown, tinytex, gratia, gamm4)
devtools::install_github('LabneuroCogDevel/LNCDR')
library(LNCDR)

## gamma avg power
gammaPCA <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_AvgPower_PCA_scores.csv')
gammaPCA$inverseAge <- 1/gammaPCA$age

ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_1)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_2)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_3)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_4)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_5)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_6)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_7)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_8)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_9)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gammaPCA[gammaPCA$visitno < 2,], aes(x=age, y=scores_10)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for component 1
lm.model <- lm(data=gammaPCA[gammaPCA$visitno < 2,], scores_1 ~ inverseAge )
anova(lm.model)

#stats for component 3
lm.model <- lm(data=gammaPCA[gammaPCA$visitno < 2,], scores_3 ~ inverseAge )
anova(lm.model)


## gamma peak power
gamma_PP_PCA <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_PeakPower_PCA_scores.csv')
gamma_PP_PCA$inverseAge <- 1/gamma_PP_PCA$age

ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_1)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_2)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_3)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_4)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_5)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_6)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_7)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_8)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_9)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], aes(x=age, y=scores_10)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for component 1
lm.model <- lm(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], scores_1 ~ inverseAge )
anova(lm.model)

#stats for component 5
lm.model <- lm(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2,], scores_3 ~ inverseAge )
anova(lm.model)

#stats for component 8
lm.model <- lm(data=gamma_PP_PCA[gamma_PP_PCA$visitno < 2 ,], scores_8 ~ inverseAge )
anova(lm.model)



## Theta avg power
Theta_AP_PCA <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_AP_PCA_scores.csv')
Theta_AP_PCA$inverseAge <- 1/Theta_AP_PCA$age

ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_1)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_2)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_3)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_4)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_5)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_6)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_7)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_8)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_9)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], aes(x=age, y=scores_10)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for component 2
lm.model <- lm(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], scores_1 ~ inverseAge )
anova(lm.model)

#stats for component 3
lm.model <- lm(data=Theta_AP_PCA[Theta_AP_PCA$visitno < 2,], scores_3 ~ inverseAge )
anova(lm.model)



## Theta peak power
Theta_PP_PCA <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_PP_PCA_scores.csv')
Theta_PP_PCA$inverseAge <- 1/Theta_AP_PCA$age

ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_1)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_2)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_3)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_4)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_5)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_6)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_7)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_8)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_9)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], aes(x=age, y=scores_10)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for component 2
lm.model <- lm(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], scores_2 ~ inverseAge )
anova(lm.model)

#stats for component 3
lm.model <- lm(data=Theta_PP_PCA[Theta_PP_PCA$visitno < 2,], scores_3 ~ inverseAge )
anova(lm.model)



## Beta avg power
Beta_AP_PCA <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_AP_PCA_scores.csv')
Beta_AP_PCA$inverseAge <- 1/Beta_AP_PCA$age

ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_1)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_2)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_3)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_4)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_5)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_6)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_7)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_8)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_9)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], aes(x=age, y=scores_10)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for component 2
lm.model <- lm(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], scores_4 ~ inverseAge )
anova(lm.model)

#stats for component 3
lm.model <- lm(data=Beta_AP_PCA[Beta_AP_PCA$visitno < 2,], scores_3 ~ inverseAge )
anova(lm.model)



## Beta peak power
Beta_PP_PCA <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_PP_PCA_scores.csv')
Beta_PP_PCA$inverseAge <- 1/Beta_PP_PCA$age

ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_1)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_2)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_3)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_4)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_5)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_6)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_7)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_8)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_9)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
ggplot(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], aes(x=age, y=scores_10)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')


#stats for component 2
lm.model <- lm(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], scores_1 ~ inverseAge )
anova(lm.model)

#stats for component 3
lm.model <- lm(data=Beta_PP_PCA[Beta_PP_PCA$visitno < 2,], scores_3 ~ inverseAge )
anova(lm.model)
