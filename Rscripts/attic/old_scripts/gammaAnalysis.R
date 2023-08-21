if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2,lme4, multcomp, tidyr, tidyverse,wesanderson, reshape2, corrplot, ggExtra, cowplot, ggsci, mgcv, ggpointdensity, devtools, ggpubr, itsadug, useful, rmarkdown, tinytex, gratia, gamm4)
devtools::install_github('LabneuroCogDevel/LnCDR')
library(LNCDR)


power <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/power_table.csv')
power$logPower <- log(power$Power)
power$epoch <- as.factor(power$epoch)
ggplot(data=power, aes(x=Age, y=logPower, color=epoch)) + geom_point() + stat_smooth(method='loess')

e24 <- power[power$epoch=='24',]
e46 <- power[power$epoch=='46',]
e <- merge(e24, e46, by='Age')
ggplot(data=e, aes(x=logPower.x, y=logPower.y)) + geom_point() + stat_smooth(method='lm') + geom_abline(intercept=0, slope=1)
power$Subjectf=as.factor(power$Subject)
m <- gam(logPower ~ s(Age, k=4)+s(Subjectf,bs="re"),
         data = power, REML=TRUE)

#m <- gam(logPower ~ s(Age, k=4),
         #data = power, REML=TRUE)
plot(m)
summary(m)


poweraverage<-power %>% group_by(Subjectf,Age) %>% summarise(meanpower=mean(logPower))
ggplot(data=poweraverage, aes(x=Age, y=meanpower)) + geom_point() + stat_smooth(method='gam') #+ geom_abline(intercept=0, slope=1)

#meanmodel<-gam(meanpower ~ s(Age, k=4),
               #data = poweraverage, REML=TRUE)


MRSI <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/subj_label_val_gm_24specs_20191102.csv')
MRSI <- MRSI[MRSI$roi_label == 'L DLPFC',]
MRSI <- MRSI %>% separate(ld8,c("luna","vdate"),sep="_")

all_data <- merge(e, MRSI, by.x = 'Subject.x', by.y = 'luna',  all = T)



ggplot(data=all_data[all_data$GABA.Cre>0 & all_data$GABA.Cre<2.5 & all_data$GABA..SD < 20,], aes(x=Age, y=GABA.Cre, color=epoch.x )) + geom_point() + stat_smooth(method='loess')
ggplot(data=all_data[all_data$Glu.Cre>0 & all_data$Glu.Cre<2.5 & all_data$Glu..SD < 20,], aes(x=Age, y=Glu.Cre, color=epoch.x)) + geom_point() + stat_smooth(method='lm')


ggplot(data=all_data[all_data$GABA.Cre>0 & all_data$GABA.Cre<2.5 & all_data$GABA..SD < 20,], aes(x=logPower.x, y=GABA.Cre, color=epoch.x)) + geom_point() + stat_smooth(method='loess')
ggplot(data=all_data[all_data$Glu.Cre>0 & all_data$Glu.Cre<2.5 & all_data$Glu..SD < 20,], aes(x=logPower.x, y=Glu.Cre, color=Age)) + geom_point() + stat_smooth(method='lm')
ggplot(data=all_data[all_data$GABA.Cre>0  & all_data$GABA.Cre<1 & all_data$GABA..SD < 20,], aes(x=logPower.x, y=GABA.Cre, color=Age)) + geom_point() + stat_smooth(method='lm')





data=all_data[all_data$Glu..SD > 20,]

all_data$ratio <- all_data$GABA.Cre/all_data$Glu.Cre
ggplot(data=all_data[all_data$ratio > 0 & all_data$ratio <.75 & all_data$GABA..SD < 20,], aes(x=logPower.x, y=ratio, color=Age)) + geom_point() + stat_smooth(method='lm')


dontHave <- all_data %>% filter(is.na(Age))
write.csv(data, file="/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/GluOver20.csv")

allSubs <- all_data$Subject.x
write.csv(allSubs, file="/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/AllSubs.csv")





