#!/usr/bin/env Rscript

BrainBehaviorStats <- function () {

# Libraries ----

library(LNCDR)
library(data.table)
library(dplyr)
library(factoextra)
library(ggplot2)
library(e1071)
library(caret)
attach(mtcars)
library(grid)
library(gridExtra)
library(plotrix)
library(mgcv)
library(readxl)
library(lme4)
library(lubridate)
library(checkmate)
library(lmerTest)
library(tidyr)
library(jtools)

# Brain vs Beh ----
fooofMRSbehavior <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/allSubjectsDLPFCfooofMRSBehaviorMeasures_20230822.csv")

### Accuracy vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior, 
               aes(x = absBestError, y = Exponent, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Exponent ~  absBestError+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior)
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(Exponent ~  absBestError* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Exponent ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Accuracy vs Offset ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Offset, x = absBestError, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Offset~  absBestError+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Offset~  absBestError* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Offset ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy vs Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = absBestError, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj ~ absBestError + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Glu_gamadj ~ absBestError* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy vs GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = absBestError, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab(" GABA ")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj ~ absBestError + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GABA_gamadj ~ absBestError* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





### Accuracy vs Imbalance ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GluGABAimbalanceABS, x = absBestError, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GluGABAimbalanceABS ~ absBestError + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GluGABAimbalanceABS ~ absBestError* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


###  Accuracy vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior, 
               aes(y = Ratio_gamadj, x = absBestError, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Ratio_gamadj~  absBestError+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  absBestError* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))




### Accuracy Var vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Exponent, x = absBestError_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Exponent ~  absBestError_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Exponent~ absBestError_sd * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Exponent ~ absBestError_sd+ age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ absBestError_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Accuracy Var vs Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = absBestError_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj ~ absBestError_sd + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Glu_gamadj ~ absBestError_sd* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy Var vs GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = absBestError_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab(" GABA ")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj ~ absBestError_sd + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GABA_gamadj ~ absBestError_sd* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





### Accuracy Var vs Offset ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Offset, x = absBestError_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Offset ~  absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Offset~  absBestError_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Offset ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Accuracy Var  vs Imbalance ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GluGABAimbalanceABS, x = absBestError_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab("GluGABAimbalanceABS")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GluGABAimbalanceABS ~ absBestError_sd+ inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GluGABAimbalanceABS ~ absBestError_sd* inverseAge + Region +GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy Var vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = absBestError_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Ratio_gamadj~  absBestError_sd+ ageage + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  absBestError_sd* age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Latency vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior , 
               aes(y = Exponent, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Exponent ~  mgsLatency+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Exponent~  mgsLatency* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Exponent ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


###  Latency  vs Offset ----

lunaize(ggplot(data = fooofMRSbehavior, 
               aes(y = Offset, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Offset~  mgsLatency+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior)
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(Offset ~ mgsLatency * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior)
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Offset ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Latency VS Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj~  mgsLatency+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Glu_gamadj~  mgsLatency* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


# visualizing age interactions 
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ 
  theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))

###  Latency VS GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("GABA")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj~  mgsLatency+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GABA_gamadj~  mgsLatency* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


# visualizing age interactions 
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ 
  theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))



###  Latency vs Imbalance ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GluGABAimbalanceABS, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("Glu GABA imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GluGABAimbalanceABS ~  mgsLatency+ inverseAge + Region + GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GluGABAimbalanceABS ~ mgsLatency* inverseAge  + Region + GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

# visualizing age interactions
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,18,Inf), labels=c("adol","adults"))), 
               aes(y = GluGABAimbalanceABS, x = mgsLatency, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30)) + scale_color_manual(values=c("gold3", "blue4", "red4"))


lm.model <- lmer(GluGABAimbalanceABS ~  mgsLatency + Region + (1|luna), data = fooofMRSbehavior%>% mutate(ageGroup = cut(age, breaks=c(0,18,Inf), labels=c("adol","adults")))  %>% filter(ageGroup == "adults"))
car::Anova(lm.model)
summ(lm.model)



###  Latency vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Ratio_gamadj~  mgsLatency+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  mgsLatency* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


# visualizing age interactions 
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ 
  theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))





### Latency Var  vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior, 
               aes(y = Exponent, x = mgsLatency_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group =1), method="lm", alpha = 0.8)) + 
  scale_color_manual(values=c("gold3", "blue4")) +
  xlab("Latency Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Exponent~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Exponent~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Exponent ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





###  Latency Var vs Offset ----
lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Offset, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group =1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency Var") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Offset~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(Offset ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Offset ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Latency Var vs Glu ----
lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group =1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency Var") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(Glu_gamadj ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





###  Latency Var vs GABA ----
lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group =1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency Var") + ylab("GABA")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(GABA_gamadj ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





### Latency Var vs Imbalance ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GluGABAimbalanceABS, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency Var") + ylab("Glu GABA imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(GluGABAimbalanceABS ~  mgsLatency_sd+ inverseAge  + Region + GMrat +(1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GluGABAimbalanceABS ~ mgsLatency_sd* inverseAge  + Region + GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Latency Var vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) +
  xlab("Latency Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Ratio_gamadj ~  mgsLatency_sd+ age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  mgsLatency_sd* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)


AIC(lmer(Ratio_gamadj ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Spatial Span Max vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Exponent, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmerTest::lmer( Exponent ~  SSP_maxSpan + inverseAge  + Region + Condition + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)

lm.model <- lmer( Exponent ~  SSP_maxSpan*age  + Region + Condition + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Exponent ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))




### Spatial Span Max vs Offset ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Offset, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmerTest::lmer( Offset ~  SSP_maxSpan + inverseAge  + Region + Condition + (1|luna), data = fooofMRSbehavior)
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmerTest::lmer( Offset ~  SSP_maxSpan*inverseAge  + Region + Condition + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Offset ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer( Glu_gamadj ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior %>% filter(ageGroup == "Adult") )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( Glu_gamadj ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("GABA")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer( GABA_gamadj ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior  %>% filter(ageGroup == "Adol"))
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( GABA_gamadj ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer( Ratio_gamadj ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior  %>% filter(ageGroup == "Adol"))
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( Ratio_gamadj ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs Imbalance ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GluGABAimbalanceABS, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GluGABAimbalanceABS ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior  %>% filter(ageGroup == "Adol"))
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( GluGABAimbalanceABS ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))



# Fooof vs Behavior ----

fooofVars <- c('Offset','Exponent')
behVars <- c('absBestError','absBestError_sd','mgsLatency', 'mgsLatency_sd', 'SSP_maxSpan')


output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    model <- paste0(fooofVar, ' ~ ', behVar, '+ inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))
    b <- model.out$coefficients[2,1]
    t <- model.out$coefficients[2,4]
    p <- model.out$coefficients[2,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(fooofVar, behVar, b, t, pcor, p))
  }
}  
output  

# Looking at age interaction
output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    model <- paste0(fooofVar, ' ~ ', behVar, '* inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))
    b <- model.out$coefficients[6,1]
    t <- model.out$coefficients[6,4]
    p <- model.out$coefficients[6,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(fooofVar, behVar, b, t, pcor, p))
  }
}  
output  

# looks at GAMs
output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    
    model <- paste0(fooofVar, ' ~ ', behVar, '+ s(age, k = 3) + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    b <- summary(model.out$gam)$p.table[2,1]
    t <- summary(model.out$gam)$p.table[2,3]
    p <- summary(model.out$gam)$p.table[2,4]
    regionb <- summary(model.out$gam)$p.table[3,1]
    regiont <- summary(model.out$gam)$p.table[3,3]
    regionp <- summary(model.out$gam)$p.table[3,4]
    condb <- summary(model.out$gam)$p.table[4,1]
    condt <- summary(model.out$gam)$p.table[4,3]
    condp <- summary(model.out$gam)$p.table[4,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    regionpcor <- p.adjust((regionp), method = "bonferroni", n = 4)
    condpcor <- p.adjust((condp), method = "bonferroni", n = 4)
    
    
    
    output <- rbind(output, data.frame(fooofVar, behVar, b, t, p, pcor, regionb, regiont, regionp, regionpcor, condb, condt, condp, condpcor))
  }
}
output  


# looks at age interactions
output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    
    model <- paste0(fooofVar, ' ~ ','s(age, k = 3, by = ', behVar,') + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    f <- summary(model.out$gam)$s.table[1,3]
    p <- summary(model.out$gam)$s.table[1,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(fooofVar, behVar, f, p, pcor))
  }
}
output  



# MRS vs Behavior ----
#'absBestError','absBestError_sd','mgsLatency', 'mgsLatency_sd', 'Glu_gamadj', 'GABA_gamadj', 'GluGABAimbalanceABS'
mrsiVars <- c('GluGABAimbalanceABS')
behVars <- c('mgsLatency')

output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    model <- paste0(mrsiVar, ' ~ ', behVar, '+ inverseAge + Region + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))  
    b <- model.out$coefficients[2,1]
    t <- model.out$coefficients[2,4]
    p <- model.out$coefficients[2,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, b, t, pcor, p))
  }
}  
output  

# Looking at age interaction
output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    model <- paste0(mrsiVar, ' ~ ', behVar, '* inverseAge + Region + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))
    b <- model.out$coefficients[5,1]
    t <- model.out$coefficients[5,4]
    p <- model.out$coefficients[5,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, b, t, pcor, p))
  }
}  
output  



# looks at GAMs
output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    
    model <- paste0(mrsiVar, ' ~ ', behVar, '+ s(age, k = 3) + Region ')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    b <- summary(model.out$gam)$p.table[2,1]
    t <- summary(model.out$gam)$p.table[2,3]
    p <- summary(model.out$gam)$p.table[2,4]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, b, t, p, pcor))
  }
}
output  

# looks at age interactions
output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    
    model <- paste0(mrsiVar, ' ~ ','s(age, k = 3, by = ', behVar,') + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    f <- summary(model.out$gam)$s.table[1,3]
    p <- summary(model.out$gam)$s.table[1,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, f, p, pcor))
  }
}
output  

}