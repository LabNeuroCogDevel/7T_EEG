#!/usr/bin/env Rscript

MSRIvsFOOOFstats <- function () {
  
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

# FOOOF vs MRS ----
fooofMRS <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectsDLPFCfooofMRSMeasures_20230822.csv")

## GLU VS Exponent ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Exponent, y = Glu_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) +
          geom_point(aes(shape=Region),alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glutamate") + xlab("Exponent") + theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Glu_gamadj ~ Exponent + age + Condition + Region + (1|luna), data = fooofMRS )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(Glu_gamadj ~ Exponent * age + Condition + Region + (1|luna), data = fooofMRS )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)


summary(lmerTest::lmer(Exponent ~ Glu_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Glu_gamadj + inverseAge + Region + Condition + (1|luna), data = fooofMRS))

summary(lmerTest::lmer(Exponent ~ Glu_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Glu_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))


# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ Glu_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS, random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)


AIC(lmer(Exponent ~ Glu_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ Glu_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))



## GLU VS Offset ----
lunaize(ggplot(data = fooofMRS, 
               aes(x = Offset, y = Glu_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(aes(shape =Region),alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glutamte") + xlab("Offset") + theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj ~ Offset + age + Condition + Region + (1|luna), data = fooofMRS )
car::Anova(lm.model)
summary(lm.model)


lm.model <- lmer(Glu_gamadj ~  Offset * age +Condition +Region +   (1|luna), data = fooofMRS )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)


summary(lmerTest::lmer(Offset ~ Glu_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
summary(lmerTest::lmer(Offset ~ Glu_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))

summary(lmerTest::lmer(Offset ~ Glu_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
summary(lmerTest::lmer(Offset ~ Glu_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))


# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Glu_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Glu_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Offset ~ Glu_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ Glu_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))




## GABA VS Exponent ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Exponent, y = GABA_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("GABA") + xlab("Exponent") + theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer(GABA_gamadj ~ Exponent + age + Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)


lm.model <- lmer(GABA_gamadj ~ Exponent *  age +Condition + Region+ (1|luna), data = fooofMRS )
summ(lm.model)


summary(lmerTest::lmer(Exponent ~ GABA_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ GABA_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

summary(lmerTest::lmer(Exponent ~ GABA_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ GABA_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ GABA_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zExp < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Exponent ~ GABA_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ GABA_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))

## GABA VS Offset ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Offset, y = GABA_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("GABA") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer(GABA_gamadj  ~ Offset + age + Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)


lm.model <- lmer(GABA_gamadj  ~ Offset * age + Condition +Region +   (1|luna), data = fooofMRS )
car::Anova(lm.model)
summ(lm.model)


summary(lmerTest::lmer(Offset ~ GABA_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GABA_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

summary(lmerTest::lmer(Offset ~ GABA_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GABA_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ GABA_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Offset ~ GABA_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ GABA_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))

## Ratio VS Exponent ----

lunaize(ggplot(data = fooofMRS, 
               aes(x = Exponent, y = Ratio_gamadj, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.1), method="lm", alpha = 0.8)  + 
          scale_color_manual(values=c("gold3", "blue4"))) +
  ylab("Glu/GABA Ratio") + xlab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Ratio_gamadj  ~ Exponent + age + Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)
car::Anova(lm.model)


lm.model <- lmer(Ratio_gamadj  ~ Exponent *age +Condition +Region+   (1|luna), data = fooofMRS )
summ(lm.model)
car::Anova(lm.model)


summary(lmerTest::lmer(Exponent ~ Ratio_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Ratio_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

summary(lmerTest::lmer(Exponent ~ Ratio_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Ratio_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ Ratio_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zExp < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Exponent ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ Ratio_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))


## Ratio VS Offset ----

lunaize(ggplot(data = fooofMRS, 
               aes(x = Offset, y = Ratio_gamadj, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glu/GABA Ratio") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')



lm.model <- lmer(Ratio_gamadj  ~  Offset + age + Condition  + Region + (1|luna), data = fooofMRS)
summ(lm.model)
car::Anova(lm.model)


lm.model <- lmer(Ratio_gamadj  ~ Offset *sex +age +Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)
car::Anova(lm.model)


summary(lmerTest::lmer(Offset ~ Ratio_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ Ratio_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

summary(lmerTest::lmer(Offset ~ Ratio_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ Ratio_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Ratio_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)


# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Ratio_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Offset ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ Ratio_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))


## Gaba glu imbalance VS offset ----

lunaize(ggplot(data = fooofMRS, 
               aes(x = Offset, y = GluGABAimbalanceABS, by =luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
  scale_color_manual(values=c("gold3", "blue4")) + 
  ylab("Glu GABA Imbalance") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ GluGABAimbalanceABS + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)


AIC(lmer(Offset ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ GluGABAimbalanceABS + inverseAge + Condition + Region + (1|luna), data = fooofMRS))


## Gaba glu imbalance VS exponent ----

lunaize(ggplot(data = fooofMRS, 
               aes(y = GluGABAimbalanceABS, x = Exponent, by =luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.01), method="lm",formula = y ~ poly(x,2),alpha=.8,size=1)) + 
  scale_color_manual(values=c("gold3", "blue4")) + 
  ylab("Glu GABA Imbalance") + xlab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')




# lmer approach, different forms of age, w/ and w/out interactions
summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + age + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + inverseage + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))

summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS*age + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))
summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS*inverseage + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))

summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + inverseage + Condition + Region + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))


# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ GluGABAimbalanceABS + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zExp < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Exponent ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ GluGABAimbalanceABS + inverseAge + Condition + Region + (1|luna), data = fooofMRS))



# Mediation ----

## Glu gaba imbalance on Exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(GluGABAimbalanceABS, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on exponent (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on imbalnce  (a)
model.M <- lme4::lmer(GluGABAimbalanceABS ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS imbalance on exponent (b)
model.Y <- lme4::lmer(Exponent ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GluGABAimbalanceABS", boot = FALSE, sims = 1000)
(summary(results))


## Glu gaba imbalance on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(GluGABAimbalanceABS, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on imbalance  (a)
model.M <- lme4::lmer(GluGABAimbalanceABS ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS imbalance on offset (b)
model.Y <- lme4::lmer(Offset ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GluGABAimbalanceABS", boot = FALSE, sims = 1000)
(summary(results))


## Glu gaba ratio on exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(Ratio_gamadj, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Ratio_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Exponent ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Ratio_gamadj", boot = FALSE, sims = 1000)
(summary(results))


## Glu gaba ratio on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(Ratio_gamadj, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Ratio_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Offset ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Ratio_gamadj", boot = FALSE, sims = 1000)
(summary(results))


## Glu on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(Glu_gamadj, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on glu  (a)
model.M <- lme4::lmer(Glu_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS glu on offset (b)
model.Y <- lme4::lmer(Offset ~ Glu_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Glu_gamadj", boot = FALSE, sims = 1000)
(summary(results))

## Glu on Exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(Glu_gamadj, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Glu_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Exponent ~ Glu_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Glu_gamadj", boot = FALSE, sims = 1000)
(summary(results))


## GABA on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(GABA_gamadj, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(GABA_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Offset ~ GABA_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GABA_gamadj", boot = FALSE, sims = 1000)
(summary(results))

## GABA on Exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(GABA_gamadj, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(GABA_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Exponent ~ GABA_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GABA_gamadj", boot = FALSE, sims = 1000)
(summary(results))

## Glu on Latency ----

mediationMatrix <- fooofMRSbehavior %>% dplyr::select(Glu_gamadj, mgsLatency, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(mgsLatency ~ age  + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Glu_gamadj ~ age  + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(mgsLatency  ~ Glu_gamadj +Region + age + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Glu_gamadj", boot = FALSE, sims = 1000)
(summary(results))



# FOOOF vs MRS Loop ----
### Linear Models ----

mrsiVars <- c('Ratio_gamadj','GluGABAimbalanceABS','Glu_gamadj','GABA_gamadj')
fooofVars <- c('Exponent', 'Offset')

# main effect
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ', mrsiVar, '+ inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRS))  
    b <- model.out$coefficients[2,1]
    t <- model.out$coefficients[2,4]
    p <- model.out$coefficients[2,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    
    regb <- model.out$coefficients[4,1]
    regt <- model.out$coefficients[4,4]
    regp <- model.out$coefficients[4,5]
    regpcor <- p.adjust((regp), method = "bonferroni", n = 2)
    
    
    conb <- model.out$coefficients[5,1]
    cont <- model.out$coefficients[5,4]
    conp <- model.out$coefficients[5,5]
    conpcor <- p.adjust((conp), method = "bonferroni", n = 2)
    
    
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor,regb, regt, regp,regpcor, conb, cont, conp,conpcor))   
  }
}  
output  

# age interaction 
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ', mrsiVar, '* inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRS))  
    b <- model.out$coefficients[6,1]
    t <- model.out$coefficients[6,4]
    p <- model.out$coefficients[6,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor))
  }
}  
output 



# region interaction 
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ', mrsiVar, '* Region + inverseAge + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRS))  
    b <- model.out$coefficients[6,1]
    t <- model.out$coefficients[6,4]
    p <- model.out$coefficients[6,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor))
  }
}  
output 


### GAMs ----
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ',mrsiVar, '+ s(age, k = 3) + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), data = fooofMRS, random=list(luna=~1)))
    b <- summary(model.out$gam)$p.table[2,1]
    t <- summary(model.out$gam)$p.table[2,3]
    p <- summary(model.out$gam)$p.table[2,4]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    
    regb <- summary(model.out$gam)$p.table[3,1]
    regt <- summary(model.out$gam)$p.table[3,3]
    regp <- summary(model.out$gam)$p.table[3,4]
    
    condb <- summary(model.out$gam)$p.table[4,1]
    condt <- summary(model.out$gam)$p.table[4,3]
    condp <- summary(model.out$gam)$p.table[4,4]
    
    regpcor <- p.adjust((regp), method = "bonferroni", n = 2)
    condpcor <- p.adjust((condp), method = "bonferroni", n = 2)
    
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor,regb, regt, regp,regpcor, condb, condt, condp,condpcor))    }
}
output 

# look at interactions 
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    
    model_formula <- as.formula(paste0(fooofVar, ' ~ ', ' s(age, k = 3, fx = T) + s(age, by = ',mrsiVar, ',k = 4, fx = T) + Region + Condition')) 
    
    model.out <- (mgcv::gamm(model_formula,
                             random=list(luna=~1),
                             data = fooofMRS))
    
    F <- summary(model.out$gam)$s.table[2,3]
    p <- summary(model.out$gam)$s.table[2,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    
    output <- rbind(output, data.frame(fooofVar, mrsiVar, F, p, pcor))
  }
}
output  

}


