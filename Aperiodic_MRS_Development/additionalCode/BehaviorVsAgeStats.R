#!/usr/bin/env Rscript
BehaviorVsAgeStats <- function (){

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

## Load Dataframes ----
behav <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/allSubjectsBehavior.csv")

## Behavior vs age ----

#### Best Saccade ----
lunaize(ggplot(data = behav, 
               aes(y = absBestError, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='black', method="lm", formula = 'y~I(1/x)', alpha = 0.4, size = 1) + 
          scale_color_manual(values=c("gold3", "blue4")))

gam.model <- gam(absBestError ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, agevar = 'age', yvar = 'absBestError', idvar = "luna" ,draw_points = T, xplotname = "Age", yplotname = "Mean MGS Accuracy (degs)"))


gam.model <-  gamm(absBestError ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)


#### Best Saccade Var ----
gam.model <- gam(absBestError_sd ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, agevar = 'age', idvar = "luna", yvar = 'absBestError_sd', draw_points = T, xplotname = "Age", yplotname = "SD MGS Accuracy (degs)"))

gam.model <-  gamm(absBestError_sd ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

#### MGS Latency ----
gam.model <- gam(mgsLatency ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'mgsLatency', draw_points = T, xplotname = "Age", yplotname = "Mean MGS Latency (s)"))

gam.model <-  gamm(mgsLatency ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

#### MGS Latency Var ----
gam.model <- gam(mgsLatency_sd ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'mgsLatency_sd', draw_points = T, xplotname = "Age", yplotname = "SD MGS Latency (s)"))

gam.model <-  gamm(mgsLatency_sd ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

#### Spatial Span Max ----

lunaize(ggplot(data = behav, 
               aes(y = SSP_maxSpan, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='blue', method="lm", formula = 'y~I(1/x)', alpha = 0.4))

gam.model <-  gamm(SSP_maxSpan ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

lm.model <- lmerTest::lmer(SSP_maxSpan ~ age + (1|luna), data = behav)
summary(lm.model)


gam.model <- gam(SSP_maxSpan ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'SSP_maxSpan', draw_points = T, xplotname = "Age", yplotname = "Max Sequence Length"))

gam.model <-  gamm(SSP_maxSpan ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

AIC(lmer(Exponent ~ age + Condition + Region + (1|luna), data = fooofLong), 
    lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooofLong))


#### Spatial Span Errors ----

lunaize(ggplot(data = behav, 
               aes(y = SSP_nErrors, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='blue', method="lm", formula = 'y~I(1/x)', alpha = 0.4))

gam.model <-  gamm(SSP_nErrors ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

lm.model <- lmerTest::lmer(SSP_nErrors ~ age + (1|luna), data = behav)
summary(lm.model)


gam.model <- gam(SSP_nErrors ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'SSP_nErrors', draw_points = T, xplotname = "Age", yplotname = "SSP nErrors"))

#### Spatial Span nTrials ----

lunaize(ggplot(data = behav, 
               aes(y = SSP_nTrials, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='blue', method="lm", formula = 'y~I(1/x)', alpha = 0.4))

gam.model <-  gamm(SSP_nTrials ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

lm.model <- lmerTest::lmer(SSP_nTrials ~ age + (1|luna), data = behav)
summary(lm.model)


gam.model <- gam(SSP_nTrials ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'SSP_nTrials', draw_points = T, xplotname = "Age", yplotname = "SSP nTrials"))


}