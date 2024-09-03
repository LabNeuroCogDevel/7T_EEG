
# Load Libraries 

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
library(eegUtils)
library(tvem)
library(interactions)
library(akima)


outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


# Load datasets ----

maxEntropy <- read.csv( '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_maxEntropy_allChansAvg.csv')
AvgChans_MSE <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_multiscaleEntropy_allChansAvg.csv')%>% separate(Subject, c('lunaID','vdate'))


AvgChans_MSE_outlier <- AvgChans_MSE %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier))


merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
behav <- merge7t[c("lunaid","visitno","eeg.age", "eeg.vgsLatency_DelayAll","eeg.BestError_DelayAll", "eeg.BestError_sd_DelayAll",
                   "eeg.mgsLatency_DelayAll", "eeg.mgsLatency_sd_DelayAll","cantab.ssp_max.span","cantab.ssp_nerrors", "cantab.ssp_ntrials", 
                   "antiET.Cor", "antiET.Err", "antiET.ErrCor", "antiET.cor.lat")]


behav <- behav %>% 
  mutate(antiET.Cor = ifelse(antiET.Cor < 16, NA, antiET.Cor))

behav$corrat <- behav$antiET.Cor / (behav$antiET.Cor + behav$antiET.Err + behav$antiET.ErrCor)

names(behav)[names(behav) == 'lunaid'] <- 'lunaID'
names(behav)[names(behav) == 'eeg.age'] <- 'age'


# Merge datasets ----

entropyBehav <- merge(maxEntropy %>% select(lunaID, visitno, maxEntropy), AvgChans_MSE_outlier  %>% 
                        select(lunaID, visitno, Var1, age), by = c('lunaID', 'visitno')) %>% 
  merge(., behav, by = c('lunaID', 'visitno', 'age'))%>%
  mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

names(entropyBehav)[names(entropyBehav) == 'Var1'] <- 'AOC'


# Age vs Anti Correct Response Rate ----

lunaize(ggplot(data = entropyBehav, aes(x = age, y = corrat)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

#Just Adults ----
AdultsEntropyBehav <- entropyBehav %>% filter(age >=18)
# AOC ----
##Anti Correct Response Rate ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue")) + 
  xlab("AUC") + ylab("Anti Saccade Correct Rate")
        

model <- lmer(corrat ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") ) + 
  xlab("AUC") + ylab("Anti Saccade Correct Rate")


model <- lmer(corrat ~ AOC + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## Anti Correct Response Rate Latency ----
### Quadratic ----

lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("Anti Saccade Correct Rate Latency")


model <- lmer(antiET.cor.lat ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----

lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("Anti Saccade Correct Rate Latency")


model <- lmer(antiET.cor.lat ~ AOC + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## MGS Best Error ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm" ,se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ AOC + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

## MGS Best Error Var ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ AOC + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

## MGS Latency ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("MGS Latency")


model <- lmer(eeg.mgsLatency_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("MGS Latency")


model <- lmer(eeg.mgsLatency_DelayAll ~ AOC + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## MGS Latency Var ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ AOC + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## Spatial Span ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = AOC, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ AOC + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

# Max Entropy ----
##Anti Correct Response Rate ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate")


model <- lmer(corrat ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate")


model <- lmer(corrat ~ maxEntropy + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## Anti Correct Response Rate Latency ----
### Quadratic ----

lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate Lat")


model <- lmer(antiET.cor.lat ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----

lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate Lat")



model <- lmer(antiET.cor.lat ~ maxEntropy + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## MGS Best Error ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm" ,se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ maxEntropy + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

## MGS Best Error Var ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ maxEntropy + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

## MGS Latency ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency")


model <- lmer(eeg.mgsLatency_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency")


model <- lmer(eeg.mgsLatency_DelayAll ~ maxEntropy + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## MGS Latency Var ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ maxEntropy + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)


## Spatial Span ----
### Quadratic ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = AdultsEntropyBehav, aes(x = maxEntropy, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ maxEntropy + age + (1|lunaID), data = AdultsEntropyBehav)
summary(model)

# All Subjects ----

# AOC ----
##Anti Correct Response Rate ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue")) + 
  xlab("AUC") + ylab("Anti Saccade Correct Rate")


model <- lmer(corrat ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") ) + 
  xlab("AUC") + ylab("Anti Saccade Correct Rate")


model <- lmer(corrat ~ AOC + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(corrat ~ AOC*age + (1|lunaID), data = entropyBehav)
summary(model)


## Anti Correct Response Rate Latency ----
### Quadratic ----

lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("Anti Saccade Correct Rate Latency")


model <- lmer(antiET.cor.lat ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----

lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("Anti Saccade Correct Rate Latency")


model <- lmer(antiET.cor.lat ~ AOC + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(antiET.cor.lat ~ AOC*age + (1|lunaID), data = entropyBehav)
summary(model)


## MGS Best Error ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav%>%
                 mutate(ageGroup = cut(age, c(0,18.24,Inf), labels = c('Adol','Adults'))), 
               aes(x = AOC, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(aes(group = ageGroup, color = ageGroup), method = "lm" ,se = FALSE) )+ 
  xlab("AUC") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ AOC + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.BestError_DelayAll ~ AOC*age + (1|lunaID), data = entropyBehav)
summary(model)

## try johnson neyman plots
johnson_neyman(model, pred = AOC, modx = age, plot = TRUE, title = "MGS Accuracy ~ AOC*age")


## MGS Best Error Var ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### linear ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+ 
  xlab("AUC") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ AOC + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.BestError_sd_DelayAll ~ AOC*age + (1|lunaID), data = entropyBehav)
summary(model)

## MGS Latency ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("MGS Latency")


model <- lmer(eeg.mgsLatency_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----


lunaize(ggplot(data = entropyBehav%>%
                 mutate(ageGroup = cut(age, c(0,28.9,Inf), labels = c('10-28','28+'))), 
               aes(x = AOC, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(aes(group = ageGroup, color = ageGroup), method = "lm" ,se = FALSE) )+ 
  xlab("AUC") + ylab("MGS Latency")



model <- lmer(eeg.mgsLatency_DelayAll ~ AOC + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.mgsLatency_DelayAll ~ AOC*age + (1|lunaID), data = entropyBehav)
summary(model)

## try johnson neyman plots
johnson_neyman(model, pred = AOC, modx = age, plot = TRUE, title = "MGS Latency ~ AOC*age")


## MGS Latency Var ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ AOC + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.mgsLatency_sd_DelayAll ~ AOC*age + (1|lunaID), data = entropyBehav)
summary(model)


## Spatial Span ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ poly(AOC, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----


lunaize(ggplot(data = entropyBehav%>%
                 mutate(ageGroup = cut(age, c(0,16.29,Inf), labels = c('10-16.29','16.29+'))), 
               aes(x = AOC, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(aes(group = ageGroup, color = ageGroup), method = "lm" ,se = FALSE) )+ 
  xlab("AUC") + ylab("Spatial Span")


lunaize(ggplot(data = entropyBehav, aes(x = AOC, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("AUC") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ AOC + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(cantab.ssp_max.span ~ AOC*age + (1|lunaID), data = entropyBehav)
summary(model)

## try johnson neyman plots
johnson_neyman(model, pred = AOC, modx = age, plot = TRUE, title = "Spatial Span ~ AOC*age")



# Max Entropy ----
##Anti Correct Response Rate ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate")


model <- lmer(corrat ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = corrat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate")


model <- lmer(corrat ~ maxEntropy + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(corrat ~ maxEntropy*age + (1|lunaID), data = entropyBehav)
summary(model)

## Anti Correct Response Rate Latency ----
### Quadratic ----

lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate Lat")


model <- lmer(antiET.cor.lat ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----

lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = antiET.cor.lat)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Anti Saccade Correct Response Rate Lat")



model <- lmer(antiET.cor.lat ~ maxEntropy + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(antiET.cor.lat ~ maxEntropy*age + (1|lunaID), data = entropyBehav)
summary(model)


## MGS Best Error ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.BestError_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm" ,se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy")

model <- lmer(eeg.BestError_DelayAll ~ maxEntropy + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.BestError_DelayAll ~ maxEntropy*age + (1|lunaID), data = entropyBehav)
summary(model)

## MGS Best Error Var ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### linear ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.BestError_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Accuracy Var")


model <- lmer(eeg.BestError_sd_DelayAll ~ maxEntropy + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.BestError_sd_DelayAll ~ maxEntropy*age + (1|lunaID), data = entropyBehav)
summary(model)

## MGS Latency ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency")


model <- lmer(eeg.mgsLatency_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency")


model <- lmer(eeg.mgsLatency_DelayAll ~ maxEntropy + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.mgsLatency_DelayAll ~ maxEntropy*age + (1|lunaID), data = entropyBehav)
summary(model)

## MGS Latency Var ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = eeg.mgsLatency_sd_DelayAll)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("MGS Latency Var")


model <- lmer(eeg.mgsLatency_sd_DelayAll ~ maxEntropy + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(eeg.mgsLatency_sd_DelayAll ~ maxEntropy*age + (1|lunaID), data = entropyBehav)
summary(model)

## Spatial Span ----
### Quadratic ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ poly(maxEntropy, 2, raw = TRUE) + age + (1|lunaID), data = entropyBehav)
summary(model)

### Linear ----
lunaize(ggplot(data = entropyBehav, aes(x = maxEntropy, y = cantab.ssp_max.span)) + geom_point(na.rm=T) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2, na.rm=T) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") )+  
  xlab("Max Entropy") + ylab("Spatial Span")


model <- lmer(cantab.ssp_max.span ~ maxEntropy + age + (1|lunaID), data = entropyBehav)
summary(model)

model <- lmer(cantab.ssp_max.span ~ maxEntropy*age + (1|lunaID), data = entropyBehav)
summary(model)
