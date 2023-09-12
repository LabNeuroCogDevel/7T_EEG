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

# Load Dataframe ----
fooofLong <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/Results/allSubjectsDLPFCfooofMeasures_20230613.csv")

# FOOOF error measures ---- 
fooofErrors <- read.csv(hera('/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsErrorMeasures_20230516.csv'))

fooof_channels <- merge(fooofErrors, chanLocs, by = c("Channel"))

fooof_regional <- rbind(
  fooof_channels %>% filter(Channel %in% c('F3', 'F5', 'F7')) %>% group_by(Subject, Condition) %>% 
    summarise(Error = mean(Error, na.rm=T), Rsqaured = mean(R.Squared, na.rm=T)) %>% mutate(Region = "LDLPFC"),
  
  fooof_channels %>% filter(Channel %in% c('F4', 'F6', 'F8')) %>% group_by(Subject, Condition) %>% 
    summarise(Error = mean(Error, na.rm=T), Rsqaured = mean(R.Squared, na.rm=T)) %>% mutate(Region = "RDLPFC")
)

avgError <- aggregate(.~Condition, data = fooof_regional[2:4], mean)
# FOOOF stats ----

## Exponent vs Offset ----

lunaize(ggplot(data = fooofLong, 
               aes(y = Exponent, x = Offset, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='black', method="gam", alpha = 0.4, size = 1) + 
          scale_color_manual(values=c("gold3", "blue4")))+ theme(legend.position = "none")+ 
  theme(text = element_text(size = 30))

corDF <- fooofLong
cor.test(corDF$Exponent , corDF$Offset, method = "pearson")

corDF <- fooofLong %>% filter(Condition == "eyesClosed")
cor.test(corDF$Exponent , corDF$Offset, method = "pearson")


lm.model <- lmerTest::lmer(Exponent ~ Offset + age + Condition + Region + (1|luna), data = fooofLong)
summary(lm.model)

lm.model.int <- lmerTest::lmer(Exponent ~ Offset*age + Condition + Region + (1|luna), data = fooofLong )
summary(lm.model.int)

AIC(lm.model, lm.model.int)


## Exponent vs age ----

lunaize(ggplot(data = fooofLong %>% filter(Condition == 'eyesOpen'), 
               aes(x = age, y = Exponent)) + 
          geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Exponent") + theme(legend.position = "none")

gam.model <-  gamm(Exponent ~ s(age, k = 3)  + Condition + Region, data = fooofLong, random=list(luna=~1))
summary(gam.model$gam)

lm.model <-  lmerTest::lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooofLong)
summary(lm.model)

AIC(lmer(Exponent ~ age + Condition + Region + (1|luna), data = fooofLong), 
    lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooofLong))


### condition interaction, controlling for region ----
fooofLong$oCondition <- ordered(fooofLong$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)


lm.model <-  lmerTest::lmer(Exponent ~ inverseAge*Condition + Region + (1|luna), data = fooofLong)
summary(lm.model)

# checking to see which condition has a stonger effect
gam.model <-  gamm(Exponent ~ s(age, k = 3) + Region, data = fooofLong %>% filter(Condition == 'eyesClosed'), random=list(luna=~1))
summary(gam.model$gam)


### region interaction, controlling for condition ----

fooofLong$oRegion <- ordered(fooofLong$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
model_formula <- as.formula("Exponent ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)



## Offset vs age ----

lunaize(ggplot(data = fooofLong%>% filter(Condition == 'eyesClosed'), 
               aes(x = age, y = Offset)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Offset") + theme(legend.position = "none")

gam.model <-  gamm(Offset ~ s(age, k = 3)  + Condition + Region, data = fooofLong, random=list(luna=~1))
summary(gam.model$gam)

AIC(lmer(Offset ~ age + Condition + Region + (1|luna), data = fooofLong), 
    lmer(Offset ~ inverseAge + Condition + Region + (1|luna), data = fooofLong))


### region interaction, controlling for condition ----

fooofLong$oRegion <- ordered(fooofLong$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
model_formula <- as.formula("Offset ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)


### condition interaction, controlling for region ----

fooofLong$oCondition <- ordered(fooofLong$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
model_formula <- as.formula("Offset ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)


## Fooof loop ----

### Fooof vs Age ----

fooofVars <- c('Offset','Exponent')

output <- c()  
for (fooofVar in fooofVars) {
  
  model <- paste0(fooofVar, ' ~ ', 's(age, k = 3) + Region + Condition')
  model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofLong))
  
  F <- summary(model.out$gam)$s.table[1,3]
  p <- summary(model.out$gam)$s.table[1,4]
  RegionF <- summary(model.out$gam)$pTerms.table[1,2]
  Regionp <- summary(model.out$gam)$pTerms.table[1,3]
  ConditionF <- summary(model.out$gam)$pTerms.table[2,2]
  Conditionp <- summary(model.out$gam)$pTerms.table[2,3]
  
  output <- rbind(output, data.frame(fooofVar, F, p, RegionF, Regionp, ConditionF, Conditionp))
}
output  


## looking at interactions 
# Interaction of condition 
fooofLong$oCondition <- ordered(fooofLong$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group

output <- c()  
for (fooofVar in fooofVars) {
  
  model_formula <- as.formula(paste0(fooofVar, ' ~ ', ' oCondition + s(age, k = 3, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region')) 
  
  model.out <- (mgcv::gamm(model_formula,
                           random = list(luna=~1),
                           data = fooofLong))
  
  ConditionF <- summary(model.out$gam)$s.table[2,3]
  Conditionp <- summary(model.out$gam)$s.table[2,4]
  
  output <- rbind(output, data.frame(fooofVar, ConditionF, Conditionp))
}
output  


# region interaction
fooofLong$oRegion <- ordered(fooofLong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC reference group

output <- c()  
for (fooofVar in fooofVars) {
  
  model_formula <- as.formula(paste0(fooofVar, ' ~ ', ' oRegion + s(age, k = 3, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition')) 
  
  model.out <- (mgcv::gamm(model_formula,
                           random = list(luna=~1),
                           data = fooofLong))
  
  RegionF <- summary(model.out$gam)$s.table[2,3]
  Regionp <- summary(model.out$gam)$s.table[2,4]
  
  output <- rbind(output, data.frame(fooofVar, RegionF, Regionp))
}
output  

