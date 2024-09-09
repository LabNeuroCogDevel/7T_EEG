
# Load Libraries 

library(LNCDR)
library(dplyr)
library(ggplot2)
library(mgcv)
library(lme4)
library(lmerTest)
library(tidyr)


outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)
# Prep Dataframe ----
## Load datasets ----

MSE <- read.csv( '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_allChans_MultiScaleEntropy.csv')

## Entropy Outlier Detection ----
entropy_outlier <- MSE %>% group_by(Subject) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) %>% ungroup()

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")
ageValues$Subject <- paste(ageValues$lunaID, ageValues$visitDate, sep = "_")

entropyAge <- merge(entropy_outlier, ageValues, by = "Subject")

## Entropy averaged across only DLPFC electrodes ----
entropyLDLPFC <- entropyAge %>% filter(labels %in% c('F3', 'F5', 'F7'))
entropyRDLPFC <- entropyAge %>% filter(labels %in% c('F4', 'F6', 'F8'))

entropyLDLPFCAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, Var1, age) ~ Subject + visitno, data = entropyLDLPFC, FUN = mean)
entropyRDLPFCAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, Var1, age) ~ Subject + visitno, data = entropyRDLPFC, FUN = mean)

entropyDLPFC <- merge(entropyLDLPFCAvg, entropyRDLPFCAvg, by = c("Subject", "visitno","age"), suffixes = c("_LDLPFC", "_RDLPFC"))
  
  
entropyAgeAvg_outlier_DLPFC <- entropyDLPFC %>%
  mutate(across(c("MSx1_LDLPFC", "MSx2_LDLPFC", 
"MSx3_LDLPFC", "MSx4_LDLPFC", "MSx5_LDLPFC", "MSx6_LDLPFC", "MSx7_LDLPFC", 
"MSx8_LDLPFC", "MSx9_LDLPFC", "MSx10_LDLPFC", "Var1_LDLPFC", 
"MSx1_RDLPFC", "MSx2_RDLPFC", "MSx3_RDLPFC", "MSx4_RDLPFC", "MSx5_RDLPFC", 
"MSx6_RDLPFC", "MSx7_RDLPFC", "MSx8_RDLPFC", "MSx9_RDLPFC", "MSx10_RDLPFC", 
"Var1_RDLPFC"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))

entropyAgeAvg_outlier_long <- entropyAgeAvg_outlier_DLPFC %>%
  pivot_longer(
    cols = starts_with("MSx"),
    names_to = c("timescale", "region"),
    names_pattern = "MSx(\\d+)_(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = region,
    values_from = value,
    names_prefix = "MSx_"
  )


# Find each subjects max entropy value on the MSE curve ----
subs = unique(entropyAgeAvg_outlier_long$lunaID)
visits = unique(entropyAgeAvg_outlier_long$visitno)


maxValues <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset
  for (visitNum in visits) {
    if (any(entropyAgeAvg_outlier_long$lunaID == subject & entropyAgeAvg_outlier_long$visitno == visitNum)) {
      subData <- entropyAgeAvg_outlier_long %>% filter(lunaID == subject) %>% filter(visitno == visitNum)
      maxSubEntropy_LDLPFC <- max(subData$MSx_LDLPFC)
      maxSubEntropy_RDLPFC <- max(subData$MSx_RDLPFC)
      
      subInfo <- data.frame(lunaID = subject, visitno = visitNum, maxEntropy_LDLPFC = maxSubEntropy_LDLPFC, maxEntropy_RDLPFC = maxSubEntropy_RDLPFC)
      maxValues <- rbind(maxValues, subInfo)
    }
  }
}

entropyValuesMax <- merge(entropyAgeAvg_outlier_DLPFC, maxValues, by = c("lunaID", "visitno"))

## Extract DLPFC Hurst from merge 7t ----

hurst <- merge7t[c("lunaid","visitno","rest.hurst.LDLPFC", "rest.hurst.RDLPFC")]
names(hurst)[names(hurst) == 'lunaid'] <- 'lunaID'
hurst$avgDLPFChurst <- rowMeans(hurst[, c("rest.hurst.LDLPFC", "rest.hurst.RDLPFC")], na.rm = TRUE)

## Merge Entropy and Hurst ----

entropyHurst <- merge(entropyValuesMax, hurst, by = c("lunaID", "visitno"))
names(entropyHurst)[names(entropyHurst) == 'Var1_LDLPFC'] <- 'AUC_LDLPFC'
names(entropyHurst)[names(entropyHurst) == 'Var1_RDLPFC'] <- 'AUC_RDLPFC'


# Max Entropy vs Hurst ----
## LDLPFC ----
### Gam ----
lunaize(ggplot(data = entropyHurst, aes(x = maxEntropy_LDLPFC, y = rest.hurst.LDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) 


gam.model <-  mgcv::gamm(rest.hurst.LDLPFC ~ maxEntropy_LDLPFC + s(age, k = 3), data = entropyHurst, random=list(lunaID=~1))
summary(gam.model$gam)

### Linear ----
lunaize(ggplot(data = entropyHurst, aes(x = maxEntropy_LDLPFC, y = rest.hurst.LDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method="lm", alpha=0.4, linewidth=2)) 


lm.model <-  lmer(rest.hurst.LDLPFC ~ maxEntropy_LDLPFC + age + (1|lunaID), data = entropyHurst)
summary(lm.model)

## RDLPFC ----
### Gam ----
lunaize(ggplot(data = entropyHurst, aes(x = maxEntropy_RDLPFC, y = rest.hurst.RDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) 


gam.model <-  mgcv::gamm(rest.hurst.RDLPFC ~ maxEntropy_RDLPFC + s(age, k = 3), data = entropyHurst, random=list(lunaID=~1))
summary(gam.model$gam)

### Linear ----
lunaize(ggplot(data = entropyHurst, aes(x = maxEntropy_RDLPFC, y = rest.hurst.RDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method="lm", alpha=0.4, linewidth=2)) 


lm.model <-  lmer(rest.hurst.RDLPFC ~ maxEntropy_RDLPFC + age + (1|lunaID), data = entropyHurst)
summary(lm.model)



# AUC vs Hurst ----
## LDLPFC ----
### Gam ----
lunaize(ggplot(data = entropyHurst, aes(x = AUC_LDLPFC, y = rest.hurst.LDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) 


gam.model <-  mgcv::gamm(rest.hurst.LDLPFC ~ AUC_LDLPFC + s(age, k = 3), data = entropyHurst, random=list(lunaID=~1))
summary(gam.model$gam)

### Linear ----
lunaize(ggplot(data = entropyHurst, aes(x = AUC_LDLPFC, y = rest.hurst.LDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method="lm", alpha=0.4, linewidth=2)) 


lm.model <-  lmer(rest.hurst.LDLPFC ~ AUC_LDLPFC + age + (1|lunaID), data = entropyHurst)
summary(lm.model)

## RDLPFC ----
### Gam ----
lunaize(ggplot(data = entropyHurst, aes(x = AUC_RDLPFC, y = rest.hurst.RDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) 


gam.model <-  mgcv::gamm(rest.hurst.RDLPFC ~ AUC_RDLPFC + s(age, k = 3), data = entropyHurst, random=list(lunaID=~1))
summary(gam.model$gam)

### Linear ----
lunaize(ggplot(data = entropyHurst, aes(x = AUC_RDLPFC, y = rest.hurst.RDLPFC)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method="lm", alpha=0.4, linewidth=2)) 


lm.model <-  lmer(rest.hurst.RDLPFC ~ AUC_RDLPFC + age + (1|lunaID), data = entropyHurst)
summary(lm.model)


