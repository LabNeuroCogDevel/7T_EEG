
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


# Load Dataframe

entropy <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_allChans_MultiScaleEntropy.csv') %>% select(c(-type))

# Entropy Outlier Detection ----
entropy_outlier <- entropy %>% group_by(Subject) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) %>% ungroup()

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")
ageValues$Subject <- paste(ageValues$lunaID, ageValues$visitDate, sep = "_")

entropyAge <- merge(entropy_outlier, ageValues, by = "Subject")

# Entropy averaged across all electrodes ----
entropyAgeAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, Var1, age) ~ Subject + visitno, data = entropyAge, FUN = mean)
write.csv(entropyAgeAvg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_multiscaleEntropy_allChansAvg.csv')



entropyAgeAvg_outlier <- entropyAgeAvg %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = MSx6)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
  geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


entropyAgeAvg_outlier_long <- entropyAgeAvg_outlier %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))


lunaize(ggplot(data = entropyAgeAvg_outlier_long, aes(x = age, y = MSx, color = timeScale)) + 
          geom_smooth(aes(group = timeScale, color = timeScale), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = entropyAgeAvg_outlier_long, aes(x = age, y = Var1)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("AOC")  


gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3), data = entropyAgeAvg_outlier_long, random=list(lunaID=~1))
summary(gam.model$gam)


# Find each subjects max entropy value on the MSE curve 
subs = unique(entropyAgeAvg_outlier_long$lunaID)
visits = unique(entropyAgeAvg_outlier_long$visitno)


maxValues <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset
  for (visitNum in visits) {
    if (any(entropyAgeAvg_outlier_long$lunaID == subject & entropyAgeAvg_outlier_long$visitno == visitNum)) {
      subData <- entropyAgeAvg_outlier_long %>% filter(lunaID == subject) %>% filter(visitno == visitNum)
      maxSubEntropy <- max(subData$MSx)
      
      subInfo <- data.frame(lunaID = subject, visitno = visitNum, maxEntropy = maxSubEntropy)
      maxValues <- rbind(maxValues, subInfo)
    }
  }
}

maxValuesAge <- merge(maxValues, ageValues, by = c('lunaID', "visitno"))
write.csv(maxValues, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_maxEntropy_allChansAvg.csv')

lunaize(ggplot(data = maxValuesAge, aes(x = age, y = maxEntropy)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge, random=list(lunaID=~1))
summary(gam.model$gam)




# Entropy averaged across only frontal electrodes ----
entropyFrontal <- entropyAge %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2','Fz', 'AF5', 'AF6'))

entropyFrontalAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, Var1, age) ~ Subject + visitno, data = entropyFrontal, FUN = mean)

entropyAgeAvg_outlier_frontal <- entropyFrontalAvg %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


lunaize(ggplot(data = entropyAgeAvg_outlier_frontal, aes(x = age, y = MSx6)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


entropyAgeAvg_outlier_frontal_long <- entropyAgeAvg_outlier_frontal %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))


lunaize(ggplot(data = entropyAgeAvg_outlier_frontal_long, aes(x = age, y = MSx, color = timeScale)) + 
          geom_smooth(aes(group = timeScale, color = timeScale), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = entropyAgeAvg_outlier_frontal_long, aes(x = age, y = Var1)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("AOC")  


gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3), data = entropyAgeAvg_outlier_frontal_long, random=list(lunaID=~1))
summary(gam.model$gam)


# Find each subjects max entropy value on the MSE curve 
subs = unique(entropyAgeAvg_outlier_frontal_long$lunaID)
visits = unique(entropyAgeAvg_outlier_frontal_long$visitno)


maxValuesFrontal <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset
  for (visitNum in visits) {
    if (any(entropyAgeAvg_outlier_frontal_long$lunaID == subject & entropyAgeAvg_outlier_frontal_long$visitno == visitNum)) {
      subData <- entropyAgeAvg_outlier_frontal_long %>% filter(lunaID == subject) %>% filter(visitno == visitNum)
      maxSubEntropy <- max(subData$MSx)
      
      subInfo <- data.frame(lunaID = subject, visitno = visitNum, maxEntropy = maxSubEntropy)
      maxValuesFrontal <- rbind(maxValuesFrontal, subInfo)
    }
  }
}

maxValuesAgeFrontal <- merge(maxValuesFrontal, ageValues, by = c('lunaID', "visitno"))


lunaize(ggplot(data = maxValuesAge, aes(x = age, y = maxEntropy)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge, random=list(lunaID=~1))
summary(gam.model$gam)


# Combine Entropy and Spontaneous by PCA ----
allSNRmeasures <- read.csv('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Lab/Projects/SNR/rMarkdown/SNRmeasures_PCAvalues.csv')

maxEntropySelect <- merge(entropyAgeAvg_outlier %>% select(lunaID, visitno, age, Var1), maxValuesAge, by = c("lunaID", "visitno", "age"))

spontaneousEntropy <- merge(maxEntropy_frontal_AOC, allSNRmeasures, by = c("lunaID", "visitno", "age")) %>% mutate(ageGroup = as.factor(cut(age, c(0,17,Inf), labels = c('Adol','Adults')))) 


lunaize(ggplot(data = spontaneousEntropy, aes(x = pc1, y = Var1)) + geom_point(aes(color = ageGroup)) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group =1, method="lm", alpha=0.4, linewidth=2)) + ylab("Area under curve") + facet_wrap(~measure)

lm.model <-  lmer(Var1 ~ age + pc1 + (1|lunaID), data = spontaneousEntropy)
summary(lm.model)

gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3), data = spontaneousEntropy, random=list(lunaID=~1))
summary(gam.model$gam)




lunaize(ggplot(data = spontaneousEntropy, aes(x = pc1, y = maxEntropy)) + geom_point(aes(color = ageGroup)) + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method = "lm", alpha=0.4, linewidth=2)) + facet_wrap(~measure) 

lm.model <-  lmer(maxEntropy ~ age + pc1 + (1|lunaID), data = spontaneousEntropy)
summary(lm.model)          

gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3) + pc1, data = spontaneousEntropy, random=list(lunaID=~1))
summary(gam.model$gam)       



spontaneousEntropy_long <- spontaneousEntropy %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>% mutate(ageGroup = as.factor(cut(age, c(0,17,Inf), labels = c('Adol','Adults')))) 
                            

lunaize(ggplot(data = spontaneousEntropy_long, aes(x = pc1, y = MSx, color = timeScale)) + 
          geom_smooth(aes(group = timeScale, color = timeScale), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


# Combine Entropy and Spontaneous by All Channel ----
## Find max entropy for every electrode ----
# Find each subjects max entropy value on the MSE curve 

entropyAge_long <- entropyAge %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))

subs = unique(entropyAge_long$Subject)
visits <- unique(na.omit(entropyAge_long$visitno))

maxValuesAllChans <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset
  
      subData <- entropyAge_long %>% filter(Subject == subject)
      chans = unique(subData$labels)
      
      for (chan in chans) {
        # Check if the subject has the visit number in the dataset
        
        subChanData <- subData %>% filter(labels == chan)
        maxSubchanEntropy <- max(subChanData$MSx)
        
        subChanInfo <- data.frame(Subject = subject, labels = chan, maxEntropy = maxSubchanEntropy)
        maxValuesAllChans <- rbind(maxValuesAllChans, subChanInfo)
      }
}



maxValuesAllChansAge <- merge(maxValuesAllChans, ageValues, by = c('Subject')) %>% merge(., entropyAge %>% select("Subject", "labels", "Var1"), by = c("Subject", "labels"))

SNRallChans_outlier <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_allChans_allfreqs.csv')

SNRallChans_outlier4040 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 40)
SNRallChans_outlier4040$SNR <- log(SNRallChans_outlier4040$Evoked/SNRallChans_outlier4040$Induced)

spontaneousEntropy <- merge(maxValuesAllChansAge, SNRallChans_outlier4040 %>% 
                              select("lunaID", "visitno", "age", "labels", "Induced", "Evoked", "SNR", "Total"), by = c("lunaID", "visitno", "age", "labels")) %>% 
  mutate(ageGroup = as.factor(cut(age, c(0,17,Inf), labels = c('Adol','Adults')))) 

hist(log(spontaneousEntropy$Induced))
hist((spontaneousEntropy$maxEntropy))

entropyAgeAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, Var1, age) ~ Subject + visitno, data = entropyAge, FUN = mean)

spontaneousEntropyAvg <- aggregate(cbind(Induced, Evoked, Total, SNR, Var1, maxEntropy, age) ~ lunaID + visitno, data = spontaneousEntropy, FUN = mean)


lunaize(ggplot(data = spontaneousEntropyAvg, aes(x = log(Induced), y = Var1)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group =1, method="lm", alpha=0.4, linewidth=2)) + ylab("Area under curve") 

lm.model <-  lmer(Var1 ~ age + log(Induced) + (1|lunaID), data = spontaneousEntropyAvg)
summary(lm.model)

lunaize(ggplot(data = spontaneousEntropyAvg, aes(x = log(Induced), y = maxEntropy)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group =1, method="lm", alpha=0.4, linewidth=2)) + ylab("Max Entropy") 

lm.model <-  lmer(maxEntropy ~ age + log(Induced) + (1|lunaID), data = spontaneousEntropyAvg)
summary(lm.model)

## Correlation between entropy and spontaneous ----
subs = unique(spontaneousEntropy$Subject)

subCors <- data.frame()

for (subject in subs) {
  
  subDataFrame <- spontaneousEntropy %>% filter(Subject == subject)
  subjectCor <- cor(x = scale(subDataFrame$maxEntropy), y = log(subDataFrame$Induced), method = "pearson", use = "complete.obs")
  
  subCorInfo <- data.frame(Subject = subject, Corr = subjectCor)
  subCors <- rbind(subCors, subCorInfo)
  
}

subCorAge <- merge(subCors, ageValues, by = c("Subject"))

lunaize(ggplot(data = subCorAge, aes(x = age, y = Corr)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method="lm", alpha=0.4, linewidth=2))


lm.model <-  lmer(Corr ~ age + (1|lunaID), data = subCorAge)
summary(lm.model)

# MSE Entropy vs FOOOF ----

entropy <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Entropy/allSubjects_allChans_MultiScaleEntropy.csv') %>% select(c(-type))

entropy_outlier <- entropy %>% group_by(Subject) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) %>% ungroup()

entropyAge <- merge(entropy_outlier, ageValues, by = "Subject")

## Entropy averaged across all electrodes ----
entropyAgeAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, Var1, age) ~ Subject + visitno, data = entropyAge, FUN = mean)

entropyAgeAvg_outlier <- entropyAgeAvg %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) 

entropyAgeAvg_outlier_long <- entropyAgeAvg_outlier %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))

## Load and avg fooof across all electrodes ----
fooofAllchans <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectsAllChannelsFooofMeasures_20230911.csv')

fooofAllchans <- fooofAllchans %>% 
  rename(labels = Channel)

fooofAvg <- aggregate(cbind(Offset, Exponent) ~ Subject + Condition, data = fooofAllchans, FUN = mean)

fooofAvg_outlier <- fooofAvg %>% group_by(Condition) %>%
  mutate(across(c("Offset", "Exponent"), naoutlier)) %>% ungroup()

## Merge MSE and fooof ----

MSEfooof <- merge(fooofAvg_outlier, entropyAgeAvg_outlier_long, by = "Subject") %>% separate(Subject, c('lunaID','vdate'))


lunaize(ggplot(data = MSEfooof, aes(x = Exponent, y = MSx))+ 
          geom_smooth(aes(group = timeScale, color = timeScale), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = MSEfooof, aes(y = Exponent, x = Var1))+ geom_point() +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = MSEfooof, aes(y = Offset, x = Var1))+ geom_point() +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


# Spectral Entropy ----

entropy <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Entropy/allSubjects_allChans_SpectralEntropy_broadband.csv') %>% select(c(-type))

entropy_outlier <- entropy %>% group_by(Subject) %>%
  mutate(across(c("gammaBandEn", "betaBandEn", "alphaBandEn", "thetaBandEn", "Spec"), naoutlier)) %>% ungroup()

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")
ageValues$Subject <- paste(ageValues$lunaID, ageValues$visitDate, sep = "_")

entropyAge <- merge(entropy_outlier, ageValues, by = "Subject")

##Averaged across all electrodes ----
entropyAgeAvg <- aggregate(cbind(gammaBandEn, betaBandEn, alphaBandEn, thetaBandEn, Spec, age) ~ Subject + visitno, data = entropyAge, FUN = mean)

entropyAgeAvg_outlier <- entropyAgeAvg %>%
  mutate(across(c("gammaBandEn", "betaBandEn", "alphaBandEn", "thetaBandEn"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = gammaBandEn)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = betaBandEn)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = alphaBandEn)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = thetaBandEn)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = Spec)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

## All channels vs FOOOF ----
entropy <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Entropy/allSubjects_allChans_SpectralEntropy_broadband.csv') %>% select(c(-type))

### Avg channels ----
entropyAvg <- aggregate(cbind(gammaBandEn, betaBandEn, alphaBandEn, thetaBandEn, Spec) ~ Subject, data = entropy, FUN = mean)

entropyAvg_outlier <- entropyAvg %>%
  mutate(across(c("gammaBandEn", "betaBandEn", "alphaBandEn", "thetaBandEn", "Spec"), naoutlier)) %>% ungroup()


fooofAllchans <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectsAllChannelsFooofMeasures_20230911.csv')

fooofAllchans <- fooofAllchans %>% 
  rename(labels = Channel)

fooofAvg <- aggregate(cbind(Offset, Exponent) ~ Subject + Condition, data = fooofAllchans, FUN = mean)

fooofAvg_outlier <- fooofAvg %>% group_by(Condition) %>%
  mutate(across(c("Offset", "Exponent"), naoutlier)) %>% ungroup()

### Merge Entropy and FOOOF ----
entropyFooof <- merge(entropyAvg_outlier, fooofAvg_outlier, by = c("Subject")) %>% separate(Subject, c('lunaID','vdate'))


lunaize(ggplot(data = entropyFooof, aes(x = gammaBandEn, y = Exponent)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = Condition), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = entropyFooof, aes(x = gammaBandEn, y = Offset)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = entropyFooof, aes(x = Spec, y = Exponent)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + facet_wrap(~Condition)


