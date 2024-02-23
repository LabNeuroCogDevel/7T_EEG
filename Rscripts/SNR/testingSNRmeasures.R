

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
library(corrplot)

# Define outlier functions ----

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

# Load in dataframes ----
merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]

colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")

snr40 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/attic/SNRdata_40Hz.csv', head = T)
snr40 <- snr40 %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))

snr40_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/attic/SNRdata_40Hz_method2.csv', head = T)
snr40_2 <- snr40_2 %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))

combineSNRmethods <- merge(snr40, snr40_2, by = c("lunaID", "visitDate","visitno", "age"), suffixes = c("","_2"))

# Find rows where itc is 0
rows_to_delete <- which(combineSNRmethods$SNR == 0 | combineSNRmethods$SNR_2 == 0)

# Delete rows with itc equal to 0
combineSNRmethods <- combineSNRmethods[-rows_to_delete, ]


### Outlier Detection Subject Level 
combineSNRmethods_naout <- combineSNRmethods %>% 
  mutate(across(c("SNR", "SNR_2"), naoutlier))


# SNR vs age ----

cowplot::plot_grid(
  ggplot(data = combineSNRmethods_naout, aes(x = age, y = SNR)) + 
    geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)),
  ggplot(data = combineSNRmethods_naout, aes(x = age, y = SNR_2)) + 
    geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
    geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T))) 

cor(combineSNRmethods_naout$SNR_2, combineSNRmethods_naout$SNR, method = "pearson", use = "complete.obs")


# Load in ITC ----
itc40 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/attic/allSubsITC_40Hz.csv')

# Find rows where itc is 0
rows_to_delete <- which(itc40$ITC == 0)

# Delete rows with itc equal to 0
itc40 <- itc40[-rows_to_delete, ]


### Outlier Detection Subject Level 
itc40_naout <- itc40 %>% group_by(lunaID) %>%
  mutate(across(c("ITC"), naoutlier)) %>% ungroup()


aggregated_itc40 <- aggregate(ITC ~ lunaID+visitDate, itc40_naout, mean)

snr_itc40 <- merge(combineSNRmethods_naout,aggregated_itc40, by = c("lunaID", "visitDate") )


# SNR vs ITC ----

ggplot(data = snr_itc40, aes(x = ITC, y = SNR_2)) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method='lm') 

summary(lmer(data = snr_itc40, SNR ~ ITC +(1|lunaID)))

cor(snr_itc40$ITC, snr_itc40$SNR, method = "pearson")

# Total Evoked Induced ----

totalEvokedInduced <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/SNRdata_40Hz_totalEvokedInduced.csv')
totalEvokedInduced <- totalEvokedInduced %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))

# Outlier Detection Subject Level 

totalEvokedInduced <- totalEvokedInduced %>% mutate(ageGroup = cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))

# Find rows where the values are 0
rows_to_delete <- which(totalEvokedInduced$TotalPower == 0 | totalEvokedInduced$EvokedPower == 0 | totalEvokedInduced$InducedPower == 0)

# Delete rows with the values equal to 0
totalEvokedInduced <- totalEvokedInduced[-rows_to_delete, ]

totalEvokedInduced <- totalEvokedInduced %>%
  mutate(across(c("TotalPower", "EvokedPower","InducedPower"), naoutlier))

# Convert back to amp
totalEvokedInduced$EvokedPower <- 10^((totalEvokedInduced$EvokedPower)/10)
totalEvokedInduced$TotalPower <- 10^((totalEvokedInduced$TotalPower)/10)

totalEvokedInduced$InducedPower <- totalEvokedInduced$TotalPower - totalEvokedInduced$EvokedPower


totalEvokedInduced$SNR <- totalEvokedInduced$EvokedPower/totalEvokedInduced$InducedPower
totalEvokedInduced$SNR2 <- (totalEvokedInduced$EvokedPower-totalEvokedInduced$InducedPower)/(totalEvokedInduced$EvokedPower+totalEvokedInduced$InducedPower)


totalEvokedInduced <- totalEvokedInduced %>%
  mutate(across(c("SNR"), naoutlier))

## Total vs age ----
total <- ggplot(data = totalEvokedInduced, aes(x = age, y = 10*log10(TotalPower))) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Total Power (dB)")

gam.model <-  gamm(TotalPower ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

## Evoked vs age ----
evoked <- ggplot(data = totalEvokedInduced, aes(x = age, y = 10*log10(EvokedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Evoked Power (dB)")

gam.model <-  gamm(EvokedPower ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

## Induced vs age ----
induced <- ggplot(data = totalEvokedInduced, aes(x = age, y = 10*log10(InducedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Induced Power (dB)")

gam.model <-  gamm(InducedPower ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

## SNR vs age ----
snr <- ggplot(data = totalEvokedInduced, aes(x = age, y = SNR)) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("SNR")

gam.model <-  gamm(SNR ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

cowplot::plot_grid(total, evoked, induced, snr)


snr_itc40 <- merge(totalEvokedInduced,aggregated_itc40, by = c("lunaID", "visitDate") )

## SNR vs ITC ----

ggplot(data = snr_itc40, aes(x = ITC, y = SNR)) + 
  geom_point() + stat_smooth(method='lm') 



# Total Evoked Induced No Hanning ----

totalEvokedInduced <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/SNRdata_40Hz_totalEvokedInduced_noHanning.csv')
totalEvokedInduced <- totalEvokedInduced %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))

# Outlier Detection Subject Level 

totalEvokedInduced <- totalEvokedInduced %>% mutate(ageGroup = cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))

# Find rows where the values are 0
rows_to_delete <- which(totalEvokedInduced$TotalPower == 0 | totalEvokedInduced$EvokedPower == 0 | totalEvokedInduced$InducedPower == 0)

# Delete rows with the values equal to 0
totalEvokedInduced <- totalEvokedInduced[-rows_to_delete, ]

totalEvokedInduced <- totalEvokedInduced %>%
  mutate(across(c("TotalPower", "EvokedPower","InducedPower"), naoutlier))

totalEvokedInduced$SNR <- totalEvokedInduced$EvokedPower/totalEvokedInduced$InducedPower
totalEvokedInduced$SNR2 <- (totalEvokedInduced$EvokedPower-totalEvokedInduced$InducedPower)/(totalEvokedInduced$EvokedPower+totalEvokedInduced$InducedPower)


totalEvokedInduced <- totalEvokedInduced %>%
  mutate(across(c("SNR", "SNR2"), naoutlier))

## Total vs age ----
total <- ggplot(data = totalEvokedInduced, aes(x = age, y = 10*log10(TotalPower))) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Total Power (dB)")

gam.model <-  gamm(TotalPower ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

## Evoked vs age ----
evoked <- ggplot(data = totalEvokedInduced, aes(x = age, y = 10*log10(EvokedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Evoked Power (dB)")

gam.model <-  gamm(EvokedPower ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

## Induced vs age ----
induced <- ggplot(data = totalEvokedInduced, aes(x = age, y = 10*log10(InducedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Induced Power (dB)")

gam.model <-  gamm(InducedPower ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

## SNR vs age ----
snr <- ggplot(data = totalEvokedInduced, aes(x = age, y = SNR)) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("SNR")

gam.model <-  gamm(SNR ~ s(age, k = 3), data = totalEvokedInduced, random=list(lunaID=~1))
summary(gam.model$gam)

cowplot::plot_grid(total, evoked, induced, snr)


snr_itc40 <- merge(totalEvokedInduced,aggregated_itc40, by = c("lunaID", "visitDate") )

## SNR vs ITC ----

ggplot(data = snr_itc40, aes(x = ITC, y = SNR)) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method='lm') 

gam.model <-  gamm(SNR ~ s(ITC, k = 3), data = snr_itc40, random=list(lunaID=~1))
summary(gam.model$gam)

# Total Evoked Induced Baseline ----

totalEvokedInduced_baseline <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/SNRdata_40Hz_totalEvokedInduced_baseline.csv')
totalEvokedInduced_baseline <- totalEvokedInduced_baseline %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))

# Outlier Detection Subject Level 

totalEvokedInduced_baseline <- totalEvokedInduced_baseline %>% mutate(ageGroup = cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))

# Find rows where the values are 0
rows_to_delete <- which(totalEvokedInduced_baseline$TotalPower == 0 | totalEvokedInduced_baseline$EvokedPower == 0 | totalEvokedInduced_baseline$InducedPower == 0)

# Delete rows with the values equal to 0
totalEvokedInduced_baseline <- totalEvokedInduced_baseline[-rows_to_delete, ]

totalEvokedInduced_baseline <- totalEvokedInduced_baseline %>%
  mutate(across(c("TotalPower", "EvokedPower","InducedPower"), naoutlier))

totalEvokedInduced_baseline$SNR <- totalEvokedInduced_baseline$EvokedPower/totalEvokedInduced_baseline$InducedPower
totalEvokedInduced_baseline$SNR2 <- (totalEvokedInduced_baseline$EvokedPower-totalEvokedInduced_baseline$InducedPower)/(totalEvokedInduced_baseline$EvokedPower+totalEvokedInduced_baseline$InducedPower)

totalEvokedInduced_baseline <- totalEvokedInduced_baseline %>%
  mutate(across(c("SNR", "SNR2"), naoutlier))

## Total vs age ----
total <- ggplot(data = totalEvokedInduced_baseline, aes(x = age, y = 10*log10(TotalPower))) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Total Power (dB)")

gam.model <-  gamm(TotalPower ~ s(age, k = 3), data = totalEvokedInduced_baseline, random=list(lunaID=~1))
summary(gam.model$gam)

## Evoked vs age ----
evoked <- ggplot(data = totalEvokedInduced_baseline, aes(x = age, y = 10*log10(EvokedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Evoked Power (dB)")

gam.model <-  gamm(EvokedPower ~ s(age, k = 3), data = totalEvokedInduced_baseline, random=list(lunaID=~1))
summary(gam.model$gam)

## Induced vs age ----
induced <- ggplot(data = totalEvokedInduced_baseline, aes(x = age, y = 10*log10(InducedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Induced Power (dB)")

gam.model <-  gamm(InducedPower ~ s(age, k = 3), data = totalEvokedInduced_baseline, random=list(lunaID=~1))
summary(gam.model$gam)

## SNR vs age ----
snr <- ggplot(data = totalEvokedInduced_baseline, aes(x = age, y = SNR)) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("SNR")

gam.model <-  gamm(SNR ~ s(age, k = 3), data = totalEvokedInduced_baseline, random=list(lunaID=~1))
summary(gam.model$gam)

cowplot::plot_grid(total, evoked, induced, snr)

# Total Evoked Induced Baseline Removed ----

totalEvokedInduced_baselineRemoved <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/SNRdata_40Hz_totalEvokedInduced_baselineRemoved.csv')
totalEvokedInduced_baselineRemoved <- totalEvokedInduced_baselineRemoved %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))

# Outlier Detection Subject Level 

totalEvokedInduced_baselineRemoved <- totalEvokedInduced_baselineRemoved %>% mutate(ageGroup = cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))

# Find rows where the values are 0
rows_to_delete <- which(totalEvokedInduced_baselineRemoved$TotalPower == 0 | totalEvokedInduced_baselineRemoved$EvokedPower == 0 | totalEvokedInduced_baselineRemoved$InducedPower == 0)

# Delete rows with the values equal to 0
totalEvokedInduced_baselineRemoved <- totalEvokedInduced_baselineRemoved[-rows_to_delete, ]

totalEvokedInduced_baselineRemoved <- totalEvokedInduced_baselineRemoved %>%
  mutate(across(c("TotalPower", "EvokedPower","InducedPower"), naoutlier))

totalEvokedInduced_baselineRemoved$SNR <- totalEvokedInduced_baselineRemoved$EvokedPower/totalEvokedInduced_baselineRemoved$InducedPower
totalEvokedInduced_baselineRemoved$SNR2 <- (totalEvokedInduced_baselineRemoved$EvokedPower-totalEvokedInduced_baselineRemoved$InducedPower)/(totalEvokedInduced_baselineRemoved$EvokedPower+totalEvokedInduced_baselineRemoved$InducedPower)

totalEvokedInduced_baselineRemoved <- totalEvokedInduced_baselineRemoved %>%
  mutate(across(c("SNR", "SNR2"), naoutlier))

## Total vs age ----
total <- ggplot(data = totalEvokedInduced_baselineRemoved, aes(x = age, y = 10*log10(TotalPower))) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Total Power (dB)")

gam.model <-  gamm(TotalPower ~ s(age, k = 3), data = totalEvokedInduced_baselineRemoved, random=list(lunaID=~1))
summary(gam.model$gam)

## Evoked vs age ----
evoked <- ggplot(data = totalEvokedInduced_baselineRemoved, aes(x = age, y = 10*log10(EvokedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Evoked Power (dB)")

gam.model <-  gamm(EvokedPower ~ s(age, k = 3), data = totalEvokedInduced_baselineRemoved, random=list(lunaID=~1))
summary(gam.model$gam)

## Induced vs age ----
induced <- ggplot(data = totalEvokedInduced_baselineRemoved, aes(x = age, y = 10*log10(InducedPower))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Induced Power (dB)")

gam.model <-  gamm(InducedPower ~ s(age, k = 3), data = totalEvokedInduced_baselineRemoved, random=list(lunaID=~1))
summary(gam.model$gam)

## SNR vs age ----
snr <- ggplot(data = totalEvokedInduced_baselineRemoved, aes(x = age, y = SNR)) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("SNR")

gam.model <-  gamm(SNR ~ s(age, k = 3), data = totalEvokedInduced_baselineRemoved, random=list(lunaID=~1))
summary(gam.model$gam)

cowplot::plot_grid(total, evoked, induced, snr)


# Total Evoked Induced Baseline Removed across Frequencies ----

TFI_40 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsTotalEvokedInduced_40hz.csv')
TFI_40 <- TFI_40 %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))

TFI_40 <- TFI_40 %>% mutate(ageGroup = cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))

TFI_40_groups <- TFI_40 %>% group_by(freqs) %>%
  summarize(meanTotal = mean(Total, na.rm=T),
            meanEvoked = mean(Evoked, na.rm=T), 
            meanInduced = mean(Induced, na.rm=T))

TFI_40_groups$evokedplusInduced <- TFI_40_groups$meanEvoked + TFI_40_groups$meanInduced

ggplot(data = TFI_40_groups%>% filter(freqs >0) , aes(x = freqs, y = evokedplusInduced)) + 
  geom_point(color = "blue") + geom_line(color = "blue") + xlab("Freqs") + ylab("Power")

ggplot(data = TFI_40 %>% filter(freqs > 0 & freqs < 75) %>% filter(lunaID == '10129'), aes(x = freqs, y = Total)) + 
  geom_point(color = "blue") + geom_line(color = "blue") + xlab("Freqs") + ylab("Power")

ggplot() +
  geom_point(data = TFI_40%>% filter(freqs > 0 & freqs < 75), aes(freqs, Total), color = "blue") +
  geom_line(data = TFI_40%>% filter(freqs > 0 & freqs < 75), aes(freqs, Total), color = "blue") +
  geom_point(data = TFI_40%>% filter(freqs > 0 & freqs < 75), aes(freqs, Evoked), color = "green") +
  geom_line(data = TFI_40%>% filter(freqs > 0 & freqs < 75), aes(freqs, Evoked), color = "green") +
  geom_point(data = TFI_40%>% filter(freqs > 0 & freqs < 75), aes(freqs, Induced), color = "red") +
  geom_line(data = TFI_40%>% filter(freqs > 0 & freqs < 75) , aes(freqs, Induced), color = "red") +
  xlab("Freqs") +
  ylab("Power") +
  theme_minimal() + facet_wrap(~ageGroup)


ggplot() +
  geom_point(data = TFI_40_groups%>% filter(freqs > 30 & freqs < 50), aes(freqs, meanTotal), color = "blue") +
  geom_line(data = TFI_40_groups%>% filter(freqs > 30 & freqs < 50), aes(freqs, meanTotal), color = "blue") +
  geom_point(data = TFI_40_groups%>% filter(freqs > 30 & freqs < 50), aes(freqs, meanEvoked), color = "green") +
  geom_line(data = TFI_40_groups%>% filter(freqs > 30 & freqs < 50), aes(freqs, meanEvoked), color = "green") +
  geom_point(data = TFI_40_groups%>% filter(freqs > 30 & freqs < 50), aes(freqs, meanInduced), color = "red") +
  geom_line(data = TFI_40_groups%>%filter(freqs > 30 & freqs < 50) , aes(freqs, meanInduced), color = "red") +
  xlab("Freqs") +
  ylab("Power") +
  theme_minimal() + facet_wrap(~ageGroup)

### Outlier Detection Subject Level 
TFI_40 <- TFI_40 %>% group_by(lunaID) %>%
  mutate(across(c("Total", "Evoked", "Induced"), naoutlier))

TFI_40_aggregate <- TFI_40 %>% filter(freqs >= 38 & freqs <= 42) %>% aggregate(. ,cbind(Total, Evoked, Induced, age) ~ lunaID + visitno + visitDate, mean)

### Outlier Detection Subject Level 
TFI_40_aggregate <- TFI_40_aggregate %>% 
  mutate(across(c("Total", "Evoked", "Induced"), naoutlier))

TFI_40_aggregate$SNR <- TFI_40_aggregate$Evoked/TFI_40_aggregate$Induced

TFI_40_aggregate <- TFI_40_aggregate %>% 
  mutate(across(c("SNR"), naoutlier))


## Total vs age ----
total <- ggplot(data = TFI_40_aggregate, aes(x = age, y = 10*log10(Total))) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Total Power (dB)")

gam.model <-  gamm(Total ~ s(age, k = 3), data = TFI_40_aggregate, random=list(lunaID=~1))
summary(gam.model$gam)

## Evoked vs age ----
evoked <- ggplot(data = TFI_40_aggregate, aes(x = age, y = 10*log10(Evoked))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Evoked Power (dB)")

gam.model <-  gamm(Evoked ~ s(age, k = 3), data = TFI_40_aggregate, random=list(lunaID=~1))
summary(gam.model$gam)

## Induced vs age ----
induced <- ggplot(data = TFI_40_aggregate, aes(x = age, y = 10*log10(Induced))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("Induced Power (dB)")

gam.model <-  gamm(Induced ~ s(age, k = 3), data = TFI_40_aggregate, random=list(lunaID=~1))
summary(gam.model$gam)

## SNR vs age ----
snr <- ggplot(data = TFI_40_aggregate, aes(x = age, y = sqrt(SNR))) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + xlab("Age") + ylab("SNR")

gam.model <-  gamm(sqrt(SNR) ~ s(age, k = 3), data = TFI_40_aggregate, random=list(lunaID=~1))
summary(gam.model$gam)

cowplot::plot_grid(total, evoked, induced, snr)


# Merge all SNR measures ----

TFI_40_aggregate_allFreqs <- TFI_40 %>% filter(freqs > 0) %>% aggregate(. ,cbind(Total, Evoked, Induced, age) ~ lunaID + visitno, mean)
TFI_40_aggregate_allFreqs$SNR <- TFI_40_aggregate_allFreqs$Evoked/TFI_40_aggregate_allFreqs$Induced

allMeasures <- merge(snr40, snr40_2, by = c("lunaID", "visitno", "age", "visitDate"), suffixes = c("_1", "_2")) %>% 
  merge(., TFI_40_aggregate_allFreqs, by = c("lunaID", "visitno", "age"))%>% merge(., aggregated_itc40, by = c("lunaID", "visitDate"))

### Outlier Detection Subject Level 
allMeasures <- allMeasures %>% 
  mutate(across(c("SNR_1", "SNR_2", "SNR"), naoutlier))

corrplot(cor(allMeasures[ ,c("SNR_1", "SNR_2", "SNR", "ITC")], use = "complete.obs"), method = 'number')

ggplot(data = allMeasures, aes(x= SNR, y = SNR_1)) + geom_point()+ stat_smooth(method = "lm")


# ITC and TFI over only 38-42 hz 
aggregated_itc40_gamma <- aggregate(ITC ~ lunaID+visitDate, itc40_naout %>% filter(freqs >= 38 & freqs <= 42), mean)
TFI_ITC <- merge(TFI_40_aggregate, aggregated_itc40_gamma,  by = c("lunaID", "visitDate"))
ggplot(data = TFI_ITC, aes(x= SNR, y = ITC)) + geom_point()+ stat_smooth(method = "lm")
cor(TFI_ITC[ ,c("SNR", "ITC")], use = "complete.obs")

ggplot(data = TFI_ITC, aes(x= age, y = ITC)) + 
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + geom_point()+ stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) 


# Finn Code ----

tfi40_long <- pivot_longer(TFI_40, cols = Total:Evoked, names_to = 'measure', values_to = 'value')
tfi40_summary <- tfi40_long %>% group_by(freqs, measure) %>%
  summarize(y = mean(value, na.rm=T),
            sd = sd(value, na.rm=T),
            n = n(),
            se = sd / sqrt(n))

ggplot(data = tfi40_summary %>% filter(freqs <= 50, measure %in% c('Evoked','Induced')), aes(x = freqs, y = y, color = measure, fill = measure)) +
  geom_line() +
  geom_ribbon(aes(ymin = y - se, ymax = y + se), alpha = 0.4, linewidth = 0) +
  scale_y_log10() + 
  theme_bw()

snr <- TFI_40 %>% mutate(rawSNR = Evoked / Induced) %>% group_by(freqs) %>%
  summarize(SNR = mean(rawSNR, na.rm=T),
            sd = sd(rawSNR, na.rm=T),
            n = n(),
            se = sd / sqrt(n))

ggplot(data = snr %>% filter(freqs <= 50), aes(x = freqs, y = SNR)) +
  geom_line() +
  geom_ribbon(aes(ymin = SNR - se, ymax = SNR + se), alpha = 0.4, linewidth = 0) +
  scale_y_log10() + 
  theme_bw()
