#!/usr/bin/env Rscript

# Libraries ----

library(LNCDR)
library(data.table)
library(dplyr)
library(factoextra)
library(ggplot2)
library(e1071)
attach(mtcars)
library(grid)
library(gridExtra)
library(mgcv)
library(readxl)
library(lme4)
library(lubridate)
library(checkmate)
library(lmerTest)
library(tidyr)

# Load data sets ----

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/resources/ChannelLocs.csv')%>% select(-type)
merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")

SNRallChans40 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/allSubjectsSNR_allChans_40Hz.csv')
SNRallChans40 <- SNRallChans40 %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))
SNRallChans40$hertz <- 40

SNRallChans30 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/allSubjectsSNR_allChans_30Hz.csv')
SNRallChans30 <- SNRallChans30 %>% separate(Subject, c('lunaID','visitDate'))%>% merge(.,ageValues, by = c("lunaID", "visitDate"))
SNRallChans30$hertz <- 30

SNRallChans20 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/allSubjectsSNR_allChans_20Hz.csv')
SNRallChans20 <- SNRallChans20 %>% separate(Subject, c('lunaID','visitDate'))%>% merge(.,ageValues, by = c("lunaID", "visitDate"))
SNRallChans20$hertz <- 20

SNRallChans <- rbind(SNRallChans40, SNRallChans30) %>% rbind(., SNRallChans20)

# Define outlier functions ----
outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


# Merge SNR with ages and Channel Info ----
SNRallChans <- merge(SNRallChans, chanLocs, by = "urchan")

# Create age groups and remove DC comp
SNRallChans <- SNRallChans %>% 
  mutate(ageGroup = as.factor(cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30'))))  %>% filter(freqs > 0) %>% filter(Total < 1000)

# SNR Outlier Detection ----
SNRallChans_outlier <- SNRallChans %>% group_by(urchan, freqs, hertz) %>%
  mutate(across(c("Evoked", "Induced", "Total"), naoutlier)) %>% ungroup()

SNRallChans_outlier$SNR <- log(SNRallChans_outlier$Evoked/SNRallChans_outlier$Induced)

SNRallChans_outlier <- SNRallChans_outlier %>% group_by(urchan, freqs) %>% mutate(across(c("SNR"), naoutlier)) %>% ungroup() %>% mutate(SNR = exp(SNR))
write.csv(SNRallChans_outlier, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/allSubjectsSNR_allChans_allfreqs.csv')
