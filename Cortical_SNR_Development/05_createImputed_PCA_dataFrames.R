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
library(mice)
library(tidyr)


# Define Outlier Function ----
outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


# 40-40 ----
## SNR ----

SNRallChans_outlier_new <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/allSubjectsSNR_allChans_allfreqs.csv')
chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/resources/ChannelLocs.csv')%>% select(-type)

SNRallChans_outlier4040 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 40)
SNRallChans_outlier4040$SNR <- log(SNRallChans_outlier4040$Evoked/SNRallChans_outlier4040$Induced)


# Identify missing values
missing_values <- is.na(SNRallChans_outlier4040$Total)

unique_subjects <- unique(SNRallChans_outlier4040$lunaID)
unique_visits <- unique(SNRallChans_outlier4040$visitno)
imputeData <- SNRallChans_outlier4040


for (lunaID in unique_subjects) {
  for (visit in unique_visits) {
    if (any(imputeData$lunaID == lunaID & imputeData$visitno == visit)) {
      
      subject_visit_data <-
        imputeData[imputeData$lunaID == lunaID &
                     imputeData$visitno == visit,]
      
      if (any(is.na(subject_visit_data$Total) | is.na(subject_visit_data$Induced) | is.na(subject_visit_data$Evoked) | is.na(subject_visit_data$SNR))) {
        
        miceData <- mice(select(subject_visit_data, c("Induced", "Evoked")), m = 20, method = 'mean')
        miceDataTotal <- mice(select(subject_visit_data, c("Total","SNR")), m = 20, method = 'mean')
        
        imputed_micedata <- cbind(subject_visit_data[, c("lunaID", "visitno")], complete(miceData), complete(miceDataTotal))
        
        imputeData[imputeData$lunaID %in% imputed_micedata$lunaID & imputeData$visitno %in% imputed_micedata$visitno, c("Total", "Induced", "Evoked", "SNR")] <- imputed_micedata[,c("Total", "Induced", "Evoked", "SNR")]
        
      }
    }
  }
}


snrImputed <- imputeData %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>% 
  filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2','Fz', 'AF5', 'AF6')) 


snr <- imputeData %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>% 
  filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2','Fz', 'AF5', 'AF6')) %>%
  pivot_wider(names_from = labels, values_from = SNR, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, -visitno)

completeidx <- which(complete.cases(snr))
snr <- snr[completeidx,] 

snr.pca <- prcomp(scale(snr))
summary(snr.pca)
snr.pca$rotation

lunaize(
  ggplot(data=melt(snr.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


#calculate total variance explained by each principal component
var_explained = snr.pca$sdev^2 / sum(snr.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

SNRchanLocs <- chanLocs 
SNRchanLocs$PC1 <- 0
SNRchanLocs$PC2 <- 0
SNRchanLocs$PC3 <- 0


SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC1")] <- -1*(snr.pca$rotation[,1])
SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC2")] <- -1*(snr.pca$rotation[,2])
SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC3")] <- -1*(snr.pca$rotation[,3])

SNRchanLocs$measure <- "SNR"

write.csv(SNRchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/SNRchanlocs.csv')

snrAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = SNR, id_cols = c("lunaID", "age", "visitno"))

snrAll$pc1 <- NA
snrAll$pc1 <- -1*(unname(unlist(snr.pca$x[,1])))
snrAll$pc2 <- NA
snrAll$pc2 <- -1*unname(unlist(snr.pca$x[,2]))
snrAll$pc3 <- NA
snrAll$pc3 <- -1*unname(unlist(snr.pca$x[,3]))

snrAll <- snrAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Evoked ----

evokedImputed <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) 

evoked <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Evoked, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, -visitno)

completeidx <- which(complete.cases(evoked))
evoked <- evoked[completeidx,]

evoked.pca <- prcomp(scale(evoked))
summary(evoked.pca)
evoked.pca$rotation


lunaize(
  ggplot(data=melt(evoked.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


var_explained = evoked.pca$sdev^2 / sum(evoked.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

evokedchanLocs <- chanLocs 
evokedchanLocs$PC1 <- 0
evokedchanLocs$PC2 <- 0
evokedchanLocs$PC3 <- 0


evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC1")] <- -1*(evoked.pca$rotation[,1])
evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC2")] <- -1*(evoked.pca$rotation[,2])
evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC3")] <- -1*(evoked.pca$rotation[,3])

evokedchanLocs$measure <- "Evoked"

write.csv(evokedchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/evokedchanLocs.csv')


evokedAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Evoked, id_cols = c("lunaID", "age", "visitno"))

evokedAll$pc1 <- NA
evokedAll$pc1 <- -1*unname(unlist(evoked.pca$x[,1]))
evokedAll$pc2 <- NA
evokedAll$pc2 <- -1*unname(unlist(evoked.pca$x[,2]))
evokedAll$pc3 <- NA
evokedAll$pc3 <- -1*unname(unlist(evoked.pca$x[,3]))

evokedAll <- evokedAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Induced ----

inducedImputed <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Induced, labels, lunaID, age, visitno)

induced <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Induced, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Induced, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, - visitno)

completeidx <- which(complete.cases(induced))
induced <- induced[completeidx,]

induced.pca <- prcomp(scale(induced))
summary(induced.pca)
induced.pca$rotation

lunaize(
  ggplot(data=melt(induced.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


#calculate total variance explained by each principal component
var_explained = induced.pca$sdev^2 / sum(induced.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

inducedchanLocs <- chanLocs 
inducedchanLocs$PC1 <- 0
inducedchanLocs$PC2 <- 0
inducedchanLocs$PC3 <- 0


inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC1")] <- 1*(induced.pca$rotation[,1])
inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC2")] <- 1*(induced.pca$rotation[,2])
inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC3")] <- 1*(induced.pca$rotation[,3])

inducedchanLocs$measure <- "Induced"

write.csv(inducedchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/inducedchanLocs.csv')

inducedAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(Induced, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Induced, id_cols = c("lunaID", "visitno", "age")) 

inducedAll$pc1 <- NA
inducedAll$pc1 <- unname(unlist(induced.pca$x[,1]))
inducedAll$pc2 <- NA
inducedAll$pc2 <- unname(unlist(induced.pca$x[,2]))
inducedAll$pc3 <- NA
inducedAll$pc3 <- unname(unlist(induced.pca$x[,3]))

inducedAll <- inducedAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Merge all measures ----

allSNRmeasures_chanLocs <- rbind(inducedchanLocs, evokedchanLocs) %>% rbind(., SNRchanLocs)
write.csv(allSNRmeasures_chanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/snrMeasures_chanLocs.csv')

allSNRmeasures_imputed <- merge(snrImputed, evokedImputed, by = c("lunaID", "visitno", "age", "labels")) %>% merge(., inducedImputed, by = c("lunaID", "visitno", "age", "labels"))
write.csv(allSNRmeasures_imputed, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/snrImputed.csv')

snrAll$measure <- "SNR"
evokedAll$measure <- "Evoked"
inducedAll$measure <- "Induced"

allSNRmeasures4040 <- rbind(snrAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure), evokedAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure)) %>% 
  rbind(., inducedAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure))

write.csv(allSNRmeasures4040, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/SNRmeasures_PCAvalues.csv')



# 30-30 ----
## SNR ----
SNRallChans_outlier3030 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 30)
SNRallChans_outlier3030$SNR <- log(SNRallChans_outlier3030$Evoked/SNRallChans_outlier3030$Induced)


# Identify missing values
missing_values <- is.na(SNRallChans_outlier3030$Total)

unique_subjects <- unique(SNRallChans_outlier3030$lunaID)
unique_visits <- unique(SNRallChans_outlier3030$visitno)
imputeData <- SNRallChans_outlier3030


for (lunaID in unique_subjects) {
  for (visit in unique_visits) {
    if (any(imputeData$lunaID == lunaID & imputeData$visitno == visit)) {
      
      subject_visit_data <-
        imputeData[imputeData$lunaID == lunaID &
                     imputeData$visitno == visit,]
      
      if (any(is.na(subject_visit_data$Total) | is.na(subject_visit_data$Induced) | is.na(subject_visit_data$Evoked) | is.na(subject_visit_data$SNR))) {
        
        miceData <- mice(select(subject_visit_data, c("Induced", "Evoked")), m = 20, method = 'mean')
        miceDataTotal <- mice(select(subject_visit_data, c("Total","SNR")), m = 20, method = 'mean')
        
        imputed_micedata <- cbind(subject_visit_data[, c("lunaID", "visitno")], complete(miceData), complete(miceDataTotal))
        
        imputeData[imputeData$lunaID %in% imputed_micedata$lunaID & imputeData$visitno %in% imputed_micedata$visitno, c("Total", "Induced", "Evoked", "SNR")] <- imputed_micedata[,c("Total", "Induced", "Evoked", "SNR")]
        
      }
    }
  }
}


snrImputed <- imputeData %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>% 
  filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2','Fz', 'AF5', 'AF6')) 


snr <- imputeData %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>% 
  filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2','Fz', 'AF5', 'AF6')) %>%
  pivot_wider(names_from = labels, values_from = SNR, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, -visitno)

completeidx <- which(complete.cases(snr))
snr <- snr[completeidx,] 

snr.pca <- prcomp(scale(snr))
summary(snr.pca)
snr.pca$rotation

lunaize(
  ggplot(data=melt(snr.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


#calculate total variance explained by each principal component
var_explained = snr.pca$sdev^2 / sum(snr.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

SNRchanLocs <- chanLocs 
SNRchanLocs$PC1 <- 0
SNRchanLocs$PC2 <- 0
SNRchanLocs$PC3 <- 0


SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC1")] <- -1*(snr.pca$rotation[,1])
SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC2")] <- -1*(snr.pca$rotation[,2])
SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC3")] <- -1*(snr.pca$rotation[,3])

SNRchanLocs$measure <- "SNR"

write.csv(SNRchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/SNRchanlocs3030.csv')

snrAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = SNR, id_cols = c("lunaID", "age", "visitno"))

snrAll$pc1 <- NA
snrAll$pc1 <- -1*(unname(unlist(snr.pca$x[,1])))
snrAll$pc2 <- NA
snrAll$pc2 <- -1*unname(unlist(snr.pca$x[,2]))
snrAll$pc3 <- NA
snrAll$pc3 <- -1*unname(unlist(snr.pca$x[,3]))

snrAll <- snrAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Evoked ----

evokedImputed <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) 

evoked <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Evoked, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, -visitno)

completeidx <- which(complete.cases(evoked))
evoked <- evoked[completeidx,]

evoked.pca <- prcomp(scale(evoked))
summary(evoked.pca)
evoked.pca$rotation


lunaize(
  ggplot(data=melt(evoked.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


var_explained = evoked.pca$sdev^2 / sum(evoked.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

evokedchanLocs <- chanLocs 
evokedchanLocs$PC1 <- 0
evokedchanLocs$PC2 <- 0
evokedchanLocs$PC3 <- 0


evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC1")] <- -1*(evoked.pca$rotation[,1])
evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC2")] <- -1*(evoked.pca$rotation[,2])
evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC3")] <- -1*(evoked.pca$rotation[,3])

evokedchanLocs$measure <- "Evoked"

write.csv(evokedchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/evokedchanLocs3030.csv')


evokedAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Evoked, id_cols = c("lunaID", "age", "visitno"))

evokedAll$pc1 <- NA
evokedAll$pc1 <- -1*unname(unlist(evoked.pca$x[,1]))
evokedAll$pc2 <- NA
evokedAll$pc2 <- -1*unname(unlist(evoked.pca$x[,2]))
evokedAll$pc3 <- NA
evokedAll$pc3 <- -1*unname(unlist(evoked.pca$x[,3]))

evokedAll <- evokedAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Induced ----

inducedImputed <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Induced, labels, lunaID, age, visitno)

induced <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Induced, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Induced, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, - visitno)

completeidx <- which(complete.cases(induced))
induced <- induced[completeidx,]

induced.pca <- prcomp(scale(induced))
summary(induced.pca)
induced.pca$rotation

lunaize(
  ggplot(data=melt(induced.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


#calculate total variance explained by each principal component
var_explained = induced.pca$sdev^2 / sum(induced.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

inducedchanLocs <- chanLocs 
inducedchanLocs$PC1 <- 0
inducedchanLocs$PC2 <- 0
inducedchanLocs$PC3 <- 0


inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC1")] <- 1*(induced.pca$rotation[,1])
inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC2")] <- 1*(induced.pca$rotation[,2])
inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC3")] <- 1*(induced.pca$rotation[,3])

inducedchanLocs$measure <- "Induced"

write.csv(inducedchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/inducedchanLocs3030.csv')

inducedAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(Induced, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Induced, id_cols = c("lunaID", "visitno", "age")) 

inducedAll$pc1 <- NA
inducedAll$pc1 <- unname(unlist(induced.pca$x[,1]))
inducedAll$pc2 <- NA
inducedAll$pc2 <- unname(unlist(induced.pca$x[,2]))
inducedAll$pc3 <- NA
inducedAll$pc3 <- unname(unlist(induced.pca$x[,3]))

inducedAll <- inducedAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Merge all measures ----

allSNRmeasures_chanLocs <- rbind(inducedchanLocs, evokedchanLocs) %>% rbind(., SNRchanLocs)
write.csv(allSNRmeasures_chanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/snrMeasures_chanLocs3030.csv')

allSNRmeasures_imputed <- merge(snrImputed, evokedImputed, by = c("lunaID", "visitno", "age", "labels")) %>% merge(., inducedImputed, by = c("lunaID", "visitno", "age", "labels"))
write.csv(allSNRmeasures_imputed, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/snrImputed3030.csv')

snrAll$measure <- "SNR"
evokedAll$measure <- "Evoked"
inducedAll$measure <- "Induced"

allSNRmeasures3030 <- rbind(snrAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure), evokedAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure)) %>% 
  rbind(., inducedAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure))

write.csv(allSNRmeasures3030, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/SNRmeasures_PCAvalues3030.csv')



# 20-20 ----
## SNR ----

SNRallChans_outlier2020 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 20)
SNRallChans_outlier2020$SNR <- log(SNRallChans_outlier2020$Evoked/SNRallChans_outlier2020$Induced)


# Identify missing values
missing_values <- is.na(SNRallChans_outlier2020$Total)

unique_subjects <- unique(SNRallChans_outlier2020$lunaID)
unique_visits <- unique(SNRallChans_outlier2020$visitno)
imputeData <- SNRallChans_outlier2020


for (lunaID in unique_subjects) {
  for (visit in unique_visits) {
    if (any(imputeData$lunaID == lunaID & imputeData$visitno == visit)) {
      
      subject_visit_data <-
        imputeData[imputeData$lunaID == lunaID &
                     imputeData$visitno == visit,]
      
      if (any(is.na(subject_visit_data$Total) | is.na(subject_visit_data$Induced) | is.na(subject_visit_data$Evoked) | is.na(subject_visit_data$SNR))) {
        
        miceData <- mice(select(subject_visit_data, c("Induced", "Evoked")), m = 20, method = 'mean')
        miceDataTotal <- mice(select(subject_visit_data, c("Total","SNR")), m = 20, method = 'mean')
        
        imputed_micedata <- cbind(subject_visit_data[, c("lunaID", "visitno")], complete(miceData), complete(miceDataTotal))
        
        imputeData[imputeData$lunaID %in% imputed_micedata$lunaID & imputeData$visitno %in% imputed_micedata$visitno, c("Total", "Induced", "Evoked", "SNR")] <- imputed_micedata[,c("Total", "Induced", "Evoked", "SNR")]
        
      }
    }
  }
}


snrImputed <- imputeData %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>% 
  filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2','Fz', 'AF5', 'AF6')) 


snr <- imputeData %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>% 
  filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2','Fz', 'AF5', 'AF6')) %>%
  pivot_wider(names_from = labels, values_from = SNR, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, -visitno)

completeidx <- which(complete.cases(snr))
snr <- snr[completeidx,] 

snr.pca <- prcomp(scale(snr))
summary(snr.pca)
snr.pca$rotation

lunaize(
  ggplot(data=melt(snr.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


#calculate total variance explained by each principal component
var_explained = snr.pca$sdev^2 / sum(snr.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

SNRchanLocs <- chanLocs 
SNRchanLocs$PC1 <- 0
SNRchanLocs$PC2 <- 0
SNRchanLocs$PC3 <- 0


SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC1")] <- -1*(snr.pca$rotation[,1])
SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC2")] <- -1*(snr.pca$rotation[,2])
SNRchanLocs[SNRchanLocs$labels %in% frontalChannels, c("PC3")] <- -1*(snr.pca$rotation[,3])

SNRchanLocs$measure <- "SNR"

write.csv(SNRchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/SNRchanlocs2020.csv')

snrAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(SNR, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = SNR, id_cols = c("lunaID", "age", "visitno"))

snrAll$pc1 <- NA
snrAll$pc1 <- -1*(unname(unlist(snr.pca$x[,1])))
snrAll$pc2 <- NA
snrAll$pc2 <- -1*unname(unlist(snr.pca$x[,2]))
snrAll$pc3 <- NA
snrAll$pc3 <- -1*unname(unlist(snr.pca$x[,3]))

snrAll <- snrAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Evoked ----

evokedImputed <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) 

evoked <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Evoked, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, -visitno)

completeidx <- which(complete.cases(evoked))
evoked <- evoked[completeidx,]

evoked.pca <- prcomp(scale(evoked))
summary(evoked.pca)
evoked.pca$rotation


lunaize(
  ggplot(data=melt(evoked.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


var_explained = evoked.pca$sdev^2 / sum(evoked.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

evokedchanLocs <- chanLocs 
evokedchanLocs$PC1 <- 0
evokedchanLocs$PC2 <- 0
evokedchanLocs$PC3 <- 0


evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC1")] <- -1*(evoked.pca$rotation[,1])
evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC2")] <- -1*(evoked.pca$rotation[,2])
evokedchanLocs[evokedchanLocs$labels %in% frontalChannels, c("PC3")] <- -1*(evoked.pca$rotation[,3])

evokedchanLocs$measure <- "Evoked"

write.csv(evokedchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/evokedchanLocs2020.csv')


evokedAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(Evoked, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Evoked, id_cols = c("lunaID", "age", "visitno"))

evokedAll$pc1 <- NA
evokedAll$pc1 <- -1*unname(unlist(evoked.pca$x[,1]))
evokedAll$pc2 <- NA
evokedAll$pc2 <- -1*unname(unlist(evoked.pca$x[,2]))
evokedAll$pc3 <- NA
evokedAll$pc3 <- -1*unname(unlist(evoked.pca$x[,3]))

evokedAll <- evokedAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Induced ----

inducedImputed <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Induced, labels, lunaID, age, visitno)

induced <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6'))  %>%
  dplyr::select(Induced, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Induced, id_cols = c("lunaID", "visitno")) %>% select(-lunaID, - visitno)

completeidx <- which(complete.cases(induced))
induced <- induced[completeidx,]

induced.pca <- prcomp(scale(induced))
summary(induced.pca)
induced.pca$rotation

lunaize(
  ggplot(data=melt(induced.pca$rotation[,1:6])) + 
    geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge')) + 
  labs(fill='RL Parameters', y='Loading', x='')


#calculate total variance explained by each principal component
var_explained = induced.pca$sdev^2 / sum(induced.pca$sdev^2)


qplot(c(1:6), var_explained[1:6]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

frontalChannels <-  (c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) 

inducedchanLocs <- chanLocs 
inducedchanLocs$PC1 <- 0
inducedchanLocs$PC2 <- 0
inducedchanLocs$PC3 <- 0


inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC1")] <- 1*(induced.pca$rotation[,1])
inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC2")] <- 1*(induced.pca$rotation[,2])
inducedchanLocs[inducedchanLocs$labels %in% frontalChannels, c("PC3")] <- 1*(induced.pca$rotation[,3])

inducedchanLocs$measure <- "Induced"

write.csv(inducedchanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/inducedchanLocs2020.csv')

inducedAll <- imputeData %>% filter(labels %in% c('F3', 'F5', 'F7', 'F1', 'F2', 'F4', 'F6', 'F8', 'AFz', 'AF1', 'AF2', 'Fp1', 'Fp2', 'Fz', 'AF5', 'AF6')) %>%
  dplyr::select(Induced, labels, lunaID, age, visitno) %>%
  pivot_wider(names_from = labels, values_from = Induced, id_cols = c("lunaID", "visitno", "age")) 

inducedAll$pc1 <- NA
inducedAll$pc1 <- unname(unlist(induced.pca$x[,1]))
inducedAll$pc2 <- NA
inducedAll$pc2 <- unname(unlist(induced.pca$x[,2]))
inducedAll$pc3 <- NA
inducedAll$pc3 <- unname(unlist(induced.pca$x[,3]))

inducedAll <- inducedAll %>%
  mutate(across(c("pc1", "pc2", "pc3"), naoutlier)) %>% ungroup()


## Merge all measures ----

allSNRmeasures_chanLocs <- rbind(inducedchanLocs, evokedchanLocs) %>% rbind(., SNRchanLocs)
write.csv(allSNRmeasures_chanLocs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/snrMeasures_chanLocs2020.csv')

allSNRmeasures_imputed <- merge(snrImputed, evokedImputed, by = c("lunaID", "visitno", "age", "labels")) %>% merge(., inducedImputed, by = c("lunaID", "visitno", "age", "labels"))
write.csv(allSNRmeasures_imputed, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/snrImputed2020.csv')

snrAll$measure <- "SNR"
evokedAll$measure <- "Evoked"
inducedAll$measure <- "Induced"

allSNRmeasures2020 <- rbind(snrAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure), evokedAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure)) %>% 
  rbind(., inducedAll %>% select(lunaID, visitno, age, pc1, pc2, pc3, measure))

write.csv(allSNRmeasures2020, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/SNRmeasures_PCAvalues2020.csv')




# Merge all Stim ----
allSNRmeasures4040$Stim <- "40 Hz"
allSNRmeasures3030$Stim <- "30 Hz"

allSNRmeasures2020$Stim <- "20 Hz"

allSNRmeasuresAllStim <- rbind(allSNRmeasures4040, allSNRmeasures3030) %>% rbind(., allSNRmeasures2020) %>% select(-pc2, -pc3) 
write.csv(allSNRmeasuresAllStim, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/SNRmeasures_PC1_allStim.csv')

