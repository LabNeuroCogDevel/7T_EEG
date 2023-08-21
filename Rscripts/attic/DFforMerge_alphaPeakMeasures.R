
# Libraries ----

library(LNCDR)
library(data.table)
require(dplyr)
library(factoextra)
library(tidyverse)
require(knitr)
library(ggplot2)
library(e1071)
library(caret)
library(sjPlot)
library(directlabels)
attach(mtcars)
library(grid)
library(gridExtra)
library(cowplot)
library(plotrix)
library(mgcv)
library(plotly)
library(readxl)
library(interactions)
library(eegUtils)
library(lme4)
library("ggpubr")
library(jtools)
library(mgcViz)
library(lubridate)
library(checkmate)
library(lmerTest)

## Define helper functions ----

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}

nsd <- function(x) {
  (abs(x - mean(x, na.rm=T)) / sd(x, na.rm=T))
}

res_with_age <- function(MRSI_input, this_roi, met_name) {
  # have columns we need
  checkmate::expect_subset(c("roi", "age", "dateNumeric", met_name), names(MRSI_input))
  
  if (length(this_roi) > 1) {
    model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5) + s(GMrat, k=3) + label')
  } else {
    model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5) + s(GMrat, k=3)')
  }
  #print(model)
  mrsi.gam <- gam(as.formula(model), 
                  data=MRSI_input %>% 
                    filter(roi %in% this_roi) %>% 
                    mutate(label = as.factor(label)), na.action = na.exclude)
  
  met_out <- MRSI_input %>% filter(roi %in% this_roi)
  met_out$met_adj <- predict(mrsi.gam, met_out %>% 
                               mutate(dateNumeric = mean(met_out$dateNumeric, na.rm=T), 
                                      GMrat = mean(met_out$GMrat, na.rm=T),
                                      label = as.factor(label))) +
    residuals(mrsi.gam)
  
  
  met_out$met_adj <- as.numeric(met_out$met_adj)
  return(met_out)
}
# Load in Data Frames ----

peakMeasures <- read.csv(hera('/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsPeakMeasures_20230710.csv'))
peakMeasures$urchan <- peakMeasures$Channel+1
peakMeasures <- subset(peakMeasures, select = c(-Channel))

chanLocs <- read.csv(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv'))
chanLocs$Channel <- chanLocs$labels

peak_channels <- merge(peakMeasures, chanLocs, by = c("urchan"))

# Select Alpha Peaks ----

# select for peaks whose center frequency is in the alpha range (8-12 Hz)

alphaPeaks <- peak_channels %>% filter(Center.Frequency > 8, Center.Frequency < 13)


# Define DLPFC Channels ----

RDLPFC <- filter(alphaPeaks, alphaPeaks$labels == 'F4' | alphaPeaks$labels == 'F6'| alphaPeaks$labels == 'F8')
RDLPFC$Region <- 'R DLPFC'

LDLPFC <- filter(alphaPeaks, alphaPeaks$labels == 'F3' | alphaPeaks$labels == 'F5'| alphaPeaks$labels == 'F7')
LDLPFC$Region <- 'L DLPFC'

bothDLPFCs <- rbind(RDLPFC, LDLPFC)

write.csv(bothDLPFCs, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/alphaPeaks_DLPFCs_20230710.csv')


# Outlier Detection ----

bothDLPFCs_new <- bothDLPFCs %>% group_by(Subject, Condition) %>% mutate(Center.Frequency = ifelse(!outliers(Center.Frequency), Center.Frequency, NA)) %>% ungroup

bothDLPFCs_new <- bothDLPFCs_new %>% group_by(Subject, Condition) %>% mutate(Bandwidth = ifelse(!outliers(Bandwidth), Bandwidth, NA)) %>% ungroup

bothDLPFCs_new <- bothDLPFCs_new %>% group_by(Subject, Condition) %>% mutate(Power = ifelse(!outliers(Power), Power, NA)) %>% ungroup

# Average Trials ----

agg_fml <- cbind(Center.Frequency, Bandwidth, Power) ~ Subject + Region + Condition
aggregateTrials    <- aggregate(bothDLPFCs_new, agg_fml, function(x) mean(x,na.rm=T), na.action=NULL)
aggregateTrials_sd <- aggregate(bothDLPFCs_new, agg_fml, function(x) sd(x,na.rm=T)  , na.action=NULL)

regionLevel <- merge(aggregateTrials, aggregateTrials_sd, by = c("Subject", "Region", "Condition"), suffixes = c("", "_sd")) 

# Subject level outlier detection ----

regionLevel_new <- regionLevel

outlier_cols <- names(regionLevel[c(4:9)])
na_outliers <- function(x) ifelse(!outliers(x), x, NA)
regionLevel_new <- regionLevel %>% group_by(Region, Condition) %>% mutate(across({outlier_cols}, na_outliers)) %>% ungroup 
regionLevel_new$Region <- as.factor(regionLevel_new$Region)

write.csv(regionLevel_new, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/alphaPeaks_DLPFCs_avgChannels_OutlierDetection_20230710.csv')

