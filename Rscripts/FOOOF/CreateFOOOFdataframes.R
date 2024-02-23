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

# Define helper functions ----

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

## Prep EEG ----

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
fooof <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "eeg.eyesOpen_LDLPFC_Offset","eeg.eyesClosed_LDLPFC_Offset","eeg.eyesOpen_RDLPFC_Offset","eeg.eyesClosed_RDLPFC_Offset","eeg.eyesClosed_LDLPFC_Exponent", "eeg.eyesOpen_LDLPFC_Exponent","eeg.eyesOpen_RDLPFC_Exponent","eeg.eyesClosed_RDLPFC_Exponent")]

fooofLong <- fooof %>% 
  select(matches('lunaid|visitno|eeg.age|(eeg).*[LR]DLPFC.*(Exponent|Offset)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.age','Region'),
              names_from=c('data','measure')) %>% 
  select(matches('lunaid|visitno|Region|eeg.age|(eeg).*(Exponent|Offset)')) %>%
  pivot_longer(cols=matches('eyes'),
               names_pattern='(.*).(eyesOpen|eyesClosed)_(.*)',
               names_to=c("data", "Condition","measure"))  %>% 
  filter(!is.na(value))  %>% 
  pivot_wider(id_cols=c('lunaid','visitno','Region','eeg.age','Condition'),
              names_from=c('measure'))

fooofLong <- fooofLong %>% 
  mutate(ageGroup = as.factor(ifelse(eeg.age <= 18, 'Adol', 'Adult')),
         lunaid = as.factor(lunaid),
         inverseAge = 1/eeg.age) 

colnames(fooofLong) <- c("luna","visitno", "Region", "age", "Condition", "Offset", "Exponent", "ageGroup", "inverseAge")


write.csv(fooofLong,"/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsDLPFCfooofMeasures_20230613.csv")





## Prep MRS ----
MRS7t <- merge7t[c("lunaid","visitno","eeg.age", "sipfc.RDLPFC_GABA_gamadj", "sipfc.RDLPFC_Glu_gamadj",  "sipfc.LDLPFC_GABA_gamadj", "sipfc.LDLPFC_Glu_gamadj")]


MRSlong <- MRS7t %>% 
  select(matches('lunaid|visitno|eeg.age|(sipfc).*[LR]DLPFC.*(gamadj)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.age','Region'),
              names_from=c('data','measure'))

colnames(MRSlong) <- c("luna","visitno","age","Region","GABA_gamadj","Glu_gamadj")

### Create MRSI derivatives ----

idx <- which(!is.na(MRSlong$Glu_gamadj) & !is.na(MRSlong$GABA_gamadj))
gabaglu.lm <- lm(Glu_gamadj ~ GABA_gamadj + Region, data = MRSlong[idx,])
MRSlong$GluGABAimbalance <- NA
MRSlong$GluGABAimbalanceABS <- NA
MRSlong[idx,]$GluGABAimbalanceABS <- abs(gabaglu.lm$residuals)
MRSlong[idx,]$GluGABAimbalance <- (gabaglu.lm$residuals)

MRSlong$Ratio_gamadj <- MRSlong$Glu_gamadj/MRSlong$GABA_gamadj 

write.csv(MRSlong,"/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/fooof/Results/allSubjectsDLPFCMRSMeasures_20240223.csv")


## Prep Behavior ----

behav <- merge7t[c("lunaid","visitno","eeg.age", "eeg.vgsLatency_DelayAll","eeg.BestError_DelayAll", "eeg.BestError_sd_DelayAll","eeg.mgsLatency_DelayAll", "eeg.mgsLatency_sd_DelayAll","cantab.ssp_max.span","cantab.ssp_nerrors", "cantab.ssp_ntrials")]

colnames(behav) <- c("luna","visitno","age","vgsLatency","absBestError","absBestError_sd","mgsLatency","mgsLatency_sd", "SSP_maxSpan", "SSP_nErrors", "SSP_nTrials")
write.csv(behav,"/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/fooof/Results/allSubjectsBehavior.csv")


## FOOOF all channels ----

allchannels <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/Results/allSubjectsFooofMeasures_20230516.csv')
allchannels <- subset(allchannels, select = -X)

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
subjectInfo <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "sex")]
subjectInfo$Subject <- paste(subjectInfo$lunaid, subjectInfo$eeg.date, sep = "_")
fooofAllChannels <- merge(allchannels, subjectInfo, 'Subject')
fooofAllChannels <- subset(fooofAllChannels, select = -c(lunaid, eeg.date))
write.csv(fooofAllChannels,'/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/Results/allSubjectsAllChannelsFooofMeasures_20230911.csv')

# Merge Data Frames ----

fooofLong$Region <- as.factor(fooofLong$Region)
MRSlong$Region <- as.factor(MRSlong$Region)

fooofLong$luna <- as.factor(fooofLong$luna)
MRSlong$luna <- as.factor(MRSlong$luna)

fooofMRS <- merge(fooofLong, MRSlong, by = c("luna", "visitno", "age", "Region"))

behav$luna <- as.factor(behav$luna)
fooofMRSbehavior <- merge(fooofMRS, behav, by = c("luna", "visitno", "age"))

write.csv(fooofMRSbehavior,"/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectsDLPFCfooofMRSBehaviorMeasures_20230822.csv")
write.csv(fooofMRS,"/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectsDLPFCfooofMRSMeasures_20230822.csv")



