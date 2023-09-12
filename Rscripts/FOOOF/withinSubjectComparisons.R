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

# Load FOOOF ----

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
fooof <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "eeg.eyesOpen_LDLPFC_Offset","eeg.eyesClosed_LDLPFC_Offset","eeg.eyesOpen_RDLPFC_Offset","eeg.eyesClosed_RDLPFC_Offset","eeg.eyesClosed_LDLPFC_Exponent", "eeg.eyesOpen_LDLPFC_Exponent","eeg.eyesOpen_RDLPFC_Exponent","eeg.eyesClosed_RDLPFC_Exponent")]
fooofLong <- fooof %>% 
  select(matches('lunaid|eeg.date|visitno|eeg.age|(eeg).*[LR]DLPFC.*(Exponent|Offset)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.date','eeg.age','Region'),
              names_from=c('data','measure')) %>% 
  select(matches('lunaid|eeg.date|visitno|Region|eeg.age|(eeg).*(Exponent|Offset)')) %>%
  pivot_longer(cols=matches('eyes'),
               names_pattern='(.*).(eyesOpen|eyesClosed)_(.*)',
               names_to=c("data", "Condition","measure"))  %>% 
  filter(!is.na(value))  %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.date','Region','eeg.age','Condition'),
              names_from=c('measure'))

fooofLong <- fooofLong %>% 
  mutate(ageGroup = as.factor(ifelse(eeg.age <= 18, 'Adol', 'Adult')),
         lunaid = as.factor(lunaid),
         inverseAge = 1/eeg.age) 

colnames(fooofLong) <- c("lunaid","visitno","eeg.date" ,"Region", "age", "Condition", "Offset", "Exponent", "ageGroup", "inverseAge")

fooofLong <- fooofLong %>%
  mutate(combination_id = paste(lunaid, Region, Condition))

# Pivot the data wider with visitno as columns
fooof_wider <- fooofLong %>%
  pivot_wider(
    id_cols = c("combination_id"),
    names_from = visitno,
    values_from = c("Offset", "Exponent", "age","eeg.date")
  )

fooof_wider <- fooof_wider %>%
  separate(combination_id, into = c("lunaid", "Region", "Condition"), sep = " ")


# Format Dates

# Convert the date strings to Date objects
fooof_wider$eeg.date_1 <- as.Date(as.character(fooof_wider$eeg.date_1), format = "%Y%m%d")
fooof_wider$eeg.date_2 <- as.Date(as.character(fooof_wider$eeg.date_2), format = "%Y%m%d")
fooof_wider$eeg.date_3 <- as.Date(as.character(fooof_wider$eeg.date_3), format = "%Y%m%d")

# Calculate the number of days between the two dates
fooof_wider$v2_v1_days <- (as.numeric(fooof_wider$eeg.date_2 - fooof_wider$eeg.date_1))/365
fooof_wider$v3_v2_days <- (as.numeric(fooof_wider$eeg.date_3 - fooof_wider$eeg.date_2))/365

## Calculate within Sub FOOOF measures 
fooof_wider$newExp_1_2 <- (as.numeric(fooof_wider$Exponent_2 - fooof_wider$Exponent_1))/fooof_wider$v2_v1_days
fooof_wider$newExp_2_3 <- (as.numeric(fooof_wider$Exponent_3 - fooof_wider$Exponent_2))/fooof_wider$v3_v2_days

fooof_wider$newOff_1_2 <- (as.numeric(fooof_wider$Offset_2 - fooof_wider$Offset_1))/fooof_wider$v2_v1_days
fooof_wider$newOff_2_3 <- (as.numeric(fooof_wider$Offset_3 - fooof_wider$Offset_2))/fooof_wider$v3_v2_days

fooof_wider$newAge_1_2 <- (as.numeric(fooof_wider$age_2 + fooof_wider$age_1)/2)
fooof_wider$newAge_2_3 <- (as.numeric(fooof_wider$age_3 + fooof_wider$age_2)/2)

# pivot back to longer to have a col with 1-2 and 2-3

foooflogitud <- fooof_wider[c("lunaid","Region","Condition","newExp_1_2","newExp_2_3","newOff_1_2","newOff_2_3","newAge_1_2","newAge_2_3")]

fooof_longer <- foooflogitud %>%
  pivot_longer(
    cols = starts_with("new"),
    names_to = c(".value", "visits"),
    names_pattern = "([^_]+)_(\\d_\\d)",
    values_to = "value"
  )

fooof_longer <- fooof_longer %>%
  mutate(visits = sub("new", "", visits))

# Graph FOOOF Rate of Change ----
## Exponent rate of change vs age ----
lunaize(ggplot(data = fooof_longer, 
                      aes(x = newAge, y = newExp)) + 
                 geom_line(aes(group=interaction(lunaid,Region,Condition), linetype = Condition, shape =Region), alpha =0.4) + 
                 geom_point(aes(shape=Region),alpha=.8) + 
                 geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Age Diff") + ylab("Exponent Rate of Change")
  theme(text = element_text(size = 30), 
        legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
        axis.title.x = element_text(vjust = -2)) 

  
  ## Offset rate of change vs age ----
  lunaize(ggplot(data = fooof_longer, 
                 aes(x = newAge, y = newOff)) + 
            geom_line(aes(group=interaction(lunaid,Region,Condition), linetype = Condition, shape =Region), alpha = 0.4) + 
            geom_point(aes(shape=Region),alpha=.8) + 
            geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Age Diff") + ylab("Offset Rate of Change")
  theme(text = element_text(size = 30), 
        legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
        axis.title.x = element_text(vjust = -2)) 
  
  
  ## Offset rate of change vs Exponent rate of change ----
  lunaize(ggplot(data = fooof_longer, 
                 aes(x = newOff, y = newExp)) + 
            geom_line(aes(group=interaction(lunaid,Region,Condition), linetype = Condition, shape =Region), alpha = 0.4) + 
            geom_point(aes(shape=Region, color = newAge),alpha=.8) + 
            geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Offset Rate of Change") + ylab("Exponent Rate of Change")
  theme(text = element_text(size = 30), 
        legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
        axis.title.x = element_text(vjust = -2)) 

  ## Exponent 1 vs exponent rate of change ----
  lunaize(ggplot(data = fooof_longer, 
                 aes(x = Exponent_1, y = newExp_1_2)) + 
            geom_line(aes(group=interaction(lunaid,Region,Condition), linetype = Condition, shape =Region, color = newAge_1_2), alpha = 0.4) + 
            geom_point(aes(shape=Region, color = newAge_1_2),alpha=.8) + 
            geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Exponent 1") + ylab("Exponent Rate of Change V2-V1")
  theme(text = element_text(size = 30), 
        legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
        axis.title.x = element_text(vjust = -2)) 

# Load MRS ----

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
MRS7t <- merge7t[c("lunaid","visitno","rest.date","rest.age", "sipfc.RDLPFC_GABA_gamadj", "sipfc.RDLPFC_Glu_gamadj",  "sipfc.LDLPFC_GABA_gamadj", "sipfc.LDLPFC_Glu_gamadj")]

MRSlong <- MRS7t %>% 
  select(matches('lunaid|visitno|rest.date|rest.age|(sipfc).*[LR]DLPFC.*(gamadj)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','rest.date','rest.age','Region'),
              names_from=c('data','measure'))

colnames(MRSlong) <- c("lunaid","visitno","rest.date","age","Region","GABA_gamadj","Glu_gamadj")

idx <- which(!is.na(MRSlong$Glu_gamadj) & !is.na(MRSlong$GABA_gamadj))
gabaglu.lm <- lm(Glu_gamadj ~ GABA_gamadj + Region, data = MRSlong[idx,])
MRSlong$GluGABAimbalance <- NA
MRSlong$GluGABAimbalanceABS <- NA
MRSlong[idx,]$GluGABAimbalanceABS <- abs(gabaglu.lm$residuals)
MRSlong[idx,]$GluGABAimbalance <- (gabaglu.lm$residuals)

MRSlong$Ratio_gamadj <- MRSlong$Glu_gamadj/MRSlong$GABA_gamadj 

MRSlong <- MRSlong %>%
  mutate(combination_id = paste(lunaid, Region))

# Pivot the data wider with visitno as columns
MRS_wider <- MRSlong %>%
  pivot_wider(
    id_cols = c("combination_id"),
    names_from = visitno,
    values_from = c("GABA_gamadj", "Glu_gamadj","Ratio_gamadj", "GluGABAimbalanceABS","age","rest.date")
  )

MRS_wider <- MRS_wider %>%
  separate(combination_id, into = c("lunaid", "Region"), sep = " ")



# Format Dates

# Convert the date strings to Date objects
MRS_wider$rest.date_1 <- as.Date(as.character(MRS_wider$rest.date_1), format = "%Y%m%d")
MRS_wider$rest.date_2 <- as.Date(as.character(MRS_wider$rest.date_2), format = "%Y%m%d")
MRS_wider$rest.date_3 <- as.Date(as.character(MRS_wider$rest.date_3), format = "%Y%m%d")

# Calculate the number of days between the two dates
MRS_wider$v2_v1_days <- (as.numeric(MRS_wider$rest.date_2 - MRS_wider$rest.date_1))/365
MRS_wider$v3_v2_days <- (as.numeric(MRS_wider$rest.date_3 - MRS_wider$rest.date_2))/365

## Calculate within Sub MRS measures 
MRS_wider$newGABA_1_2 <- (as.numeric(MRS_wider$GABA_gamadj_2 - MRS_wider$GABA_gamadj_1))/MRS_wider$v2_v1_days
MRS_wider$newGABA_2_3 <- (as.numeric(MRS_wider$GABA_gamadj_3 - MRS_wider$GABA_gamadj_2))/MRS_wider$v3_v2_days

MRS_wider$newRatio_1_2 <- (as.numeric(MRS_wider$Ratio_gamadj_2 - MRS_wider$Ratio_gamadj_1))/MRS_wider$v2_v1_days
MRS_wider$newRatio_2_3 <- (as.numeric(MRS_wider$Ratio_gamadj_3 - MRS_wider$Ratio_gamadj_2))/MRS_wider$v3_v2_days

MRS_wider$newImb_1_2 <- (as.numeric(MRS_wider$GluGABAimbalanceABS_2 - MRS_wider$GluGABAimbalanceABS_1))/MRS_wider$v2_v1_days
MRS_wider$newImb_2_3 <- (as.numeric(MRS_wider$GluGABAimbalanceABS_3 - MRS_wider$GluGABAimbalanceABS_2))/MRS_wider$v3_v2_days

MRS_wider$newGlu_1_2 <- (as.numeric(MRS_wider$Glu_gamadj_2 - MRS_wider$Glu_gamadj_1))/MRS_wider$v2_v1_days
MRS_wider$newGlu_2_3 <- (as.numeric(MRS_wider$Glu_gamadj_3 - MRS_wider$Glu_gamadj_2))/MRS_wider$v3_v2_days

MRS_wider$newAge_1_2 <- (as.numeric(MRS_wider$age_2 + MRS_wider$age_1)/2)
MRS_wider$newAge_2_3 <- (as.numeric(MRS_wider$age_3 + MRS_wider$age_2)/2)

# pivot back to longer to have a col with 1-2 and 2-3

MRSlogitud <- MRS_wider[c("lunaid","Region","newGlu_1_2","newGlu_2_3","newGABA_1_2","newGABA_2_3","newRatio_1_2", "newRatio_2_3","newImb_1_2","newImb_2_3","newAge_1_2","newAge_2_3")]

MRS_longer <- MRSlogitud %>%
  pivot_longer(
    cols = starts_with("new"),
    names_to = c(".value", "visits"),
    names_pattern = "([^_]+)_(\\d_\\d)",
    values_to = "value"
  )

MRS_longer <- MRS_longer %>%
  mutate(visits = sub("new", "", visits))


# Graph MRS Rate of Change ----
## Glu rate of change vs age ----
lunaize(ggplot(data = MRS_longer, 
               aes(x = newAge, y = newGlu)) + 
          geom_line(aes(group=interaction(lunaid, Region),  shape =Region), alpha = 0.4) + 
          geom_point(aes(shape=Region),alpha=.8) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Avg Age") + ylab("Glu Rate of Change")
theme(text = element_text(size = 30), 
      legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
      axis.title.x = element_text(vjust = -2)) 


## GABA rate of change vs age ----
lunaize(ggplot(data = MRS_longer, 
               aes(x = newAge, y = newGABA)) + 
          geom_line(aes(group=interaction(lunaid,Region),  shape =Region), alpha = 0.4) + 
          geom_point(aes(shape=Region),alpha=.8) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Avg Age") + ylab("GABA Rate of Change")
theme(text = element_text(size = 30), 
      legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
      axis.title.x = element_text(vjust = -2)) 

## Ratio rate of change vs age ----
lunaize(ggplot(data = MRS_longer, 
               aes(x = newAge, y = newRatio)) + 
          geom_line(aes(group=interaction(lunaid,Region),  shape =Region), alpha = 0.4) + 
          geom_point(aes(shape=Region),alpha=.8) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Avg Age") + ylab("Ratio Rate of Change")
theme(text = element_text(size = 30), 
      legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
      axis.title.x = element_text(vjust = -2)) 

## Imbalance rate of change vs age ----
lunaize(ggplot(data = MRS_longer, 
               aes(x = newAge, y = newImb)) + 
          geom_line(aes(group=interaction(lunaid,Region),  shape =Region), alpha = 0.4) + 
          geom_point(aes(shape=Region),alpha=.8) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1) + geom_abline(slope=0,intercept=0, linetype=2)) + xlab("Avg Age") + ylab("Imbalance Rate of Change")
theme(text = element_text(size = 30), 
      legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
      axis.title.x = element_text(vjust = -2)) 


# Combine FOOOF and MRS ----

fooofMRSwider <- merge(MRS_longer, fooof_longer, by = c("lunaid","Region","visits"), suffixes = c("_MRS", "_fooof"))

## Exponent Rate of Change vs Imbalance Rate of Change ----
lunaize(ggplot(data = fooofMRSwider, 
               aes(y = newExp, x = newImb)) + 
          geom_line(aes(group=interaction(lunaid,Region),  shape =Region, linetype=Region), alpha = 0.4) + 
          geom_point(aes(shape=Region),alpha=.8) + 
          geom_smooth(aes(group = 1), method="lm",alpha=0.8,size=1) +   facet_wrap(~Condition)+
          geom_abline(slope=0,intercept=0, linetype=2)) + ylab("Exponent Rate of Change") + xlab("Imbalance Rate of Change")
theme(text = element_text(size = 30), 
      legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
      axis.title.x = element_text(vjust = -2))


## Offset Rate of Change vs Glu Rate of Change ----
lunaize(ggplot(data = fooofMRSwider, 
               aes(y = newExp, x = newGlu)) + 
          geom_line(aes(group=interaction(lunaid,Region),  shape =Region, linetype=Region), alpha = 0.4) + 
          geom_point(aes(shape=Region),alpha=.8) + 
          geom_smooth(aes(group = 1), method="lm",alpha=0.8,size=1) +   facet_wrap(~Condition)+
          geom_abline(slope=0,intercept=0, linetype=2)) + ylab("Offset Rate of Change") + xlab("Glu Rate of Change")
theme(text = element_text(size = 30), 
      legend.position = "none", plot.margin = margin(0.5, 0.5, .5, 0.5, "cm") , aspect.ratio = 1, axis.title.y = element_text(vjust = 1.5), 
      axis.title.x = element_text(vjust = -2)) 