
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

# Define outlier functions ----
outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

# Load in Merge7T ----
merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')

snr <- merge7t[c("lunaid","eeg.date","visitno","eeg.age","eeg.snr.ERSP_40Hz_LDLPFC", "eeg.snr.ERSP_40Hz_RDLPFC", "eeg.snr.ITC_40Hz_LDLPFC", 
                 "eeg.snr.ITC_40Hz_RDLPFC", "eeg.snr.BaselinePower_40Hz_LDLPFC", "eeg.snr.BaselinePower_40Hz_RDLPFC", "eeg.snr.Induced_40Hz_LDLPFC", 
                 "eeg.snr.Induced_40Hz_RDLPFC", "eeg.snr.Evoked_40Hz_RDLPFC", "eeg.snr.Evoked_40Hz_LDLPFC")]

snrLong <- snr %>% 
  select(matches('lunaid|visitno|eeg.age|(eeg).snr.*(ERSP|Induced|BaselinePower|ITC|Evoked).40Hz.*[LR]DLPFC')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*)_(40Hz).([LR]DLPFC)',
               names_to=c("measure","Freq", "Region"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.age','Region','Freq'),
              names_from=c('measure')) %>% 
  select(matches('lunaid|visitno|Region|Freq|eeg.age|(eeg).*(ERSP|Induced|BaselinePower|ITC|Evoked)'))

colnames(snrLong) <- c("luna", "visitno", "age", "Region", "Freq", "ERSP", 
                       "ITC", "BaselinePower", "Induced", "Evoked")

# Outlier by age group ----
snrLong <- snrLong %>% mutate(ageGroup = cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))

snrLong <- snrLong %>% 
  group_by(ageGroup) %>% 
  mutate(across(c("ERSP", "ITC", "Induced", "BaselinePower", "Evoked"), naoutlier))%>% ungroup()

# DLPFC from Merge7T ----
## ERSP  ----

lunaize(ggplot(data = snrLong, 
               aes(x = age, y = ERSP)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(aes(shape=Region, color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("ERSP") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) +
  scale_shape_manual(values = c(16, 17), labels = c("L DLPFC", "R DLPFC"))  + 
  guides(shape = guide_legend(title = "Region", override.aes = list(size = 2,fill=NA,alpha=0.8), title.theme = element_text(size = 15)), linetype = guide_legend(title = "Condition",override.aes = list(linewidth = 1,color ="black", alpha=0.8), title.theme = element_text(size = 15)), color="none", size ="none")


gam.model <-  gamm(ERSP ~ s(age, k = 3)  + Region, data = snrLong, random=list(luna=~1))
summary(gam.model$gam)


## ITC  ----

lunaize(ggplot(data = snrLong, 
               aes(x = age, y = ITC)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(aes(shape=Region, color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("ITC") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) +
  scale_shape_manual(values = c(16, 17), labels = c("L DLPFC", "R DLPFC"))  + 
  guides(shape = guide_legend(title = "Region", override.aes = list(size = 2,fill=NA,alpha=0.8), title.theme = element_text(size = 15)), linetype = guide_legend(title = "Condition",override.aes = list(linewidth = 1,color ="black", alpha=0.8), title.theme = element_text(size = 15)), color="none", size ="none")


gam.model <-  gamm(ITC ~ s(age, k = 3)+ Region, data = snrLong, random=list(luna=~1))
summary(gam.model$gam)


## Induced  ----

lunaize(ggplot(data = snrLong, 
               aes(x = age, y = Induced)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(aes(shape=Region, color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("Induced") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) +
  scale_shape_manual(values = c(16, 17), labels = c("L DLPFC", "R DLPFC"))  + 
  guides(shape = guide_legend(title = "Region", override.aes = list(size = 2,fill=NA,alpha=0.8), title.theme = element_text(size = 15)), linetype = guide_legend(title = "Condition",override.aes = list(linewidth = 1,color ="black", alpha=0.8), title.theme = element_text(size = 15)), color="none", size ="none")


gam.model <-  gamm(Induced ~ s(age, k = 3)+ Region, data = snrLong, random=list(luna=~1))
summary(gam.model$gam)


## Baseline Power   ----

lunaize(ggplot(data = snrLong, 
               aes(x = age, y = BaselinePower)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(aes(shape=Region, color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("Baseline Power") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) +
  scale_shape_manual(values = c(16, 17), labels = c("L DLPFC", "R DLPFC"))  + 
  guides(shape = guide_legend(title = "Region", override.aes = list(size = 2,fill=NA,alpha=0.8), title.theme = element_text(size = 15)), color="none", size ="none")


gam.model <-  gamm(BaselinePower ~ s(age, k = 3)+ Region, data = snrLong, random=list(luna=~1))
summary(gam.model$gam)


## Evoked Power   ----

lunaize(ggplot(data = snrLong, 
               aes(x = age, y = Evoked)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(aes(shape=Region, color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("Evoked Power") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) +
  scale_shape_manual(values = c(16, 17), labels = c("L DLPFC", "R DLPFC"))  + 
  guides(shape = guide_legend(title = "Region", override.aes = list(size = 2,fill=NA,alpha=0.8), title.theme = element_text(size = 15)), color="none", size ="none")


gam.model <-  gamm(Evoked ~ s(age, k = 3)+ Region, data = snrLong, random=list(luna=~1))
summary(gam.model$gam)



## ERSP / Baseline ----

snrLong$ERSPbaselineRatio <- snrLong$ERSP / snrLong$BaselinePower

lunaize(ggplot(data = snrLong, 
               aes(x = age, y = ERSPbaselineRatio)) + 
          geom_line(aes(group=interaction(luna)), alpha = 0.2) + 
          geom_point(aes(color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("ERSP/Baseline") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) 

gam.model <-  gamm(ERSPbaselineRatio ~ s(age, k = 3), data = SNR_CP, random=list(luna=~1))
summary(gam.model$gam)

# Load in Central Parietal Channel Data ----
SNR_CP <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectCentralParietal_SNRMeasures_20231114.csv')

# Outlier by age group ----
SNR_CP <- SNR_CP %>% mutate(ageGroup = cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))

SNR_CP <- SNR_CP %>% 
  group_by(ageGroup) %>% 
  mutate(across(c("ERSP", "ITC", "Induced", "BaselinePower"), naoutlier))%>% ungroup()

## ERSP  ----

lunaize(ggplot(data = SNR_CP, 
               aes(x = age, y = ERSP)) + 
          geom_line(aes(group=interaction(luna)), alpha = 0.2) + 
          geom_point(aes(color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("ERSP") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) 

gam.model <-  gamm(ERSP ~ s(age, k = 3) + labels, data = SNR_CP, random=list(luna=~1))
summary(gam.model$gam)

## ITC  ----

lunaize(ggplot(data = SNR_CP, 
               aes(x = age, y = ITC)) + 
          geom_line(aes(group=interaction(luna)), alpha = 0.2) + 
          geom_point(aes(color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("ITC") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) 


gam.model <-  gamm(ITC ~ s(age, k = 3) + labels, data = SNR_CP, random=list(luna=~1))
summary(gam.model$gam)

## Induced  ----

lunaize(ggplot(data = SNR_CP, 
               aes(x = age, y = Induced)) + 
          geom_line(aes(group=interaction(luna)), alpha = 0.2) + 
          geom_point(aes(color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("Induced") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) 


gam.model <-  gamm(Induced ~ s(age, k = 3) + labels, data = SNR_CP, random=list(luna=~1))
summary(gam.model$gam)


## Baseline Power   ----

lunaize(ggplot(data = SNR_CP, 
               aes(x = age, y = BaselinePower)) + 
          geom_line(aes(group=interaction(luna)), alpha = 0.2) + 
          geom_point(aes(color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("Baseline Power") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) 


gam.model <-  gamm(BaselinePower ~ s(age, k = 3) + labels, data = SNR_CP, random=list(luna=~1))
summary(gam.model$gam)

## ERSP / Baseline ----

SNR_CP$ERSPbaselineRatio <- SNR_CP$ERSP / SNR_CP$BaselinePower


lunaize(ggplot(data = SNR_CP, 
               aes(x = age, y = ERSPbaselineRatio)) + 
          geom_line(aes(group=interaction(luna)), alpha = 0.2) + 
          geom_point(aes(color=ageGroup), size=2.5,alpha=.8) + 
          geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.4,size=2)) + xlab("Age") + ylab("ERSP/Baseline") +
  scale_color_manual(values=c("#E3B23C", "plum4", "#2E041F","black"))+
  theme(text = element_text(size = 25), legend.text=element_text(size=15),
        legend.position = c(0.8, 0.9), legend.text.align = 0,legend.box = "horizontal", legend.key.width = unit(2, "line"), aspect.ratio = 1,
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2), plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) 

gam.model <-  gamm(ERSPbaselineRatio ~ s(age, k = 3), data = SNR_CP, random=list(luna=~1))
summary(gam.model$gam)







