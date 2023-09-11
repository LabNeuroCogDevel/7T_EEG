
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

# Load dataframes ----

fooofMRSbehavior <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectsDLPFCfooofMRSBehaviorMeasures_20230822.csv")
fooofLong <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/Results/allSubjectsDLPFCfooofMeasures_20230613.csv")
MRSlong <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/fooof/Results/allSubjectsDLPFCMRSMeasures_20230613.csv")
fooofMRS<- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/allSubjectsDLPFCfooofMRSMeasures_20230822.csv")

# Exponent vs Age ----

lunaize(ggplot(data = fooofLong, 
               aes(x = age, y = Exponent)) + 
          geom_line(aes(group=interaction(luna,Region,Condition), linetype = Condition, shape =Region), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.2) + 
          geom_smooth(aes(group = 1, alpha = 0.5), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Exponent") + theme(legend.position = "none")


# Offset vs Age ----

lunaize(ggplot(data = fooofLong, 
               aes(x = age, y = Offset)) + 
          geom_line(aes(group=interaction(luna,Region,Condition), linetype = Condition, shape =Region), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.2) + 
          geom_smooth(aes(group = 1, alpha = 0.5), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Offset") + theme(legend.position = "none")

# Glu vs Age ----

lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = Glu_gamadj, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.2) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
          scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glutamate")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


# GABA vs Age ----

lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = GABA_gamadj, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.2) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
          scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("GABA")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

# Ratio vs Age ----

lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = Ratio_gamadj, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.2) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
          scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


# Imbalance vs Age ----

lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = GluGABAimbalanceABS, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.2) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
          scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


# Glu VS Exponent ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Exponent, y = Glu_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) +
          geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glutamate") + xlab("Exponent") + theme(text = element_text(size = 30))+ theme(legend.position='none')


# Glu VS Offset ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Offset, y = Glu_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) +
          geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glutamate") + xlab("Offset") + theme(text = element_text(size = 30))+ theme(legend.position='none')


# Imbalance VS Exponent ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Exponent, y = GluGABAimbalanceABS, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) +
          geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glu GABA Imbalance") + xlab("Exponent") + theme(text = element_text(size = 30))+ theme(legend.position='none')


# Imbalance VS Offset ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Offset, y = GluGABAimbalanceABS, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) +
          geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glu GABA Imbalance") + xlab("Offset") + theme(text = element_text(size = 30))+ theme(legend.position='none')


### Latency vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior , 
               aes(x = Exponent, y = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Latency") + xlab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

###  Latency  vs Offset ----

lunaize(ggplot(data = fooofMRSbehavior, 
               aes(x = Offset, y = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Latency") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

### Spatial Span Max vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(x = Exponent, y = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Span Max") + xlab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')



### Spatial Span Max vs Offset ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(x = Offset, y = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Span Max") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
