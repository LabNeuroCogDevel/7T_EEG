

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


# Load in dataframes ----

alphaDLPFC <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/alphaPeaks_DLPFCs_avgChannels_OutlierDetection_20230710.csv')
agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')

alpha_age <- merge(alphaDLPFC, agefile, by = "Subject")
alpha_age <- alpha_age %>% separate(Subject,c("luna","vdate"), remove=F)


# Alpha Measures vs age ----
## Center Frequency ----

lunaize(ggplot(data = alpha_age, 
               aes(x = age, y = Center.Frequency)) + 
          geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  facet_wrap(~Condition) + theme(text = element_text(size = 30)) + xlab("Age") + ylab("Center Frequency") + theme(legend.position = "none")

gam.model <-  gamm(Center.Frequency ~ s(age, k = 3)  + Condition + Region, data = alpha_age, random=list(luna=~1))
summary(gam.model$gam)

## Bandwidth ----

lunaize(ggplot(data = alpha_age, 
               aes(x = age, y = Bandwidth)) + 
          geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  facet_wrap(~Condition) + theme(text = element_text(size = 30)) + xlab("Age") + ylab("Bandwidth") + theme(legend.position = "none")


## Power ----

lunaize(ggplot(data = alpha_age, 
               aes(x = age, y = Power)) + 
          geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  facet_wrap(~Condition) + theme(text = element_text(size = 30)) + xlab("Age") + ylab("Power") + theme(legend.position = "none")


