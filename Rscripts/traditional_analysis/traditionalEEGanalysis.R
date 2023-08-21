
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
library(interactions)

traditionalPower_wholeBrain <- function() {

  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  traditionalGammaPowerDelay <- traditionalEEG()
  
(lunaize(ggplot(data = traditionalGammaPowerDelay, aes(x = age, y = (logPower))) + geom_point() + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Whole Brain Gamma Power vs Age") + xlab("Age") + ylab("Gamma Power")) +theme(plot.title = element_text(hjust = 0.5)))
  
lm.model <- lm(data = traditionalGammaPowerDelay[], ((logPower)) ~ (inverseAge))
anova(lm.model) 
summary(lm.model)

}


traditionalPower_allChannels <- function() {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  traditionalGammaPower_allChannels <- traditionalEEG_IndividualChannels()
  traditionalGammaPower_allChannels <- merge(traditionalGammaPower_allChannels, Behavior, by = c("Subject", "age", "inverseAge"))
  traditionalGammaPower_allChannels$Epoch <- as.factor(traditionalGammaPower_allChannels$Epoch)
  

  wholeBrain <- aggregate(.~Subject + Epoch, data = traditionalGammaPower_allChannels, mean)

  (lunaize(ggplot(data = wholeBrain2, aes(x = age, y = (logPower), color = Epoch)) + geom_point() + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Mean Channel Power vs Age") + xlab("Age") + ylab("Gamma Power")) +theme(plot.title = element_text(hjust = 0.5)))
  
  lm.model <- lm(data = wholeBrain[], ((logPower)) ~ (inverseAge))
  anova(lm.model) 
  summary(lm.model)
  
  
  Gamma_Frontal <- aggregate(.~Subject + Epoch, data = filter(traditionalGammaPower_allChannels, str_detect(traditionalGammaPower_allChannels$Channel, "F")), mean)
  
  (lunaize(ggplot(data = Gamma_Frontal, aes(x = age, y = (logPower))) + geom_point() + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Frontal Gamma Power vs Age") + xlab("Age") + ylab("Gamma Power")) +theme(plot.title = element_text(hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Frontal[], ((logPower)) ~ (inverseAge))
  anova(lm.model) 
  summary(lm.model)
  
  
  Gamma_Parietal <- aggregate(.~Subject + Epoch, data = filter(traditionalGammaPower_allChannels, str_detect(traditionalGammaPower_allChannels$Channel, "P")), mean)
  
  (lunaize(ggplot(data = Gamma_Parietal, aes(x = age, y = (logPower))) + geom_point() + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Parietal Gamma Power vs Age") + xlab("Age") + ylab("Gamma Power")) +theme(plot.title = element_text(hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Parietal[], ((logPower)) ~ (inverseAge))
  anova(lm.model) 
  summary(lm.model)
  
  
  Gamma_Occipital <- aggregate(.~Subject + Epoch, data = filter(traditionalGammaPower_allChannels, str_detect(traditionalGammaPower_allChannels$Channel, "O")), mean)
  
  (lunaize(ggplot(data = Gamma_Occipital, aes(x = age, y = (logPower))) + geom_point() + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Occipital Gamma Power vs Age") + xlab("Age") + ylab("Gamma Power")) +theme(plot.title = element_text(hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Occipital[], ((logPower)) ~ (inverseAge))
  anova(lm.model) 
  summary(lm.model)
  
  
  
  gamma_long <- rbind(
    wholeBrain %>% dplyr::select(Subject, inverseAge, age, Epoch, Power,mgsLatency, mgsLatency_sd, vgsLatency, vgsLatency_sd, absBestError, absBestError_sd) %>% mutate(type='Whole Brain'),
    Gamma_Frontal %>% dplyr::select(Subject, inverseAge, age, Epoch,  Power,  vgsLatency, vgsLatency_sd, mgsLatency, mgsLatency_sd,  absBestError, absBestError_sd) %>% mutate(type='Frontal')
  ) %>% rbind(., Gamma_Parietal %>% dplyr::select(Subject, inverseAge, age, Epoch,  Power,   vgsLatency, vgsLatency_sd, mgsLatency, mgsLatency_sd,  absBestError, absBestError_sd) %>% mutate(type='Parietal')) %>%  rbind(., Gamma_Occipital %>% dplyr::select(Subject, inverseAge, age, Epoch,  Power,  vgsLatency, vgsLatency_sd, mgsLatency, mgsLatency_sd,  absBestError, absBestError_sd) %>% mutate(type='Occipital'))
  
  gamma_long$logPower <- 20*log10(gamma_long$Power)
  
  
  lunaize(ggplot(data = gamma_long, 
                 aes(x = age, y = Power, color = Epoch)) + 
            stat_smooth(method='lm', alpha=0.15, formula = 'y~I(1/x)')) + facet_wrap( . ~ type)  + xlab("Age") + ylab("Gamma Power")
  
  lm.model <- lm(data = filter(gamma_long, str_detect(gamma_long$type, "Whole Brain")), ((logPower)) ~ (inverseAge*Epoch))
  anova(lm.model) 

  
  # Vs Behavior
  
  lunaize(ggplot(data = gamma_long, 
                 aes(x = logPower, y = vgsLatency_sd, color = Epoch)) + 
            geom_point() + stat_smooth(method='lm', alpha=0.15)) + facet_wrap( . ~ type)  + xlab("Gamma Power") + ylab("vgs Latency var (s)")
  
  
  lm.model <- lm(data = filter(gamma_long, str_detect(gamma_long$type, "Parietal")), ((vgsLatency_sd)) ~ (logPower))
  anova(lm.model) 
  summary(lm.model)
  
  lm.model <- lm(data = filter(gamma_long, str_detect(gamma_long$type, "Parietal")), ((vgsLatency_sd)) ~ (logPower*age))
  anova(lm.model) 
  summary(lm.model)
  
  
  lunaize(ggplot(data = gamma_long %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))), 
                 aes(x = logPower, y = vgsLatency_sd, color = ageGroup)) + 
            geom_point() + stat_smooth(method='lm') + xlab("Power") + ylab("vgs Latency var (s)")+ facet_wrap( . ~ type)) 
  
  
  gam.model <- gam(vgsLatency_sd ~ s(logPower, k=4, by=ageGroup) + s(age, k=4),
                    data = filter(gamma_long, str_detect(gamma_long$type, "Parietal")) %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult'))))
  summary(gam.model)  

  
## try johnson neyman plots
  lm.model <- lm(data = filter(gamma_long, str_detect(gamma_long$type, "Parietal")), ((mgsLatency_sd)) ~ (logPower*age))
  johnson_neyman(lm.model, pred = logPower, modx = age, plot = TRUE, title = "MGS Latency Var ~ logPower*age : Parietal")
  
}

