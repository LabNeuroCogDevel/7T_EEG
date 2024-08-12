
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
library(eegUtils) #remotes::install_github("craddm/eegUtils")
library(lme4)
library("ggpubr")

createDataframes <- function() {


agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')

chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
chanLocs$Channel <- chanLocs$labels

criticality <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Criticality/allDFA.csv')
criticality = subset(criticality, select = -c(6:16) )

colnames(criticality)[colnames(criticality) == "lunaid"] <- "luna"
colnames(criticality)[colnames(criticality) == "ld8"] <- "Subject"
colnames(criticality)[colnames(criticality) == "chan"] <- "Channel"


criticality_ages <- merge(agefile, criticality, by = "Subject") %>% merge(.,chanLocs, by = "Channel")


regionDF <- criticality_ages %>% filter(., Channel == 'F3' | Channel == 'F5'| Channel == 'F7') %>% group_by(Subject, eyes) %>% 
  summarise(fEIbeta = mean(fEIbeta, na.rm=T), 
            Ampbeta = mean(Ampbeta, na.rm=T), 
            DFAbeta = mean(DFAbeta, na.rm=T)) %>% mutate(Region = "LDLPFC") %>% 
  rbind(.,criticality_ages %>% filter(., Channel == 'F4' | Channel == 'F6'| Channel == 'F8') %>% group_by(Subject, eyes) %>% 
          summarise(fEIbeta = mean(fEIbeta, na.rm=T), 
                    Ampbeta = mean(Ampbeta, na.rm=T), 
                    DFAbeta = mean(DFAbeta, na.rm=T)) %>% mutate(Region = "RDLPFC")) %>%
  rbind(.,criticality_ages %>% filter(., str_detect(Channel, "F")) %>% group_by(Subject, eyes) %>% 
          summarise(fEIbeta = mean(fEIbeta, na.rm=T), 
                    Ampbeta = mean(Ampbeta, na.rm=T), 
                    DFAbeta = mean(DFAbeta, na.rm=T)) %>% mutate(Region = "MPFC")) %>%
  merge(., agefile, by="Subject")  %>% separate(Subject,c("luna","vdate"), remove=F)

write.csv(regionDF, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Criticality/RegionalCriticalityMeasures.csv')

}

criticalityMRS <- function() {
  
  MRS <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj.csv')
  criticalityDF <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Criticality/RegionalCriticalityMeasures.csv')
  
  MRS_criticality <- merge(MRS, criticalityDF, by = c("luna", "visitno", "Region"))
  
  lunaize(ggplot(data = MRS_criticality , aes(x = Glu.Cr.adj, y = DFAbeta, by =luna, color =age.x)) + geom_point() + geom_line(aes(group=luna))+ geom_smooth(aes(group = 1), method = "lm")) + facet_wrap(.~eyes+Region)+ scale_color_gradient2(midpoint=20, low="gold3", mid="white", high="blue4")
  
  gam.model <-  gam(Spontaneous ~ s(Offset, k=3)  + s(age.x, k = 3) + s(age.x, k=3, by=Offset) + Condition + s(age.x, k=3, by=Condition), data = MRS_criticality  %>% mutate(ageGroup = as.factor(ifelse(age.x <= 18, 'Adol', 'Adult'))) , random = ~(1|luna))
  
  
  
  
}