

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

CreateRegionsfromElectrodes <- function (){
 
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')
  chanLocs$Channel <- chanLocs$labels
  
  fooof <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsFooofMeasures_20230111.csv')
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
  }
  
  fooof_new <- fooof
  
  #outliers per channel
  cols = names(fooof[4:5])
  for (col in cols) {
    
    fooof_grouped <- fooof %>% group_by(Subject) %>% group_by(Condition)
    
    indx <- outliers(fooof_grouped[[col]])
    
    fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
    fooof_new[[col]] <- as.numeric(fooof_new[[col]])
    
  }  
  
  fooof_age <- merge(fooof_new, agefile, by = 'Subject')
  
  fooof_channels <- merge(fooof_age, chanLocs, by = c("Channel"))
  
  fooof_regional <- fooof_channels %>% filter(., Channel == 'F3' | Channel == 'F5'| Channel == 'F7') %>% group_by(Subject, Condition) %>% 
    summarise(Offset = mean(Offset, na.rm=T), 
              Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "LDLPFC") %>% 
    rbind(.,fooof_channels %>% filter(., Channel == 'F4' | Channel == 'F6'| Channel == 'F8') %>% group_by(Subject, Condition) %>% 
            summarise(Offset = mean(Offset, na.rm=T), 
                      Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "RDLPFC")) 
  
  #write.csv(fooof_regional, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/DLPFCFOOOFMeasures_20230111.csv')
  
  
}


CreateMRSDF <- function() {
  
  LDLPFC_GABA <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/LDLPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
  LDLPFC_Glu <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/LDLPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
  LDLPFC_Ratio <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/LDLPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
  
  LDLPFC <- LDLPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
    merge(., LDLPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"),all.x = T, all.y = T) %>%
    merge(., LDLPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
    mutate(Region = "LDLPFC")
  
  
  RDLPFC_GABA <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/RDLPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
  RDLPFC_Glu <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/RDLPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
  RDLPFC_Ratio <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/RDLPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
  
  RDLPFC <- RDLPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
    merge(., RDLPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"),all.x = T, all.y = T) %>%
    merge(., RDLPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
    mutate(Region = "RDLPFC")
  
  
  MPFC_GABA <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/MPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
  MPFC_Glu <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/MPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
  MPFC_Ratio <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/Data/MPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
  
  MPFC <- MPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
    merge(., MPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"), all.x = T, all.y = T) %>%
    merge(., MPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
    mutate(Region = "MPFC")
  
  MRSregions <- rbind(LDLPFC, RDLPFC) %>% rbind(., MPFC) %>% mutate(luna = subjID)%>% mutate(visitno = visitnum)
  #write.csv(MRSregions, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj.csv')
  
  gluRes <- lm(Glu.Cr.adj~GMrat, data = MRSregions) 
  Subject <- dplyr::select(MRSregions[!is.na(MRSregions$GMrat & MRSregions$Glu.Cr.adj),],Subject,Region)
  gluResDF <- data.frame(Subject, gluRes$residuals)
  
  gabaRes <- lm(GABA.Cr.adj~GMrat, data = MRSregions) 
  Subject <- dplyr::select(MRSregions[!is.na(MRSregions$GMrat & MRSregions$GABA.Cr.adj),],Subject,Region)
  gabaResDF <- data.frame(Subject, gabaRes$residuals)
  
  ratioRes <- lm(Ratio.adj~GMrat, data = MRSregions) 
  Subject <- dplyr::select(MRSregions[!is.na(MRSregions$GMrat & MRSregions$Ratio.adj),],Subject,Region)
  ratioResDF <- data.frame(Subject, ratioRes$residuals)
  
  MRSregionsRes <- merge(MRSregions, gluResDF, by = c("Subject", "Region"), all = T) %>% merge(.,gabaResDF, by = c("Subject", "Region"), all = T) %>% merge(.,ratioResDF, by = c("Subject", "Region"), all = T)
 
   #write.csv(MRSregionsRes, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_frontalRegions_adj_resids.csv')
  
  
  #create imbalance measure
  idx <- which(!is.na(MRSregionsRes$GABA.Cr.adj) & !is.na(MRSregionsRes$Glu.Cr.adj))
  gabaglu.lm <- lm(Glu.Cr.adj ~ GABA.Cr.adj + Region, data = MRSregionsRes[idx,])
  MRSregionsRes$GluGABAimbalance <- NA
  MRSregionsRes$GluGABAimbalanceABS <- NA
  
  
  MRSregionsRes[idx,]$GluGABAimbalanceABS <- abs(gabaglu.lm$residuals)
  MRSregionsRes[idx,]$GluGABAimbalance <- (gabaglu.lm$residuals)
  
  
  MRSregionsRes$GluMinusGABA <- MRSregionsRes$gluRes.residuals - MRSregionsRes$gabaRes.residuals
  
  MRSregionsRes <- MRSregionsRes %>% filter(Region == "LDLPFC" | Region == "RDLPFC")
  
  #write.csv(MRSregionsRes, file = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_DLPFCRegions_adj_resids.csv')
}


FooofParametersAcrossAge <- function() {
  
  #fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/DLPFCFOOOFMeasures_20230111.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  
  fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F) %>% mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult')))
  
  
  ##Exponent vs Offset
  
  lunaize(ggplot(data = fooof_regional_age, 
                 aes(y = Exponent, x = Offset, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.4) + 
            geom_point(alpha=.4) + geom_smooth(aes(group = Condition, linetype = Condition), method="lm", alpha = 0.9, size = 1) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  corDF <- fooof_regional_age %>% filter(Condition == "eyesOpen")
  
  cor.test(corDF$Exponent , corDF$Offset, method = "pearson")
  
  lm.model <- lmer(Exponent ~ Offset + age + Condition + Region + (1|luna), data = fooof_regional_age)
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model) 
  
  lm.model <- lmer(Exponent ~ Offset*age + Condition + Region + (1|luna), data = fooof_regional_age )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model) 
  
  ## Exponent vs age
  
  lunaize(ggplot(data = fooof_regional_age, 
                 aes(x = age, y = Exponent)) + 
            geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
            geom_point(aes(shape=Region),alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
    facet_wrap(~Condition) + theme(text = element_text(size = 30)) + xlab("Age") + ylab("Exponent")
  
  gam.model <-  gam(Exponent ~ s(age, k = 3)  + Condition + Region, data = fooof_regional_age  , random = ~(1|luna))
  summary(gam.model)
  
  
  # condition interaction, controlling for region
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
  # region interaction, controlling for condition
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Exponent ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
  
  ## offset vs age
  
  lunaize(ggplot(data = fooof_regional_age, 
                 aes(x = age, y = Offset)) + 
            geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
            geom_point(aes(shape=Region),alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
    facet_wrap(~Condition)+ theme(text = element_text(size = 30)) + xlab("Age") + ylab("Offset")
  
  gam.model <-  gam(Offset ~ s(age, k = 3)  + Condition + Region, data = fooof_regional_age  , random = ~(1|luna))
  summary(gam.model)
  
  
  # region interaction, controlling for condition
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Offset ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
  #condition interaction, controlling for region
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Offset ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
}


MRSmeasuresAcrossAge <- function() {
  
  #MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_DLPFCRegions_adj_resids.csv')

  ## GLU VS AGE 
  
  lunaize(ggplot(data = MRSregionsRes, 
                 aes(x = age, y = gluRes.residuals, by = luna, shape = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glutamate")+ theme(text = element_text(size = 30))
  
  gam.model <-  gam(Glu.Cr.adj ~ s(age, k = 3) + Region + GMrat, data = MRSregionsRes , random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # LDLFPC will be the reference group
  model_formula <- as.formula("Glu.Cr.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes )
  summary(model$gam)
  
  
  ## GABA VS AGE
  lunaize(ggplot(data = MRSregionsRes,
                 aes(x = age, y = gabaRes.residuals, by =luna, shape = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
            xlab("Age") +ylab("GABA") + theme(text = element_text(size = 30)))
  
  gam.model <-  gam(GABA.Cr.adj ~ s(age, k = 3) + Region + GMrat, data = MRSregionsRes , random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
  model_formula <- as.formula("GABA.Cr.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes )
  summary(model$gam)
  
  
  ## RATIO VS AGE
  lunaize(ggplot(data = MRSregionsRes, 
                 aes(x = age, y = ratioRes.residuals, by =luna, shape = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1)) + 
    xlab("Age") +ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))
  
  
  gam.model <-  gam(Ratio.adj ~ s(age, k = 3) + Region + GMrat, data = MRSregionsRes , random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
  model_formula <- as.formula("Ratio.adj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)+GMrat")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes )
  summary(model$gam)
  
  
  
  ## GLU GABA Imbalance VS AGE 
  lunaize(ggplot(data = MRSregionsRes, 
                 aes(x = age, y = GluGABAimbalanceABS, by = luna, shape = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
            xlab("Age") +ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30)))
  
  gam.model <-  gam(GluGABAimbalanceABS ~ s(age, k = 3) + GMrat + Region , data = MRSregionsRes , random = ~(1|luna))
  summary(gam.model)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
  model_formula <- as.formula("GluGABAimbalanceABS ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + GMrat")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes )
  summary(model$gam)
  
 
  
}


FOOOFvsMRS <- function() {
  
  #fooof_regional <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/DLPFCFOOOFMeasures_20230111.csv')
  #agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
  #agefile$inverseAge <- 1/agefile$age
  #MRSregionsRes <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/MRS_DLPFCRegions_adj_resids.csv')
  
  #fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject')
  
  #fooof_regional_age <- fooof_regional_age %>% separate(Subject,c("luna","vdate"), remove=F)
  
  ## MERGE MRS AND FOOOF
  fooofMRS <- merge(MRSregionsRes, fooof_regional_age, by = c("luna", "Region", "visitno", "sex")) 
  
  
  ## GLU VS Exponent
  lunaize(ggplot(data = fooofMRS , 
                 aes(x = Exponent, y = gluRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    ylab("Glutamate") + xlab("Exponent") + theme(text = element_text(size = 30))
  
  
  lm.model <- lmer(gluRes.residuals ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(gluRes.residuals ~ Exponent * age.x + Condition + Region + (1|luna), data = fooofMRS )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## GLU VS Offset 
  lunaize(ggplot(data = fooofMRS, 
                 aes(x = Offset, y = gluRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    ylab("Glutamte") + xlab("Offset") + theme(text = element_text(size = 30))
  
  
  lm.model <- lmer( gluRes.residuals ~ Offset + age.x + Condition + Region + (1|luna), data = fooofMRS )
  car::Anova(lm.model)
  summary(lm.model)
  
  
  lm.model <- lmer(gluRes.residuals ~  Offset * age.x +Condition +Region +   (1|luna), data = fooofMRS )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## GABA VS Exponent
  lunaize(ggplot(data = fooofMRS , 
                 aes(x = Exponent, y = gabaRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    ylab("GABA") + xlab("Exponent") + theme(text = element_text(size = 30))
  
  lm.model <- lmer(gabaRes.residuals ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS )
  summ(lm.model)
  
  
  lm.model <- lmer(gabaRes.residuals ~ Exponent *  age.x +Condition + Region+ (1|luna), data = fooofMRS )
  summ(lm.model)
  
  
  ## GABA VS Offset
  lunaize(ggplot(data = fooofMRS , 
                 aes(x = Offset, y = gabaRes.residuals, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    ylab("GABA") + xlab("Offset")+ theme(text = element_text(size = 30))
  
  lm.model <- lmer(gabaRes.residuals  ~ Offset + age.x + Condition + Region + (1|luna), data = fooofMRS )
  summ(lm.model)
  
  
  lm.model <- lmer(gabaRes.residuals  ~ Offset * age.x + Condition +Region +   (1|luna), data = fooofMRS )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  ## Ratio VS Exponent
  lunaize(ggplot(data = fooofMRS, 
                 aes(x = Exponent, y = ratioRes.residuals, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.1), method="lm", alpha = 0.8)  + 
            scale_color_manual(values=c("gold3", "blue4"))) +
            ylab("Glu/GABA Ratio") + xlab("Exponent")+ theme(text = element_text(size = 30))
  
  
  lm.model <- lmer(ratioRes.residuals  ~ Exponent + age.x + Condition + Region + (1|luna), data = fooofMRS )
  summ(lm.model)
  car::Anova(lm.model)
  
  
  lm.model <- lmer(ratioRes.residuals  ~ Exponent *age.x +Condition +Region+   (1|luna), data = fooofMRS )
  summ(lm.model)
  car::Anova(lm.model)
  
  
  ## Ratio VS Offset 
  lunaize(ggplot(data = fooofMRS, 
                 aes(x = Offset, y = ratioRes.residuals, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    ylab("Glu/GABA Ratio") + xlab("Offset")+ theme(text = element_text(size = 30))
  
  
  
  lm.model <- lmer(ratioRes.residuals  ~  Offset + age.x + Condition  + Region + (1|luna), data = fooofMRS)
  summ(lm.model)
  car::Anova(lm.model)
  
  
  lm.model <- lmer(ratioRes.residuals  ~ Offset *sex +age.x +Condition + Region + (1|luna), data = fooofMRS )
  summ(lm.model)
  car::Anova(lm.model)
  
  
  ## gaba glu imbalance VS offset
  lunaize(ggplot(data = fooofMRS, 
                 aes(x = Offset, y = GluGABAimbalanceABS, by =luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4")) + 
    ylab("Glu GABA Imbalance") + xlab("Offset")+ theme(text = element_text(size = 30))
  
  summary(lmerTest::lmer(GluGABAimbalanceABS ~ Offset + age.x + GMrat + Region + Condition + (1|subjID), data = fooofMRS ))
  
  summary(lmerTest::lmer(GluGABAimbalanceABS ~ Offset*age.x + GMrat + Region + Condition + (1|subjID), data = fooofMRS ))
  
  
  ## gaba glu imbalance VS exponent 
  lunaize(ggplot(data = fooofMRS, 
                 aes(y = GluGABAimbalanceABS, x = Exponent, by =luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4")) + 
    ylab("Glu GABA Imbalance") + xlab("Exponent")+ theme(text = element_text(size = 30))
  
  summary(lmerTest::lmer(GluGABAimbalanceABS ~ Exponent + age.x + GMrat + Region + Condition +(1|subjID), data = fooofMRS ))
  
  summary(lmerTest::lmer(GluGABAimbalanceABS ~ Exponent*age.x + GMrat + Region + Condition +(1|subjID), data = fooofMRS ))
  
  
  #effect size of imbalance-exponent on brain electrode map
  
  fooofallChannels <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allChannelsFOOOFMeasures_20230214.csv')
  fooofallChannels <- fooofallChannels %>% separate(Subject,c("luna","vdate"), remove=F)
  
  fooofallChannels_MRS <- merge(fooofallChannels, MRSregionsRes, by =c("luna", "visitno")) 
  
  exponentImbalanceCoef <- data.frame()
  for (chan in unique(fooofallChannels_MRS$Channel)) {
    ageRes <- lm(z~age.x + e + age.x + GMrat + Condition+ Region , data = fooofallChannels_MRS %>% filter(Channel == chan) %>% mutate(z = scale(GluGABAimbalanceABS)[,1],e = scale(Exponent)[,1] ,invage=1/age.x))   
    exponentImbalanceCoef <- rbind(exponentImbalanceCoef, data.frame(chan,ageRes$coefficients[2]))
  }
  
  exponentImbalanceCoef$labels <- exponentImbalanceCoef$chan
  exponentImbalanceCoef <- merge(exponentImbalanceCoef, chanLocs, by = "labels")
  
  lunaize(ggplot(exponentImbalanceCoef, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.03, low="white", mid="red", high="firebrick4") + ggtitle("Effect size of imbalance~exponent controlling for age, condition, region, and GMrat")  + theme(text = element_text(size = 30)))
  
  
  offsetImbalanceCoef <- data.frame()
  for (chan in unique(fooofallChannels_MRS$Channel)) {
    ageRes <- lm(z~age.x + e + age.x + GMrat + Condition+ Region , data = fooofallChannels_MRS %>% filter(Channel == chan) %>% mutate(z = scale(GluGABAimbalanceABS)[,1],e = scale(Offset)[,1] ,invage=1/age.x))   
    offsetImbalanceCoef <- rbind(offsetImbalanceCoef, data.frame(chan,ageRes$coefficients[2]))
  }
  
  offsetImbalanceCoef$labels <- offsetImbalanceCoef$chan
  offsetImbalanceCoef <- merge(offsetImbalanceCoef, chanLocs, by = "labels")
  
  lunaize(ggplot(offsetImbalanceCoef, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.025, low="white", mid="red", high="firebrick4") + ggtitle("Effect size of imbalance~offset controlling for age, condition, region, and GMrat") + theme(text = element_text(size = 30))) 
  
  
  ## Mediation

  # Glu gaba imbalance on Exponent
  mediationMatrix <- fooofMRS %>% dplyr::select(GluGABAimbalanceABS, Exponent, Condition, Region, luna, visitno, age.x, GMrat) 
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  GluGABAimbalanceABSRes <- lm(GluGABAimbalanceABS~GMrat, data = mediationMatrix) 
  residuals <- GluGABAimbalanceABSRes$residuals
  mediationMatrix$GluGABAimbalanceABSRes <- residuals
  
  
  #the effect of age on offset (c)
  model.0 <- lme4::lmer(Exponent ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  summ(model.0)
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  #the effect of age on Ratio  (a)
  model.M <- lme4::lmer(GluGABAimbalanceABSRes ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  summ(model.M)
  
  #the effect of MRS ratio on offset (b)
  model.Y <- lme4::lmer(Exponent ~ GluGABAimbalanceABSRes + age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "GluGABAimbalanceABSRes", boot = FALSE, sims = 1000)
  (summary(results))
  
  # Glu gaba imbalance on Offset
  mediationMatrix <- fooofMRS %>% dplyr::select(GluGABAimbalanceABS, Offset, Condition, Region, luna, visitno, age.x, GMrat) 
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  GluGABAimbalanceABSRes <- lm(GluGABAimbalanceABS~GMrat, data = mediationMatrix) 
  residuals <- GluGABAimbalanceABSRes$residuals
  mediationMatrix$GluGABAimbalanceABSRes <- residuals
  
  
  #the effect of age on offset (c)
  model.0 <- lme4::lmer(Offset ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  summ(model.0)
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  #the effect of age on Ratio  (a)
  model.M <- lme4::lmer(GluGABAimbalanceABSRes ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  summ(model.M)
  
  #the effect of MRS ratio on offset (b)
  model.Y <- lme4::lmer(Offset ~ GluGABAimbalanceABSRes + age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "GluGABAimbalanceABSRes", boot = FALSE, sims = 1000)
  (summary(results))
  
 
  # Glu gaba ratio on exponent
  mediationMatrix <- fooofMRS %>% dplyr::select(ratioRes.residuals, Exponent, Condition, Region, luna, visitno, age.x, GMrat) 
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  #the effect of age on offset (c)
  model.0 <- lme4::lmer(Exponent ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  summ(model.0)
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  #the effect of age on Ratio  (a)
  model.M <- lme4::lmer(ratioRes.residuals ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  summ(model.M)
  
  #the effect of MRS ratio on offset (b)
  model.Y <- lme4::lmer(Exponent ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "ratioRes.residuals", boot = FALSE, sims = 1000)
  (summary(results))
  
  
  # Glu gaba ratio on Offset
  mediationMatrix <- fooofMRS %>% dplyr::select(ratioRes.residuals, Offset, Condition, Region, luna, visitno, age.x, GMrat) 
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  #the effect of age on offset (c)
  model.0 <- lme4::lmer(Offset ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  summ(model.0)
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  #the effect of age on Ratio  (a)
  model.M <- lme4::lmer(ratioRes.residuals ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  summ(model.M)
  
  #the effect of MRS ratio on offset (b)
  model.Y <- lme4::lmer(Offset ~ ratioRes.residuals + age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "ratioRes.residuals", boot = FALSE, sims = 1000)
  (summary(results))
}


fooofvsMRSvsBehavior <- function() {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/03.BehaviorAcrossAge.R", envir = knitr::knit_global(), chdir = TRUE)
  Behavior <- Behavior_Sublevel_Maria() %>% separate(Subject,c("luna","vdate"), remove=F) %>% dplyr::select(Subject, luna, visitno, absBestError, absBestError_sd, mgsLatency, mgsLatency_sd)

  fooofMRSbehavior <- merge(fooofMRS, Behavior, by = c("luna", "visitno"))
  fooofMRSbehavior$inverseAge <- 1/fooofMRSbehavior$age.x
  
  #write.csv(fooofMRSbehavior, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/fooofMRSbehavior_20230313.csv')
  
  ## accuracy VS Exponent
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(x = absBestError, y = Exponent, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  absBestError+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior)
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Exponent ~  absBestError* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## accuracy var VS Exponent
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = Exponent, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  absBestError_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Exponent~ absBestError_sd * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## accuracy var VS offset
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = Offset, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy Var") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Offset ~  absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Offset~  absBestError_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)

  ## accuracy  VS offset
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = Offset, x = absBestError, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Offset~  absBestError+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Offset~  absBestError* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)

  
  ## latency  VS offset
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(y = Offset, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Latency") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Offset~  mgsLatency+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior)
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Offset ~ mgsLatency * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior)
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## latency var  VS offset
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = Offset, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group =1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Latency Var") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Offset~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer(Offset ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)

  
  ## latency  VS exponent
  lunaize(ggplot(data = fooofMRSbehavior , 
                 aes(y = Exponent, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Latency") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  mgsLatency+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Exponent~  mgsLatency* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
 
  
  ## looking at significant age interactions 
  
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
                 aes(y = Exponent, x = mgsLatency, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.1)) + 
    geom_point(alpha=.2) + 
    geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.5) + 
    xlab("Latency") + ylab("Exponent") + 
    theme(text = element_text(size = 30)) + 
    theme(legend.position='none') + 
    scale_color_manual(values=c("gold3", "blue4", "red4"))
  
  lm.model <- lmer(Exponent ~  mgsLatency + Condition + Region + (1|luna), data = fooofMRSbehavior%>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30")))  %>% filter(ageGroup == "17-22"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  
  ## latency var  VS exponent
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(y = Exponent, x = mgsLatency_sd, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group =1), method="lm", alpha = 0.8)) + 
    xlab("Latency Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( Exponent~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( Exponent~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)

  
  ## accuarcy VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = GluGABAimbalanceABS, x = absBestError, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy") + ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( GluGABAimbalanceABS ~ absBestError + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( GluGABAimbalanceABS ~ absBestError* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)

  ## accuracy var  VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = GluGABAimbalanceABS, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy Var") + ylab("GluGABAimbalanceABS")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( GluGABAimbalanceABS ~ absBestError_sd+ inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( GluGABAimbalanceABS ~ absBestError_sd* inverseAge + Region +GMrat + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## Lat VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = GluGABAimbalanceABS, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Latency") + ylab("Glu GABA imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( GluGABAimbalanceABS ~  mgsLatency+ inverseAge + Region + GMrat + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( GluGABAimbalanceABS ~ mgsLatency* inverseAge  + Region + GMrat + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
 
  # visualizing age interactions
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age.x, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
                 aes(y = GluGABAimbalanceABS, x = mgsLatency, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
    geom_point(alpha=.5) + 
    geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
    xlab("Latency") + ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ 
    theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))
  
  
  ## Lat Var VS imbalance
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = GluGABAimbalanceABS, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Latency Var") + ylab("Glu GABA imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(GluGABAimbalanceABS ~  mgsLatency_sd+ inverseAge  + Region + GMrat +(1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( GluGABAimbalanceABS ~ mgsLatency_sd* inverseAge  + Region + GMrat + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)

  
  ## accuracy VS ratio
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(y = ratioRes.residuals, x = absBestError, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(ratioRes.residuals~  absBestError+ age.x + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError* age.x  + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  ## accuracy variability  VS ratio
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = ratioRes.residuals, x = absBestError_sd, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError_sd+ age.xage.x + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  absBestError_sd* age.x + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  
  ## Latency VS ratio
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = ratioRes.residuals, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer( ratioRes.residuals~  mgsLatency+ age.x + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  mgsLatency* age.x  + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  
  # visualizing age interactions 
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
                 aes(y = ratioRes.residuals, x = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
    geom_point(alpha=.5) + 
    geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
    xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ 
    theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))
  
  
  ## Lat Var VS ratio
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(y = ratioRes.residuals, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) +
    xlab("Latency Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(ratioRes.residuals ~  mgsLatency_sd+ age.x  + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  
  
  lm.model <- lmer( ratioRes.residuals~  mgsLatency_sd* age.x  + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
 

  
}

