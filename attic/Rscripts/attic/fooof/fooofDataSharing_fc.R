# 1.0 Libraries ----

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

## 1.1 Define helper functions ----

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}

nsd <- function(x) {
  (abs(x - mean(x, na.rm=T)) / sd(x, na.rm=T))
}


# 2.0 Prep EEG ----
  
  agefile <- read.csv('~/scratch/fooof/agefile_20230112.csv')
  chanLocs <- read.csv('~/scratch/fooof/ChannelLocs.csv')
  chanLocs$Channel <- chanLocs$labels
  
  fooof <- read.csv('~/scratch/fooof/allSubjectsFooofMeasures_20230111.csv')
  
## 2.1 Remove EEG outlier per channel ----

  # REPLACING THIS (doesn't work; need mutate to respect grouping):
  # fooof_new <- fooof
  # cols = names(fooof[4:5])
  # for (col in cols) {
  #   fooof_grouped <- fooof %>% group_by(Subject, Condition)
  #   indx <- outliers(fooof_grouped[[col]])
  #   fooof_new[[paste0(col, '_nsd')]] <- nsd(fooof_grouped[[col]])
  #   fooof_new[[paste0(col, '_outlier')]] <- indx
  #   fooof_new[[col]] <- Map(replace, fooof_new[[col]], indx, NA)
  #   fooof_new[[col]] <- as.numeric(fooof_new[[col]])
  # }  
  
  fooof_new <- fooof
  fooof_new <- fooof %>% group_by(Subject, Condition) %>% mutate(Exponent = ifelse(!outliers(Exponent), Exponent, NA)) %>% ungroup
  fooof_new <- fooof_new %>% group_by(Subject, Condition) %>% mutate(Offset = ifelse(!outliers(Offset), Offset, NA)) %>% ungroup
  
## 2.2 Merge demographic and channel info ----
  
  fooof_age <- merge(fooof_new, agefile, by = 'Subject')
  
  fooof_channels <- merge(fooof_age, chanLocs, by = c("Channel"))
  
  fooof_regional <- rbind(
    fooof_channels %>% filter(Channel %in% c('F3', 'F5', 'F7')) %>% group_by(Subject, Condition) %>% 
      summarise(Offset = mean(Offset, na.rm=T), Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "LDLPFC"),
    
    fooof_channels %>% filter(Channel %in% c('F4', 'F6', 'F8')) %>% group_by(Subject, Condition) %>% 
            summarise(Offset = mean(Offset, na.rm=T), Exponent = mean(Exponent, na.rm=T)) %>% mutate(Region = "RDLPFC")
  )
  
  fooof_regional_age <- merge(fooof_regional, agefile, by = 'Subject') %>% 
    separate(Subject,c("luna","vdate"), remove=F) %>% 
    mutate(ageGroup = as.factor(ifelse(age <= 18, 'Adol', 'Adult')),
           luna = as.factor(luna),
           inverseAge = 1/age)

  str(fooof_regional_age)

  
# 3.0 Prep MRSI ----
  
## 3.1 Load MRSI data ----
  
  LDLPFC_GABA <- read.csv('~/scratch/fooof/LDLPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
  LDLPFC_Glu <- read.csv('~/scratch/fooof/LDLPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
  LDLPFC_Ratio <- read.csv('~/scratch/fooof/LDLPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
  
  LDLPFC <- LDLPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
    merge(., LDLPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), 
          by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"),all.x = T, all.y = T) %>%
    merge(., LDLPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), 
          by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
    mutate(Region = "LDLPFC")
  
  
  RDLPFC_GABA <- read.csv('~/scratch/fooof/RDLPFC_longGABA_adj_df.csv') %>% mutate(GABA.Cr.adj = met_adj)
  RDLPFC_Glu <- read.csv('~/scratch/fooof/RDLPFC_longGlu_adj_df.csv') %>% mutate(Glu.Cr.adj = met_adj)
  RDLPFC_Ratio <- read.csv('~/scratch/fooof/RDLPFC_longRatio_adj_df.csv') %>% mutate(Ratio.adj = met_adj)
  
  RDLPFC <- RDLPFC_GABA %>% dplyr::select(ld8,subjID, sex, age, GABA.Cr.adj, visitnum, GMrat) %>% mutate(Subject = ld8) %>% 
    merge(., RDLPFC_Glu %>% dplyr::select(ld8,subjID, sex, age, Glu.Cr.adj,visitnum, GMrat) %>% mutate(Subject = ld8), 
          by = c("Subject", "subjID", "sex", "age", "ld8", "visitnum", "GMrat"),all.x = T, all.y = T) %>%
    merge(., RDLPFC_Ratio %>% dplyr::select(ld8,subjID, sex, age, Ratio.adj,visitnum, GMrat) %>% mutate(Subject = ld8), 
          by = c("Subject", "subjID", "sex", "age", "ld8","visitnum", "GMrat"),all.x = T, all.y = T) %>%
    mutate(Region = "RDLPFC")
  
  MRSregions <- rbind(LDLPFC, RDLPFC)
  
## 3.2 Compute GM residuals ----

  gluRes <- lm(Glu.Cr.adj~GMrat, data = MRSregions)
  gluResDF <- cbind(MRSregions %>% filter(!is.na(GMrat) & !is.na(Glu.Cr.adj)), gluRes.residuals=gluRes$residuals) %>% 
    dplyr::select(Subject, Region, gluRes.residuals)
  
  gabaRes <- lm(GABA.Cr.adj~GMrat+Region, data = MRSregions) 
  gabaResDF <- cbind(MRSregions %>% filter(!is.na(GMrat) & !is.na(GABA.Cr.adj)), gabaRes.residuals=gabaRes$residuals) %>% 
    dplyr::select(Subject, Region, gabaRes.residuals)
  
  ratioRes <- lm(Ratio.adj~GMrat+Region, data = MRSregions) 
  ratioResDF <- cbind(MRSregions %>% filter(!is.na(GMrat) & !is.na(Ratio.adj)), ratioRes.residuals=ratioRes$residuals) %>% 
    dplyr::select(Subject, Region, ratioRes.residuals)
  
  MRSregionsRes <- merge(MRSregions, gluResDF, by = c("Subject", "Region"), all = T) %>% 
    merge(.,gabaResDF, by = c("Subject", "Region"), all = T) %>% 
    merge(.,ratioResDF, by = c("Subject", "Region"), all = T) %>% 
    mutate(subjID = as.factor(subjID),
           inverseAge = 1/age) %>%
    separate(ld8,c("luna","vdate"), remove=F) %>%
    mutate(luna = as.factor(luna))
 
## 3.3 Create MRSI derivatives ----
  
  idx <- which(!is.na(MRSregionsRes$gluRes.residuals) & !is.na(MRSregionsRes$gabaRes.residuals))
  gabaglu.lm <- lm(gluRes.residuals ~ gabaRes.residuals + Region, data = MRSregionsRes[idx,])
  MRSregionsRes$GluGABAimbalance <- NA
  MRSregionsRes$GluGABAimbalanceABS <- NA
  MRSregionsRes[idx,]$GluGABAimbalanceABS <- abs(gabaglu.lm$residuals)
  MRSregionsRes[idx,]$GluGABAimbalance <- (gabaglu.lm$residuals)
  
  MRSregionsRes$GluMinusGABA <- MRSregionsRes$gluRes - MRSregionsRes$gabaRes
  
  str(MRSregionsRes)
  
  if (F) {
    ggplot(data = MRSregionsRes, aes(x=age, y=GluGABAimbalanceABS)) + geom_point() + stat_smooth()
    summary(lmerTest::lmer(GluGABAimbalanceABS ~ age + (1|subjID), data = MRSregionsRes))
  }
  
## 4.0 FOOOF stats ----

## 4.1 Exponent vs Offset ----
  
  lunaize(ggplot(data = fooof_regional_age, 
                 aes(y = Exponent, x = Offset, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
            geom_point(alpha=.4) + 
            stat_smooth(aes(group = Condition, linetype = Condition), color='black', method="gam", alpha = 0.4, size = 1) + 
            scale_color_manual(values=c("gold3", "blue4")))
  
  corDF <- fooof_regional_age %>% filter(Condition == "eyesOpen")
  cor.test(corDF$Exponent , corDF$Offset, method = "pearson")
  
  corDF <- fooof_regional_age %>% filter(Condition == "eyesClosed")
  cor.test(corDF$Exponent , corDF$Offset, method = "pearson")
  
  
  lm.model <- lmerTest::lmer(Exponent ~ Offset + age + Condition + Region + (1|luna), data = fooof_regional_age)
  summary(lm.model)
  
  lm.model.int <- lmerTest::lmer(Exponent ~ Offset*age + Condition + Region + (1|luna), data = fooof_regional_age )
  summary(lm.model.int)
  
  AIC(lm.model, lm.model.int)

  
## 4.2 Exponent vs age ----
  
  lunaize(ggplot(data = fooof_regional_age, 
                 aes(x = age, y = Exponent)) + 
            geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
            geom_point(aes(shape=Region),alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
    facet_wrap(~Condition) + theme(text = element_text(size = 30)) + xlab("Age") + ylab("Exponent")
  
  gam.model <-  gam(Exponent ~ s(age, k = 3)  + Condition + Region, data = fooof_regional_age  , random = ~(1|luna))
  summary(gam.model)
  
  lm.model <-  lmerTest::lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooof_regional_age)
  summary(lm.model)
  
  AIC(lmer(Exponent ~ age + Condition + Region + (1|luna), data = fooof_regional_age), 
      lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooof_regional_age))
  
  
## 4.2.1 condition interaction, controlling for region ----
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
  lm.model <-  lmerTest::lmer(Exponent ~ inverseAge*Condition + Region + (1|luna), data = fooof_regional_age)
  summary(lm.model)
  
  
## 4.2.2 region interaction, controlling for condition ----
  
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Exponent ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
  
## 4.3 offset vs age ----
  
  lunaize(ggplot(data = fooof_regional_age, 
                 aes(x = age, y = Offset)) + 
            geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
            geom_point(aes(shape=Region),alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
    facet_wrap(~Condition)+ theme(text = element_text(size = 30)) + xlab("Age") + ylab("Offset")
  
  gam.model <-  gam(Offset ~ s(age, k = 3)  + Condition + Region, data = fooof_regional_age  , random = ~(1|luna))
  summary(gam.model)
  
  AIC(lmer(Offset ~ age + Condition + Region + (1|luna), data = fooof_regional_age), 
      lmer(Offset ~ inverseAge + Condition + Region + (1|luna), data = fooof_regional_age))
  
  
## 4.3.1 region interaction, controlling for condition ----
  
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Offset ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
## 4.3.2 condition interaction, controlling for region ----
  
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Offset ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = fooof_regional_age  )
  summary(model$gam)
  
  
# 5.0 MRSI stats ----
  
## 5.1 GLU VS AGE ----
  
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
  
  
## 5.2 GABA VS AGE ----
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
  
  
## 5.3 RATIO VS AGE ----
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
  
  
  
## 5.4 GLU GABA Imbalance VS AGE ----
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
  
 
# 6.0 FOOOF vs MRS ----
  

## 6.1 MERGE MRS AND FOOOF ----
  
  fooofMRS <- merge(MRSregionsRes %>% rename(visitno = visitnum), 
                    fooof_regional_age, 
                    by = c("luna", "Region", "visitno", "sex")) %>%
    mutate(zExp = abs(scale(Exponent)[,1]),
           zOffset = abs(scale(Offset)[,1])) %>%
    mutate(Region = as.factor(Region))
  
  
  
## 6.2 GLU VS Exponent ----
  
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
  
  
  summary(lmerTest::lmer(Exponent ~ gluRes.residuals + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  summary(lmerTest::lmer(Exponent ~ gluRes.residuals + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  
  summary(lmerTest::lmer(Exponent ~ gluRes.residuals*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  summary(lmerTest::lmer(Exponent ~ gluRes.residuals*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Exponent ~ gluRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zExp < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)

  
  
    
## 6.3 GLU VS Offset ----
  lunaize(ggplot(data = fooofMRS %>% filter(zOffset < 2), 
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
  
  
  summary(lmerTest::lmer(Offset ~ gluRes.residuals + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
  summary(lmerTest::lmer(Offset ~ gluRes.residuals + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
  
  summary(lmerTest::lmer(Offset ~ gluRes.residuals*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
  summary(lmerTest::lmer(Offset ~ gluRes.residuals*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
  
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Offset ~ gluRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zOffset < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Offset ~ gluRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zOffset < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  
    
## 6.4 GABA VS Exponent ----
  
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
  
  
  summary(lmerTest::lmer(Exponent ~ gabaRes.residuals + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  summary(lmerTest::lmer(Exponent ~ gabaRes.residuals + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  
  summary(lmerTest::lmer(Exponent ~ gabaRes.residuals*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  summary(lmerTest::lmer(Exponent ~ gabaRes.residuals*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Exponent ~ gabaRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zExp < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  
## 6.5 GABA VS Offset ----
  
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
  
  
  summary(lmerTest::lmer(Offset ~ gabaRes.residuals + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  summary(lmerTest::lmer(Offset ~ gabaRes.residuals + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  
  summary(lmerTest::lmer(Offset ~ gabaRes.residuals*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  summary(lmerTest::lmer(Offset ~ gabaRes.residuals*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Offset ~ gabaRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zOffset < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  

## 6.6 Ratio VS Exponent ----
  
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
  
  
  summary(lmerTest::lmer(Exponent ~ ratioRes.residuals + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  summary(lmerTest::lmer(Exponent ~ ratioRes.residuals + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  
  summary(lmerTest::lmer(Exponent ~ ratioRes.residuals*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  summary(lmerTest::lmer(Exponent ~ ratioRes.residuals*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Exponent ~ ratioRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zExp < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  
## 6.7 Ratio VS Offset ----
  
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
  
  
  summary(lmerTest::lmer(Offset ~ ratioRes.residuals + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  summary(lmerTest::lmer(Offset ~ ratioRes.residuals + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  
  summary(lmerTest::lmer(Offset ~ ratioRes.residuals*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  summary(lmerTest::lmer(Offset ~ ratioRes.residuals*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Offset ~ ratioRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zOffset < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Offset ~ ratioRes.residuals + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zOffset < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
## 6.8 gaba glu imbalance VS offset ----
  
  lunaize(ggplot(data = fooofMRS, 
                 aes(x = Offset, y = GluGABAimbalanceABS, by =luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4")) + 
    ylab("Glu GABA Imbalance") + xlab("Offset")+ theme(text = element_text(size = 30))
  

  summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  
  summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
  summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

  # gam approach, with non-linear form of age
  gam.model <- gam(Offset ~ GluGABAimbalanceABS + s(age.x, k=3) + Region + Condition, 
                   data = fooofMRS %>% filter(zOffset < 2),
                   random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  
## 6.9 gaba glu imbalance VS exponent ----
  
  lunaize(ggplot(data = fooofMRS %>% mutate(zExp = abs(scale(Exponent)[,1])) %>% filter(zExp < 2), 
                 aes(x = GluGABAimbalanceABS, y = Exponent, by =luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4")) + 
    xlab("Glu GABA Imbalance") + ylab("Exponent")+ theme(text = element_text(size = 30))
  
  
  # lmer approach, different forms of age, w/ and w/out interactions
  summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2) ))
  summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2) ))
  
  summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS*age.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2) ))
  summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS*inverseAge.x + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2) ))
  
  # gam approach, with non-linear form of age
  gam.model <- gam(Exponent ~ GluGABAimbalanceABS + s(age.x, k=3) + Region + Condition, 
              data = fooofMRS %>% filter(zExp < 2),
              random = ~(1|luna))
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  
  
# 7.0 effect size of imbalance-exponent on brain electrode map ----
  
  fooofallChannels <- read.csv('~/scratch/fooof/allChannelsFOOOFMeasures_20230214.csv')
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
  
  
# 8.0 Mediation ----

## 8.1 Glu gaba imbalance on Exponent ----
  
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
  
  
## 8.2 Glu gaba imbalance on Offset ----
  
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
  
 
## 8.3 Glu gaba ratio on exponent ----
  
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
  
  
## 8.4 Glu gaba ratio on Offset ----
  
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

  


# 9.0 FOOOF vs Beh ----
  
## 9.1 Load and prep behavior ----
  
  z_thresh = 2
  data <- merge(
    read.csv('~/scratch/fooof/Behavior_20220819.csv') %>% mutate(Subject = paste0(LunaID, '_', ScanDate)), 
    read.csv('~/scratch/fooof/agefile_20220914.csv'), 
    by='Subject')

  Behavior <- data %>%
    mutate(absPositionError = abs(PositionError)) %>%
    mutate(absBestError = abs(BestError)) %>%
    mutate(
      absPositionError = ifelse(absPositionError > 23, NA, absPositionError),
      absBestError = ifelse(absBestError > 23, NA, absBestError))
  
  
  behaviorClean <- Behavior %>%
    mutate(
      mgsLatency = ifelse(mgsLatency < .1, NA, mgsLatency), 
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      vgsLatency = ifelse(vgsLatency < .1, NA, vgsLatency), 
      vgsLatency = ifelse(abs(scale(vgsLatency)[,1])<z_thresh, vgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA)
    )
  
  
  m_mean <- behaviorClean %>%
    group_by(Subject, age, visitno) %>%
    summarize_at(vars(absPositionError, absBestError, mgsLatency, vgsLatency), mean, na.rm=T) %>%  mutate(inverseAge = 1/age)
  
  m_sd <- behaviorClean %>%
    group_by(Subject, age, visitno) %>%
    summarize_at(vars(absPositionError, absBestError, mgsLatency, vgsLatency), sd, na.rm=T) %>%  mutate(inverseAge = 1/age)
  
  m <- merge(m_mean, m_sd, by = "Subject", suffixes = c("","_sd"))
  
  # Check for subject-level outliers
  # MP - added best error
  # SM - added VGS latency
  m_clean <- m %>% ungroup() %>%
    dplyr::mutate(
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA),
      vgsLatency = ifelse(abs(scale(vgsLatency)[,1])<z_thresh, vgsLatency, NA),
      mgsLatency_sd = ifelse(abs(scale(mgsLatency_sd)[,1])<z_thresh, mgsLatency_sd, NA),
      vgsLatency_sd = ifelse(abs(scale(vgsLatency_sd)[,1])<z_thresh, vgsLatency_sd, NA),
      absPositionError_sd = ifelse(abs(scale(absPositionError_sd)[,1])<z_thresh, absPositionError_sd, NA), 
      absBestError_sd = ifelse(abs(scale(absBestError_sd)[,1])<z_thresh, absBestError_sd, NA)
    )
  
  str(m_clean)

  ggplot(data = pivot_longer(m_clean, cols = c(absPositionError:vgsLatency, absPositionError_sd:vgsLatency_sd), names_to = 'measure', values_to = 'val'),
         aes(x = age, y = val, group = measure)) + geom_point() + stat_smooth() + facet_wrap(. ~ measure, scales = 'free')
  
  
## 9.2 Compare FOOOF & Behavior ----
  
  fooofMRSbehavior <- merge(fooofMRS, m_clean %>% separate(Subject,c("luna","vdate"), remove=F), by = c("luna", "visitno"))
  
  
## 9.2.1 accuracy VS Exponent ----
  
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
  
  
## 9.2.2 accuracy var VS Exponent ----
  
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
  
  
## 9.2.3 accuracy var VS offset ----
  
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

## 9.2.4 accuracy  VS offset ----
  
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

  
## 9.2.5 latency  VS offset ----
  
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
  
  
## 9.2.6 latency var  VS offset ----
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

  
## 9.2.7 latency  VS exponent ----
  
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
 
  
## 9.2.8 looking at significant age interactions ----
  
  lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
                 aes(y = Offset, x = mgsLatency, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
    geom_point(alpha=.5) + 
    geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
    xlab("Latency") + ylab("Offset") + 
    theme(text = element_text(size = 30)) + 
    theme(legend.position='none') + 
    scale_color_manual(values=c("gold3", "blue4", "red4"))
  
  lm.model <- lmer(Exponent ~  mgsLatency + Condition + Region + (1|luna), data = fooofMRSbehavior%>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30")))  %>% filter(ageGroup == "23-30"))
  car::Anova(lm.model)
  summ(lm.model)
  
  
  
## 9.2.9 latency var  VS exponent ----
  
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

  
## 9.2.10 accuarcy VS imbalance ----
  
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

## 9.2.11 accuracy var  VS imbalance ----
  
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
  
  
## 9.2.12 Lat VS imbalance ----
  
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
  
  
## 9.2. 13 Lat Var VS imbalance ----
  
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

  
## 9.2.14 accuracy VS ratio ----
  
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
  
## 9.2.15 accuracy variability  VS ratio ----
  
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
  
  
## 9.2.16 Latency VS ratio ----
  
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
  
  
## 9.2.17 Lat Var VS ratio ----
  
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
 
## 9.3 generic ----
  
  names(fooofMRSbehavior)
  fooofVar <- 'Offset'
  behVar <- 'mgsLatency'
  
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes_string(x = fooofVar, y = behVar, by = 'luna', color = 'ageGroup'))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) +
    theme(text = element_text(size = 16))#+ theme(legend.position='none')

  model <- paste0(fooofVar, ' ~ ', behVar, '*inverseAge + Region + Condition + (1|luna)')
  summary(lmerTest::lmer(model, data = fooofMRSbehavior))  
  print(model)

  ## 9.4 loop ----
  fooofVars <- c('Offset','Exponent')
  behVars <- c('absBestError','mgsLatency', 'vgsLatency', 'absBestError_sd', 'mgsLatency_sd')

  output <- c()  
  for (fooofVar in fooofVars) {
    for (behVar in behVars) {
      model <- paste0(fooofVar, ' ~ ', behVar, '*inverseAge + Region + Condition + (1|luna)')
      model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))  
      t <- model.out$coefficients[2,4]
      p <- model.out$coefficients[2,5]
      t_int <- model.out$coefficients[6,4]
      p_int <- model.out$coefficients[6,5]
      output <- rbind(output, data.frame(fooofVar, behVar, t, p, t_int, p_int))
    }
  }  
  output  
  
  
  ## 9.5 mrsi loop ----
  
  names(fooofMRSbehavior)
  
  mrsiVars <- c('ratioRes.residuals','GluGABAimbalance','gluRes.residuals','gabaRes.residuals')
  behVars <- c('absBestError','mgsLatency', 'vgsLatency', 'absBestError_sd', 'mgsLatency_sd')
  
  output <- c()  
  for (mrsiVar in mrsiVars) {
    for (behVar in behVars) {
      model <- paste0(mrsiVar, ' ~ ', behVar, '*inverseAge + Region + Condition + (1|luna)')
      model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))  
      t <- model.out$coefficients[2,4]
      p <- model.out$coefficients[2,5]
      t_int <- model.out$coefficients[6,4]
      p_int <- model.out$coefficients[6,5]
      output <- rbind(output, data.frame(mrsiVar, behVar, t, p, t_int, p_int))
    }
  }  
  output  
  
  ## 9.6 imbalance-accuracy ----
  
  lunaize(ggplot(data = fooofMRSbehavior  , 
                 aes(x = ratioRes.residuals, y = vgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) +
    theme(text = element_text(size = 14))#+ theme(legend.position='none')
  