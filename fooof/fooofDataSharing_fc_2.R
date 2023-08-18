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
#library(eegUtils)
library(lme4)
library("ggpubr")
library(jtools)
library(mgcViz)
library(lubridate)
library(checkmate)

## 1.1 Define helper functions ----

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}

nsd <- function(x) {
  (abs(x - mean(x, na.rm=T)) / sd(x, na.rm=T))
}

res_with_age <- function(MRSI_input, this_roi, met_name) {
  # have columns we need
  expect_subset(c("roi", "age", "dateNumeric", met_name), names(MRSI_input))
  
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

# 2.0 Prep EEG ----
  
# agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
# chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')

agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20230112.csv')
chanLocs <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/AuditorySS/ChannelLocs.csv')


chanLocs$Channel <- chanLocs$labels

#fooof <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Resul/allSubjectsFooofMeasures_20230111.csv')

fooof <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsFooofMeasures_20230111.csv')


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
           inverseAge = 1/age) %>%
    group_by(Condition, Region) %>%
    mutate(zExp = abs(scale(Exponent)[,1]),
           zOffset = abs(scale(Offset)[,1])) %>%
    ungroup()

  str(fooof_regional_age)

  
# 3.0 Prep MRSI ----
  
## 3.1 Load MRSI data ----
  
  #MRS <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/MRSI_slice-PFC_longitudinal.csv") 
  
  MRS <- read.csv("H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/MRSI_slice-PFC_longitudinal.csv") 
  
  
  # remove garbage data based on visual inspection of LCModel fits/spectra
  # load in lcm spreadsheet that contains IDs with bad data that needs to be removed 
  #lcm <- read_excel("/Volumes/Hera/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/2_lcm_OR_120622.xlsx", col_names = FALSE)
  
  lcm <- read_excel("H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/2_lcm_OR_120622.xlsx", col_names = FALSE) 
  
  lcm <- separate(lcm, "...1", c("ld8", "junk","y","x"),extra="merge", sep = "[-.]")
  lcm <- dplyr::select(lcm, -junk)
  lcm$bad <- TRUE
  MRS <- MRS %>% mutate(x=216+1-x,y=216+1-y)
  MRS <- merge(MRS, lcm, by=c("ld8", "x", "y"), all=T) 
  MRS <- filter(MRS, is.na(bad))
  MRS <- dplyr::select(MRS, -bad)
  
  #keep people's correct coordinates
  MRS <- MRS %>% filter(!is.na(roi))
  # get rid of junk data - glu can't be zero
  MRS<- MRS %>% filter(Glu.Cr != 0)
  
  # get rid of peole who have bad data for 3 major metabolite peaks - GPC+Cho, NAA+NAAG, Cr
  MRS<- filter(MRS, GPC.Cho.SD <= 10 | is.na(GPC.Cho.SD))
  MRS <- filter(MRS, NAA.NAAG.SD <= 10 | is.na(NAA.NAAG.SD))
  MRS <- filter(MRS, Cr.SD <= 10 | is.na(Cr.SD))
  
  #  get rid of people who have lots of macromolecule in their spectra
  MRS <- filter(MRS, MM20.Cr <= 3 | is.na(MM20.Cr))
  
  #make inverse age column
  MRS$invage <- 1/MRS$age
  #make age^2 column
  MRS$age2 <- (MRS$age - mean(MRS$age))^2
  
  z_thres = 2 # standard deviation outlier detection threshold
  
  # set up some stuff for the adjustment scripts
  MRS$vdate <-ymd(gsub(".*_","", MRS$ld8))
  MRS$dateNumeric <- as.numeric(as.POSIXct(MRS$vdate, format="%m-%d-%Y"))
  MRS$relDate <- pmax(0, MRS$dateNumeric - 1560000000)
  
  # make good glu dataframe
  MRS_glu <- MRS %>%
    filter(Glu.SD <= 20) %>% #crlb 20 
    group_by(roi) %>% 
    mutate(zscore=scale(Glu.Cr, center=T, scale=T)) %>%
    filter(abs(zscore) <= z_thres) #standard deviation 
  
  # good gaba dataframe
  MRS_gaba <- MRS %>%
    filter(GABA.SD <= 20) %>%
    group_by(roi) %>%
    mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
    filter(abs(zscore) <= z_thres)
  
  # option 1 for ratio dataframe - make ratio then adjust
  MRS_ratio_1 <- MRS_glu %>%
    filter(GABA.SD <= 20) %>% # this will use only people that have both good glu and good gaba data
    group_by(roi) %>%
    mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
    filter(abs(zscore) <= z_thres) 
  
  MRS_ratio_1$glugabarat <- MRS_ratio_1$Glu.Cr/MRS_ratio_1$GABA.Cr 
  
  # this is all unadjusted data
  MRS_glu_longitudinal <- separate(MRS_glu, col=ld8, into=c('subjID', 'date'), sep='_', remove = FALSE)
  MRS_glu_longitudinal <- separate(MRS_glu_longitudinal, col=vdate, into=c('year', 'month', 'day'), sep='-', remove = FALSE)
  
  MRS_gaba_longitudinal <- separate(MRS_gaba, col=ld8, into=c('subjID', 'date'), sep='_', remove = FALSE)
  MRS_gaba_longitudinal <- separate(MRS_gaba_longitudinal, col=vdate, into=c('year', 'month', 'day'), sep='-', remove = FALSE)
  
  MRS_ratio_longitudinal <- separate(MRS_ratio_1, col=ld8, into=c('subjID', 'date'), sep='_', remove = FALSE)
  MRS_ratio_longitudinal <- separate(MRS_ratio_longitudinal, col=vdate, into=c('year', 'month', 'day'), sep='-', remove = FALSE)
  
  # adjust for date & %gm while controlling for age
  glu_adj <- res_with_age(MRS_glu_longitudinal, this_roi = c(9,10), met_name = 'Glu.Cr') %>% ungroup() %>% mutate(met_adj_new = met_adj)
  gaba_adj <- res_with_age(MRS_gaba_longitudinal, this_roi = c(9,10), met_name = 'GABA.Cr') %>% ungroup() %>% mutate(met_adj_new = met_adj)
  ratio_adj <- res_with_age(MRS_ratio_longitudinal, this_roi = c(9,10), met_name = 'glugabarat') %>% ungroup() %>% mutate(met_adj_new = met_adj)
  
  # merge together
  MRSregionsRes <- merge(
    merge(glu_adj %>% dplyr::select(ld8, subjID, age, sex, visitno=visitnum, Region=label, gluRes.residuals = met_adj),
          gaba_adj %>% dplyr::select(ld8, subjID, age, sex, visitno=visitnum, Region=label, gabaRes.residuals = met_adj),
          by=c('ld8','subjID','visitno','Region','age','sex')),
    ratio_adj %>% dplyr::select(ld8, subjID, age, sex, visitno=visitnum, Region=label, ratioRes.residuals = met_adj),
    by=c('ld8','subjID','visitno','Region','age','sex')) %>% 
    mutate(Region = ifelse(Region == 'L DLPFC','LDLPFC', ifelse(Region == 'R DLPFC','RDLPFC','err'))) %>% 
    mutate(subjID = as.factor(subjID),
           inverseAge = 1/age) %>%
    separate(ld8,c("luna","vdate"), remove=F) %>%
    mutate(luna = as.factor(luna))
    
  
 
## 3.2 Create MRSI derivatives ----
  
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
  
  str(MRSregionsRes)
  
## 4.0 FOOOF stats ----

## 4.1 Exponent vs Offset ----
  
  lunaize(ggplot(data = fooof_regional_age %>% filter(zOffset < 2), 
                 aes(y = Exponent, x = Offset, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
            geom_point(alpha=.4) + 
            stat_smooth(aes(group = Condition, linetype = Condition), color='black', method="gam", alpha = 0.4, size = 1) + scale_color_manual(values=c("gold3", "blue4")))+ theme(legend.position='none')+ theme(text = element_text(size = 30))
  
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
  
  lunaize(ggplot(data = fooof_regional_age %>% filter(zExp < 2), 
                 aes(x = age, y = Exponent)) + 
            geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
            geom_point(aes(shape=Region),alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
    facet_wrap(~Condition) + theme(text = element_text(size = 30)) + xlab("Age") + ylab("Exponent")
  
  gam.model <-  gam(Exponent ~ s(age, k = 3)  + Condition + Region + s(luna, bs='re'), 
                    data = fooof_regional_age %>% filter(zExp < 2))
  summary(gam.model)
  
  lm.model <-  lmerTest::lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooof_regional_age %>% filter(zExp < 2))
  summary(lm.model)
  
  AIC(lmer(Exponent ~ age + Condition + Region + (1|luna), data = fooof_regional_age %>% filter(zExp < 2)), 
      lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooof_regional_age %>% filter(zExp < 2)))
  
  
## 4.2.1 condition interaction, controlling for region ----
  
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region + s(luna, bs='re')")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                data = fooof_regional_age %>% filter(zExp < 2) )
  summary(model$gam)
  
  
  lm.model <-  lmerTest::lmer(Exponent ~ inverseAge*Condition + Region + (1|luna), data = fooof_regional_age)
  summary(lm.model)
  
  
## 4.2.2 region interaction, controlling for condition ----
  
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Exponent ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition + s(luna, bs='re')")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                data = fooof_regional_age  )
  summary(model$gam)
  
  
  
## 4.3 offset vs age ----
  
  lunaize(ggplot(data = fooof_regional_age %>% filter(zOffset < 2), 
                 aes(x = age, y = Offset)) + 
            geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
            geom_point(aes(shape=Region),alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
    facet_wrap(~Condition)+ theme(text = element_text(size = 30)) + xlab("Age") + ylab("Offset")
  
  gam.model <-  gam(Offset ~ s(age, k = 3)  + Condition + Region + s(luna, bs='re'), 
                    data = fooof_regional_age %>% filter(zOffset < 2))
  summary(gam.model)
  
  AIC(lmer(Offset ~ age + Condition + Region + (1|luna), data = fooof_regional_age %>% filter(zOffset < 2)), 
      lmer(Offset ~ inverseAge + Condition + Region + (1|luna), data = fooof_regional_age %>% filter(zOffset < 2)))
  
  
## 4.3.1 region interaction, controlling for condition ----
  
  fooof_regional_age$oRegion <- ordered(fooof_regional_age$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
  model_formula <- as.formula("Offset ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition + s(luna, bs='re')")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                data = fooof_regional_age  )
  summary(model$gam)
  
  
## 4.3.2 condition interaction, controlling for region ----
  
  fooof_regional_age$oCondition <- ordered(fooof_regional_age$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
  model_formula <- as.formula("Offset ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region + s(luna, bs='re')")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                data = fooof_regional_age  )
  summary(model$gam)
  
  
# 5.0 MRSI stats ----
  
## 5.1 GLU VS AGE ----
  
  lunaize(ggplot(data = MRSregionsRes %>% mutate(Region = as.factor(Region)), 
                 aes(x = age, y = gluRes.residuals, by = luna, shape = Region, color = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            #geom_smooth(aes(group = NA), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
            geom_smooth(aes(group = Region), method='lm',alpha=.8,size=1) + 
            scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glutamate")+ theme(text = element_text(size = 30))
  
  gam.model <-  mgcv::gamm(gluRes.residuals ~ s(age, k = 3) + Region,
                    random=list(luna=~1),
                    data = MRSregionsRes)
  summary(gam.model$gam)
  
  summary(mgcv::gamm(gluRes.residuals ~ s(age, k = 3) + Region, 
                     data = MRSregionsRes, random=list(luna=~1))$gam)
  

## 5.2 GABA VS AGE ----
  lunaize(ggplot(data = MRSregionsRes,
                 aes(x = age, y = gabaRes.residuals, by =luna, shape = Region, color = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
            xlab("Age") +ylab("GABA") + theme(text = element_text(size = 30)))
  
  gam.model <-  mgcv::gamm(gabaRes.residuals ~ s(age, k = 3) + Region,
                           random=list(luna=~1),
                           data = MRSregionsRes)
  summary(gam.model$gam)
  
  summary(mgcv::gamm(gabaRes.residuals ~ s(age, k = 3) + Region, 
                     data = MRSregionsRes, random=list(luna=~1))$gam)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # LDLFPC will be the reference group
  model_formula <- as.formula("gabaRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes )
  summary(model$gam)
  
  
  
## 5.3 RATIO VS AGE ----
  
  lunaize(ggplot(data = MRSregionsRes, 
                 aes(x = age, y = ratioRes.residuals, by =luna, shape = Region, color = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = Region, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1)) + 
    xlab("Age") +ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))
  
  lunaize(ggplot(data = MRSregionsRes, 
                 aes(x = age, y = ratioRes.residuals, by =luna, shape = Region)) + 
            geom_line(aes(group=interaction(luna, Region)), alpha = 0.2)) + geom_point(alpha=.5) + 
    geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
    xlab("Age") +ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))
    
  gam.model <-  mgcv::gamm(ratioRes.residuals ~ s(age, k = 3) + Region,
                           random=list(luna=~1),
                           data = MRSregionsRes)
  summary(gam.model$gam)
  
  summary(mgcv::gamm(ratioRes.residuals ~ s(age, k = 3) + Region, 
                     data = MRSregionsRes, random=list(luna=~1))$gam)
  
  MRSregionsRes$oRegion <- ordered(MRSregionsRes$Region, levels = c("LDLPFC","RDLPFC")) # LDLFPC will be the reference group
  model_formula <- as.formula("ratioRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
  # Note we keep fx = T for reliable p-values.
  model <- gamm(model_formula,
                random = list(luna=~1),
                data = MRSregionsRes )
  summary(model$gam)
  
  ggsave("dlpfc_ratio.pdf", dpi=300, height = 6, width = 8)
## 5.4 GLU GABA Imbalance VS AGE ----
  
  LNCDR::lunaize(ggplot(data = MRSregionsRes, 
                 aes(x = age, y = GluGABAimbalanceABS, by = luna, shape = Region)) + 
            geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.8,size=1) + 
            xlab("Age") +ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30)))
  
  gam.model <-  mgcv::gamm(GluGABAimbalanceABS ~ s(age, k = 3) + Region,
                           random=list(luna=~1),
                           data = MRSregionsRes)
  summary(gam.model$gam)
  
  summary(mgcv::gamm(GluGABAimbalanceABS ~ s(age, k = 3) + Region, 
                     data = MRSregionsRes, random=list(luna=~1))$gam)
  
  ggsave("dlpfc_imbalance.pdf", dpi=300, height = 6, width = 8)
  

  
# 6.0 FOOOF vs MRS ----
  

## 6.1 MERGE MRS AND FOOOF ----
  
  fooofMRS <- merge(MRSregionsRes, 
                    fooof_regional_age, 
                    by = c("luna", "Region", "visitno", "sex"), all=T) %>%
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
  gam.model <- gam(Exponent ~ gluRes.residuals + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
                   data = fooofMRS %>% filter(zExp < 2))
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
  gam.model <- gam(Offset ~ gluRes.residuals + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
                   data = fooofMRS %>% filter(zOffset < 2))
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
  gam.model <- gam(Exponent ~ gabaRes.residuals + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
                   data = fooofMRS %>% filter(zExp < 2))
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
  gam.model <- gam(Offset ~ gabaRes.residuals + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
                   data = fooofMRS %>% filter(zOffset < 2))
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
  gam.model <- gam(Exponent ~ ratioRes.residuals + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
                   data = fooofMRS %>% filter(zExp < 2))
  
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  
## 6.7 Ratio VS Offset ----
  
  lunaize(ggplot(data = fooofMRS %>% filter(abs(zOffset) < 2), 
                 aes(x = ratioRes.residuals, y = Offset, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.1), method="lm", formula=y~poly(x,2), alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Glu/GABA Ratio") + ylab("Offset")+ theme(text = element_text(size = 30))
  
  
  
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
  gam.model <- gamm(Offset ~ s(ratioRes.residuals, k=3, fx=T) + s(age.x, k=3) + Region + Condition,
                   random=list(luna=~1),
                   data = fooofMRS %>% filter(zOffset < 2))
  summary(gam.model$gam)
  print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)
  
  # quadratic
  poly.df <- fooofMRS %>% mutate(ratio.sq = scale(ratioRes.residuals^2)[,1]) %>%
    dplyr::select(Offset, ratioRes.residuals, ratio.sq, inverseAge.x, Region, Condition, luna, zOffset) %>%
    filter(complete.cases(.))
  
  lmer.base <- lmer(Offset ~ inverseAge.x + Region + Condition + (1|luna), 
                    data=poly.df %>% filter(abs(zOffset) < 2))
  lmer.quad <- lmer(Offset ~ ratioRes.residuals + ratio.sq + inverseAge.x + Region + Condition + (1|luna), 
                    data=poly.df %>% filter(abs(zOffset) < 2))
  anova(lmer.base, lmer.quad)
  summary(lmer.quad)  
  
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
#  gam.model <- gam(Offset ~ GluGABAimbalanceABS + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
  #                 data = fooofMRS %>% filter(zOffset < 2))
#  summary(gam.model)
#  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  
  gamm.model <- mgcv::gamm(Offset ~ s(GluGABAimbalance, k=3, fx=T) + 
                             s(age.x, k=3, fx=T) + Region + Condition,
                           random=list(luna=~1),
                           data = fooofMRS %>% filter(zOffset < 2))
  summary(gamm.model$gam)
  
  
  
## 6.9 gaba glu imbalance VS exponent ----
  
  lunaize(ggplot(data = fooofMRS %>% mutate(zExp = abs(scale(Exponent)[,1])) %>% filter(zExp < 2), 
                 aes(x = GluGABAimbalanceABS, y = Exponent, by =luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
    scale_color_manual(values=c("gold3", "blue4")) + 
    xlab("Glu GABA Imbalance") + ylab("Exponent")+ theme(text = element_text(size = 30))
  
  # run this one
  gam.model <-  mgcv::gamm(Exponent ~ s(age.x, k = 3) + GluGABAimbalanceABS + Condition + Region,
                           random=list(luna=~1),
                           data = fooofMRS %>% filter(zExp < 2))
  summary(gam.model$gam)
  
  
  
  
  lunaize(ggplot(data = fooofMRS %>% mutate(zExp = abs(scale(Exponent)[,1])) %>% filter(zExp < 2), 
                 aes(x = GluGABAimbalanceABS, y = Exponent, by =luna)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
    xlab("Glu GABA Imbalance") + ylab("Exponent")+ theme(text = element_text(size = 20))
  
  ggsave("imbalance_exp.pdf", dpi=300, height = 6, width = 8)
  lunaize(ggplot(data = fooofMRS %>% mutate(zExp = abs(scale(Exponent)[,1])) %>% filter(zExp < 2), 
                 aes(x = GluGABAimbalanceABS, y = Exponent, by = luna, color = age.x)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), alpha = 0.1, method="gam",formula=y~s(x,k=3),alpha=.8,size=1)) + 
    xlab("Glu GABA Imbalance") + ylab("Exponent")+ theme(text = element_text(size = 20))
  
  
  # lmer approach, different forms of age, w/ and w/out interactions
  summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + inverseAge.x + Region + Condition + (1|subjID), 
                         data = fooofMRS %>% filter(zExp < 2) ))
  
  summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS*inverseAge.x + Region + Condition + (1|subjID), 
                         data = fooofMRS %>% filter(zExp < 2) ))
  

  # gam approach, with non-linear form of age
  gam.model <- gam(Exponent ~ GluGABAimbalanceABS + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
              data = fooofMRS %>% filter(zExp < 2))
  summary(gam.model)
  
  gam.model <- gam(Exponent ~ GluGABAimbalanceABS + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
                   data = fooofMRS,
                   na.action=na.exclude)
  fooofMRS$imb_exp_resid <- residuals(gam.model)
  gam.model <- gam(Exponent ~ GluGABAimbalanceABS + s(age.x, k=3) + Region + Condition + s(luna, bs='re'), 
                   data = fooofMRS %>% mutate(z_resid = scale(imb_exp_resid)[,1]) %>% filter(abs(z_resid) < 2),
                   na.action=na.exclude)
  summary(gam.model)
  print(plot(getViz(gam.model), allTerms = T), pages = 1)
  str(fooofMRS)
  
  
  
  #quadratic
  
  lunaize(ggplot(data = fooofMRS %>% mutate(GluGABAimbalance_mod = sign(GluGABAimbalance)*atan(abs(GluGABAimbalance)), zExp = abs(scale(Exponent)[,1])) %>% filter(zExp < 2), 
                 aes(x = GluGABAimbalance, y = Exponent, by =luna, color = age.x)) + 
            geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), alpha = 0.2, method="lm", formula=y~poly(x,2),alpha=.8,size=1)) + 
    xlab("Glu GABA Imbalance") + ylab("Exponent")+ theme(text = element_text(size = 20))
  
  poly.df <- fooofMRS %>% mutate(imbal.sq = scale(GluGABAimbalance^2)[,1]) %>%
    dplyr::select(Exponent, GluGABAimbalance, imbal.sq, inverseAge.x, Region, Condition, luna) %>%
    filter(complete.cases(.))
  
  lmer.base <- lmer(Exponent ~ inverseAge.x + Region + Condition + (1|luna), data=poly.df)
  lmer.quad <- lmer(Exponent ~ GluGABAimbalance + imbal.sq + inverseAge.x + Region + Condition + (1|luna), data=poly.df)
  anova(lmer.base, lmer.quad)
  summary(lmer.quad)  
  
## 6.10 fooof-mrs loop ----
  fooofVars <- c('Offset','Exponent')
  mrsiVars <- c('ratioRes.residuals','GluGABAimbalanceABS','GluGABAimbalance','gluRes.residuals','gabaRes.residuals')
  
  output <- c()  
  for (fooofVar in fooofVars) {
    for (mrsiVar in mrsiVars) {
      
      if (fooofVar == 'Exponent') {
        fooofMRS$z <- fooofMRS$zExp
      } else {
        fooofMRS$z <- fooofMRS$zOffset
      }
        
      model <- paste0(fooofVar, ' ~ ', mrsiVar, '*inverseAge.x + Region + Condition + (1|subjID)')
      model.out <- summary(lmerTest::lmer(as.formula(model), 
                                          data = fooofMRS %>% filter(abs(z) < 2)))  
      
      t <- model.out$coefficients[2,4]
      p <- model.out$coefficients[2,5]
      t_int <- model.out$coefficients[6,4]
      p_int <- model.out$coefficients[6,5]
      
      model <- paste0(fooofVar, ' ~ ', mrsiVar, ' + inverseAge.x + Region + Condition + (1|subjID)')
      model.out <- summary(lmerTest::lmer(as.formula(model), data = fooofMRS %>% filter(abs(z) < 2)))  
            
      summary(model.out)
      
      t_meonly <- model.out$coefficients[2,4]
      p_meonly <- model.out$coefficients[2,5]
      
      gamm.model <- paste0(fooofVar, ' ~ s(', mrsiVar, ', k=3, fx=T) + s(age.x, k=3, fx=T) + Region + Condition')
      gamm.model.out <- summary(mgcv::gamm(as.formula(gamm.model),
                                   random=list(luna=~1),
                                   data = fooofMRS %>% filter(abs(z) < 2))$gam)
      p_gamm <- gamm.model.out$s.table[1,4]
      f_gamm <- gamm.model.out$s.table[1,3]
      
      output <- rbind(output, data.frame(fooofVar, mrsiVar, t, p, t_int, p_int, t_meonly, p_meonly, f_gamm, p_gamm))
    }
  }  
  output  
  
  gamm.model <- mgcv::gamm(Exponent ~ s(ratioRes.residuals, k=3, fx=T) + s(GluGABAimbalance, k=3, fx=T) + 
               s(age.x, k=3, fx=T) + Region + Condition,
             random=list(luna=~1),
             data = fooofMRS %>% filter(abs(zExp) < 2))
  summary(gamm.model$gam)
  
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
  
  mediationMatrix <- fooofMRS %>% dplyr::select(GluGABAimbalanceABS, Exponent, Condition, Region, luna, visitno, age.x) 
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  #the effect of age on offset (c)
  model.0 <- lme4::lmer(Exponent ~ age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  summ(model.0)
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  #the effect of age on Ratio  (a)
  model.M <- lme4::lmer(GluGABAimbalanceABS ~ age.x + Region + (1|luna) + Condition, data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  summ(model.M)
  
  #the effect of MRS ratio on offset (b)
  model.Y <- lme4::lmer(Exponent ~ GluGABAimbalanceABS + age.x + Condition + Region + (1|luna), data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  summ(model.Y)
  
  
  results <- mediation::mediate(model.M, model.Y, treat = "age.x", mediator = "GluGABAimbalanceABS", boot = FALSE, sims = 1000)
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
    ) %>% mutate(latencyDiff = mgsLatency - vgsLatency)
  
  str(m_clean)

  ggplot(data = pivot_longer(m_clean, cols = c(absPositionError:vgsLatency, absPositionError_sd:vgsLatency_sd, latencyDiff), names_to = 'measure', values_to = 'val'),
         aes(x = age, y = val, group = measure)) + geom_point() + stat_smooth() + facet_wrap(. ~ measure, scales = 'free')
  
  
## 9.2 Compare FOOOF & Behavior ----
  
  fooofMRSbehavior <- merge(fooofMRS, 
                            m_clean %>% separate(Subject,c("luna","vdate"), remove=F), 
                            by = c("luna", "visitno"), 
                            all=T)
  
  
## 9.2.1 accuracy VS Exponent ----
  
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(x = absBestError, y = Exponent, by = luna, color = ageGroup)) + 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior)
  car::Anova(lm.model)
  summ(lm.model)
  summary(lm.model)
  
  
  lm.model <- lmer(Exponent ~  absBestError*inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  summary(lm.model)
  
## 9.2.2 accuracy var VS Exponent ----
  
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(y = Exponent, x = absBestError_sd, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) + 
    xlab("Accuracy Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')
  
  
  lm.model <- lmer(Exponent ~  absBestError_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summ(lm.model)
  summary(lm.model)
  
  lm.model <- lmer( Exponent~ absBestError_sd * inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
  car::Anova(lm.model)
  summary(lm.model)
  summ(lm.model)
  summary(lm.model)  
  
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
  summary(lm.model)
  
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
  summary(lm.model)
  
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
  behVar <- 'latencyDiff'
  
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

  ## 9.4 fooof loop ----
  
  fooofVars <- c('Offset','Exponent')
  behVars <- c('absBestError','mgsLatency', 'vgsLatency', 'absBestError_sd', 'mgsLatency_sd','latencyDiff')

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
  
  mrsiVars <- c('ratioRes.residuals','GluGABAimbalanceABS','gluRes.residuals','gabaRes.residuals')
  behVars <- c('absBestError','mgsLatency', 'vgsLatency', 'absBestError_sd', 'mgsLatency_sd','latencyDiff')
  
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
                 aes(x = ratioRes.residuals, y = mgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) +
    theme(text = element_text(size = 14))#+ theme(legend.position='none')
  
  
  
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(x = ratioRes.residuals, y = vgsLatency, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.5) + 
            geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
            scale_color_manual(values=c("gold3", "blue4"))) +
    theme(text = element_text(size = 14))#+ theme(legend.position='none')
  
  
  ## 10.0 Residualized correlations
  
  behVars <- c('absBestError','mgsLatency', 'vgsLatency', 'absBestError_sd', 'mgsLatency_sd','latencyDiff')
  mrsiVars <- c('ratioRes.residuals','GluGABAimbalanceABS','gluRes.residuals','gabaRes.residuals')
  fooofVars <- c('Offset','Exponent')
  
  for (behVar in behVars) {
      model <- paste0(behVar, ' ~ ', 's(age, k=3)')
      gam.model <- gam(as.formula(model), data = fooofMRSbehavior, na.action = na.exclude)
      fooofMRSbehavior[[paste0(behVar, '_ageRes')]] <- residuals(gam.model)
      #print(plot(getViz(gam.model), allTerms = T), pages = 1)
  }  
  for (mrsiVar in mrsiVars) {
    model <- paste0(mrsiVar, ' ~ ', 's(age, k=3) + Region')
    gam.model <- gam(as.formula(model), data = fooofMRSbehavior, na.action = na.exclude)
    fooofMRSbehavior[[paste0(mrsiVar, '_ageRes')]] <- residuals(gam.model)
    #print(plot(getViz(gam.model), allTerms = T), pages = 1)
  }  
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ', 's(age, k=3) + Region + Condition')
    gam.model <- gam(as.formula(model), data = fooofMRSbehavior, na.action = na.exclude)
    fooofMRSbehavior[[paste0(fooofVar, '_ageRes')]] <- residuals(gam.model)
  }  
  
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(x = Offset_ageRes, y = latencyDiff_ageRes, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.3) + 
            geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.2) + 
            scale_color_manual(values=c("gold3", "blue4"))) +
    theme(text = element_text(size = 14))#+ theme(legend.position='none')
  
  
  lunaize(ggplot(data = fooofMRSbehavior, 
                 aes(x = Exponent_ageRes, y = GluGABAimbalanceABS_ageRes, by = luna, color = ageGroup))+ 
            geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
            geom_point(alpha=.3) + 
            geom_smooth(aes(group = 1),color='black', method="lm", alpha = 0.2) + 
            scale_color_manual(values=c("gold3", "blue4"))) +
    theme(text = element_text(size = 14))#+ theme(legend.position='none')
  