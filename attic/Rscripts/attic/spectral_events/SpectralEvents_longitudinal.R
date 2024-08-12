
# Libraries ----
library(LNCDR)
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
library(lme4)
library(mgcViz)
library(tidymv)
library(jtools)

# Load in merge 7T DF ----

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')

## Prep Spectral Event DF ----

SE <- merge7t[c("lunaid","visitno","eeg.age", "eeg.Delay_LDLPFC_GammaTrialPower_mean", "eeg.Delay_LDLPFC_logGammaPower_mean", 
                "eeg.Delay_RDLPFC_GammaEventDuration_mean", "eeg.Delay_RDLPFC_GammaEventNumber_mean","eeg.Delay_LDLPFC_GammaEventNumber_mean", 
                "eeg.Delay_RDLPFC_GammaTrialPower_mean", "eeg.Delay_LDLPFC_GammaEventDuration_mean",
                "eeg.Delay_RDLPFC_logGammaPower_mean", "eeg.Delay_LDLPFC_GammaEventDuration_sd", "eeg.Delay_LDLPFC_GammaEventNumber_sd",
                "eeg.Delay_LDLPFC_GammaTrialPower_sd", "eeg.Delay_LDLPFC_logGammaPower_sd",        
                "eeg.Delay_RDLPFC_GammaEventDuration_sd", "eeg.Delay_RDLPFC_GammaEventNumber_sd",    
                "eeg.Delay_RDLPFC_GammaTrialPower_sd", "eeg.Delay_RDLPFC_logGammaPower_sd")]

SElong <- SE %>% 
  select(matches('lunaid|visitno|eeg.age|(eeg).*[LR]DLPFC.*(Gamma)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure")) %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.age','Region', 'data'),
              names_from=c('measure'))

colnames(SElong)[5] <- "Epoch"
colnames(SElong)[3] <- "age"
colnames(SElong)[1] <- "luna"


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


## Prep Behavior ----

behav <- merge7t[c("lunaid","visitno","eeg.age", "eeg.BestError_DelayAll", "eeg.BestError_sd_DelayAll","eeg.mgsLatency_DelayAll", "eeg.mgsLatency_sd_DelayAll")]

colnames(behav) <- c("luna","visitno","age","absBestError","absBestError_sd","mgsLatency","mgsLatency_sd")


# Merge Data Frames ----

SElong$Region <- as.factor(SElong$Region)
MRSlong$Region <- as.factor(MRSlong$Region)

SElong$luna <- as.factor(SElong$luna)
MRSlong$luna <- as.factor(MRSlong$luna)

SE_MRS <- merge(SElong, MRSlong, by = c("luna", "visitno", "age", "Region"), all.x = T, all.y = T)

behav$luna <- as.factor(behav$luna)
SE_MRS_behavior <- merge(SE_MRS, behav, by = c("luna", "visitno", "age"), all.x = T, all.y = T)

SE_MRS_behavior$inverseAge <- 1/SE_MRS_behavior$age

# Spectral Event vs Age ----

## Trial Power ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = age, y = logGammaPower_mean)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Power (log)") + theme(legend.position = "none")


gam.model <-  gamm(logGammaPower_mean  ~ s(age, k=3) + Region, data = SE_MRS_behavior,  random=list(luna=~1))
summary(gam.model$gam)

# interaction with region
SE_MRS_behavior$oRegion <- ordered(SE_MRS_behavior$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("logGammaPower_mean ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- mgcv::gamm(model_formula,
                    random = list(luna=~1),
                    data = SE_MRS_behavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)


## Trial Power Variability----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = age, y = logGammaPower_sd)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Power Var(log)") + theme(legend.position = "none")


gam.model <-  gamm(logGammaPower_sd  ~ s(age, k=3) + Region, data = SE_MRS_behavior,  random=list(luna=~1))
summary(gam.model$gam)

# interaction with region
SE_MRS_behavior$oRegion <- ordered(SE_MRS_behavior$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("logGammaPower_sd ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- mgcv::gamm(model_formula,
                    random = list(luna=~1),
                    data = SE_MRS_behavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)



## Event Number ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = age, y = GammaEventNumber_mean)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Event Number") + theme(legend.position = "none")


gam.model <-  gamm(GammaEventNumber_mean  ~ s(age, k=3) + Region, data = SE_MRS_behavior,  random=list(luna=~1))
summary(gam.model$gam)

# interaction with region
SE_MRS_behavior$oRegion <- ordered(SE_MRS_behavior$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("GammaEventNumber_mean ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- mgcv::gamm(model_formula,
                    random = list(luna=~1),
                    data = SE_MRS_behavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)

## Event Number Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = age, y = GammaEventNumber_sd)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Event Number Var") + theme(legend.position = "none")


gam.model <-  gamm(GammaEventNumber_sd  ~ s(age, k=3) + Region, data = SE_MRS_behavior,  random=list(luna=~1))
summary(gam.model$gam)

# interaction with region
SE_MRS_behavior$oRegion <- ordered(SE_MRS_behavior$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("GammaEventNumber_sd ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- mgcv::gamm(model_formula,
                    random = list(luna=~1),
                    data = SE_MRS_behavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)

## Event Duration ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = age, y = GammaEventDuration_mean)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Event Duration") + theme(legend.position = "none")


gam.model <-  gamm(GammaEventDuration_mean  ~ s(age, k=3) + Region, data = SE_MRS_behavior,  random=list(luna=~1))
summary(gam.model$gam)

# interaction with region
SE_MRS_behavior$oRegion <- ordered(SE_MRS_behavior$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("GammaEventDuration_mean ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- mgcv::gamm(model_formula,
                    random = list(luna=~1),
                    data = SE_MRS_behavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)

## Event Duration Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = age, y = GammaEventDuration_sd)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Event Duration Var") + theme(legend.position = "none")


gam.model <-  gamm(GammaEventDuration_sd  ~ s(age, k=3) + Region, data = SE_MRS_behavior,  random=list(luna=~1))
summary(gam.model$gam)

# interaction with region
SE_MRS_behavior$oRegion <- ordered(SE_MRS_behavior$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("GammaEventDuration_sd ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- mgcv::gamm(model_formula,
                    random = list(luna=~1),
                    data = SE_MRS_behavior %>% filter(Region == "LDLPFC" | Region == "RDLPFC"))
summary(model$gam)


# MRS vs Spectral Events ----

#SE_MRS_behavior <- SE_MRS_behavior[complete.cases(SE_MRS_behavior), ]
## Imbalance ----

### Trial Power ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_mean, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior), 
    lmer(GluGABAimbalanceABS ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior))


### Trial Power Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_sd, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power Var") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ logGammaPower_sd + inverseAge + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Number ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_mean, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventNumber_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Number Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_sd, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventNumber_sd + inverseAge + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Duration ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_mean, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventDuration_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Duration Var ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_sd, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventDuration_sd + inverseAge + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

## Imbalance (Adults Only) ----

### Trial Power ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = logGammaPower_mean, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior), 
    lmer(GluGABAimbalanceABS ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior))


### Trial Power Variability ----

lunaize(ggplot(data = SE_MRS_behavior%>% filter(age>18), 
               aes(x = logGammaPower_sd, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power Var") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ logGammaPower_sd + inverseAge + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)

### Event Number ----

lunaize(ggplot(data = SE_MRS_behavior%>% filter(age>18), 
               aes(x = GammaEventNumber_mean, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventNumber_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)


### Event Number Variability ----

lunaize(ggplot(data = SE_MRS_behavior%>% filter(age>18), 
               aes(x = GammaEventNumber_sd, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventNumber_sd + inverseAge + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)

### Event Duration ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = GammaEventDuration_mean, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventDuration_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)


### Event Duration Var ----

lunaize(ggplot(data = SE_MRS_behavior%>% filter(age>18), 
               aes(x = GammaEventDuration_sd, y = GluGABAimbalanceABS)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var") + ylab("Glu GABA Imbalance") + theme(legend.position = "none")

lm.model <- lmer(GluGABAimbalanceABS ~ GammaEventDuration_sd + inverseAge + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)




## Ratio ----

### Trial Power ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_mean, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


AIC(lmer(Ratio_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior), 
    lmer(Ratio_gamadj ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior))


### Trial Power Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_sd, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power Var") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ logGammaPower_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Number ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_mean, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventNumber_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Number Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_sd, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventNumber_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Duration ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_mean, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventDuration_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Duration Var ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_sd, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventDuration_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


## Ratio (Adults Only) ----

### Trial Power ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = logGammaPower_mean, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior %>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)


AIC(lmer(Ratio_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior), 
    lmer(Ratio_gamadj ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior))


### Trial Power Variability ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = logGammaPower_sd, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power Var") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ logGammaPower_sd + age + Region + (1|luna), data = SE_MRS_behavior %>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)

### Event Number ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = GammaEventNumber_mean, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventNumber_mean + age + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)


### Event Number Variability ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = GammaEventNumber_sd, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventNumber_sd + age + Region + (1|luna), data = SE_MRS_behavior%>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)

### Event Duration ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = GammaEventDuration_mean, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventDuration_mean + age + Region + (1|luna), data = SE_MRS_behavior %>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)


### Event Duration Var ----

lunaize(ggplot(data = SE_MRS_behavior %>% filter(age>18), 
               aes(x = GammaEventDuration_sd, y = Ratio_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var") + ylab("Glu GABA Ratio") + theme(legend.position = "none")

lm.model <- lmer(Ratio_gamadj ~ GammaEventDuration_sd + age + Region + (1|luna), data = SE_MRS_behavior %>% filter(age>18))
car::Anova(lm.model)
summ(lm.model)



## Glu ----

### Trial Power ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_mean, y = Glu_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power") + ylab("Glu") + theme(legend.position = "none")

lm.model <- lmer(Glu_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior), 
    lmer(Glu_gamadj ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior))



### Trial Power Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_sd, y = Glu_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power Var") + ylab("Glu") + theme(legend.position = "none")

lm.model <- lmer(Glu_gamadj ~ logGammaPower_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Number ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_mean, y = Glu_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number") + ylab("Glu") + theme(legend.position = "none")

lm.model <- lmer(Glu_gamadj ~ GammaEventNumber_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Number Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_sd, y = Glu_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var") + ylab("Glu") + theme(legend.position = "none")

lm.model <- lmer(Glu_gamadj ~ GammaEventNumber_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Duration ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_mean, y = Glu_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration") + ylab("Glu") + theme(legend.position = "none")

lm.model <- lmer(Glu_gamadj ~ GammaEventDuration_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Duration Var ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_sd, y = Glu_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var") + ylab("Glu") + theme(legend.position = "none")

lm.model <- lmer(Glu_gamadj ~ GammaEventDuration_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

## GABA ----

### Trial Power ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_mean, y = GABA_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power") + ylab("GABA") + theme(legend.position = "none")

lm.model <- lmer(GABA_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ logGammaPower_mean + age + Region + (1|luna), data = SE_MRS_behavior), 
    lmer(GABA_gamadj ~ logGammaPower_mean + inverseAge + Region + (1|luna), data = SE_MRS_behavior))



### Trial Power Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = logGammaPower_sd, y = GABA_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Trial Power Var") + ylab("GABA") + theme(legend.position = "none")

lm.model <- lmer(GABA_gamadj ~ logGammaPower_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Number ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_mean, y = GABA_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number") + ylab("GABA") + theme(legend.position = "none")

lm.model <- lmer(GABA_gamadj ~ GammaEventNumber_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Number Variability ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventNumber_sd, y = GABA_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var") + ylab("GABA") + theme(legend.position = "none")

lm.model <- lmer(GABA_gamadj ~ GammaEventNumber_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)

### Event Duration ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_mean, y = GABA_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration") + ylab("GABA") + theme(legend.position = "none")

lm.model <- lmer(GABA_gamadj ~ GammaEventDuration_mean + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


### Event Duration Var ----

lunaize(ggplot(data = SE_MRS_behavior, 
               aes(x = GammaEventDuration_sd, y = GABA_gamadj)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var") + ylab("GABA") + theme(legend.position = "none")

lm.model <- lmer(GABA_gamadj ~ GammaEventDuration_sd + age + Region + (1|luna), data = SE_MRS_behavior)
car::Anova(lm.model)
summ(lm.model)


# Predicting Trajectories ----
 
## Create Wide Data Frame ---- 


MRSlong <- MRS7t %>% 
  select(matches('lunaid|visitno|eeg.age|(sipfc).*[LR]DLPFC.*(gamadj)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.age','Region'),
              names_from=c('data','measure'))



wider_data <- SE_MRS_behavior %>% 
  pivot_wider(id_cols=c('luna','Region'), names_from = visitno, values_from = c("age","GammaTrialPower_mean","logGammaPower_mean","GammaEventDuration_mean","GammaEventNumber_mean",   
                                                    "GammaEventDuration_sd", "GammaEventNumber_sd", "logGammaPower_sd", "GammaTrialPower_sd",
                                                    "GABA_gamadj", "Glu_gamadj", "GluGABAimbalance", "GluGABAimbalanceABS","Ratio_gamadj"))

wider_data$ageDiff <- wider_data$age_2 - wider_data$age_1
wider_data$ageMean <- (wider_data$age_2 + wider_data$age_1)/2



## Imbalance ----
wider_data$GluGABAimbalance_diff <- wider_data$GluGABAimbalance_2 - wider_data$GluGABAimbalance_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_sd_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_sd_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_sd_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


## Imbalance (Under 18) ----
wider_data$GluGABAimbalance_diff <- wider_data$GluGABAimbalance_2 - wider_data$GluGABAimbalance_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data %>% filter(age_1<18), 
               aes(x = logGammaPower_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = logGammaPower_sd_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data %>% filter(age_1<18), 
               aes(x = GammaEventNumber_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventNumber_sd_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventDuration_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventDuration_sd_diff, y = GluGABAimbalance_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("Imbalance Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ GluGABAimbalance_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)



## Ratio ----
wider_data$ratio_diff <- wider_data$Ratio_gamadj_2 - wider_data$Ratio_gamadj_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_sd_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_sd_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_sd_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)

## Ratio (Under 18) ----
wider_data$ratio_diff <- wider_data$Ratio_gamadj_2 - wider_data$Ratio_gamadj_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data  %>% filter(age_1<18), 
               aes(x = logGammaPower_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data  %>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data  %>% filter(age_1<18), 
               aes(x = logGammaPower_sd_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data %>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data  %>% filter(age_1<18), 
               aes(x = GammaEventNumber_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data %>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data %>% filter(age_1<18), 
               aes(x = GammaEventNumber_sd_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data %>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data  %>% filter(age_1<18), 
               aes(x = GammaEventDuration_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data %>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data %>% filter(age_1<18), 
               aes(x = GammaEventDuration_sd_diff, y = ratio_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("Ratio Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ ratio_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data %>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


## GABA ----
wider_data$GABA_diff <- wider_data$GABA_gamadj_2 - wider_data$GABA_gamadj_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_sd_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_sd_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_sd_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)

## GABA (under 18) ----
wider_data$GABA_diff <- wider_data$GABA_gamadj_2 - wider_data$GABA_gamadj_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = logGammaPower_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = logGammaPower_sd_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventNumber_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.5) + 
          geom_point(aes(shape=Region), alpha = 0.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventNumber_sd_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventDuration_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventDuration_sd_diff, y = GABA_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("GABA Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ GABA_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


## Glu ----
wider_data$Glu_diff <- wider_data$Glu_gamadj_2 - wider_data$Glu_gamadj_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = logGammaPower_sd_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventNumber_sd_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data, 
               aes(x = GammaEventDuration_sd_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)

## Glu (under 18) ----
wider_data$Glu_diff <- wider_data$Glu_gamadj_2 - wider_data$Glu_gamadj_1

### Trial Power ----
wider_data$logGammaPower_diff <- wider_data$logGammaPower_mean_2 - wider_data$logGammaPower_mean_1

diffpredictModeldiff <- lmer(logGammaPower_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data) # does the difference from v2-v1 in eeg predict the difference in mrs 

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = logGammaPower_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Trial Power Var ----
wider_data$logGammaPower_sd_diff <- wider_data$logGammaPower_sd_2 - wider_data$logGammaPower_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = logGammaPower_sd_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Gamma Power Var Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(logGammaPower_sd_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number ----
wider_data$GammaEventNumber_diff <- wider_data$GammaEventNumber_mean_2 - wider_data$GammaEventNumber_mean_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventNumber_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Number Var ----
wider_data$GammaEventNumber_sd_diff <- wider_data$GammaEventNumber_sd_2 - wider_data$GammaEventNumber_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventNumber_sd_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Number Var Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventNumber_sd_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration ----
wider_data$GammaEventDuration_diff <- wider_data$GammaEventDuration_mean_2 - wider_data$GammaEventDuration_mean_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventDuration_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)


### Event Duration Var ----
wider_data$GammaEventDuration_sd_diff <- wider_data$GammaEventDuration_sd_2 - wider_data$GammaEventDuration_sd_1

lunaize(ggplot(data = wider_data%>% filter(age_1<18), 
               aes(x = GammaEventDuration_sd_diff, y = Glu_diff)) + 
          geom_line(aes(group=interaction(luna,Region)),alpha = 0.2) + 
          geom_point(aes(shape = Region),alpha=.2) +
          geom_smooth(aes(group = Region, alpha = 0.2, linetype = Region), method="lm",alpha=0.5,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Event Duration Var Diff") + ylab("Glu Diff") + theme(legend.position = "none")

diffpredictModeldiff <- lmer(GammaEventDuration_sd_diff ~ Glu_diff + ageDiff + ageMean + Region + (1|luna), data = wider_data%>% filter(age_1<18)) # does the difference from v2-v1 in eeg predict the difference in mrs 
car::Anova(diffpredictModeldiff)

