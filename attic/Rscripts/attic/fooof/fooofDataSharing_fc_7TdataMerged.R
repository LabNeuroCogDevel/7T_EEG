# 1.0 Libraries ----

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

## 1.1 Define helper functions ----

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}

nsd <- function(x) {
  (abs(x - mean(x, na.rm=T)) / sd(x, na.rm=T))
}

res_with_age <- function(MRSI_input, this_roi, met_name) {
  # have columns we need
  checkmate::expect_subset(c("roi", "age", "dateNumeric", met_name), names(MRSI_input))
  
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
# 2.0 Load in Data Frames ----

## 2.1 Prep EEG ----

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
fooof <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "eeg.eyesOpen_LDLPFC_Offset","eeg.eyesClosed_LDLPFC_Offset","eeg.eyesOpen_RDLPFC_Offset","eeg.eyesClosed_RDLPFC_Offset","eeg.eyesClosed_LDLPFC_Exponent", "eeg.eyesOpen_LDLPFC_Exponent","eeg.eyesOpen_RDLPFC_Exponent","eeg.eyesClosed_RDLPFC_Exponent")]

fooofLong <- fooof %>% 
  select(matches('lunaid|visitno|eeg.age|(eeg).*[LR]DLPFC.*(Exponent|Offset)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','eeg.age','Region'),
              names_from=c('data','measure')) %>% 
  select(matches('lunaid|visitno|Region|eeg.age|(eeg).*(Exponent|Offset)')) %>%
  pivot_longer(cols=matches('eyes'),
               names_pattern='(.*).(eyesOpen|eyesClosed)_(.*)',
               names_to=c("data", "Condition","measure"))  %>% 
  filter(!is.na(value))  %>% 
  pivot_wider(id_cols=c('lunaid','visitno','Region','eeg.age','Condition'),
              names_from=c('measure'))

fooofLong <- fooofLong %>% 
  mutate(ageGroup = as.factor(ifelse(eeg.age <= 18, 'Adol', 'Adult')),
         lunaid = as.factor(lunaid),
         inverseAge = 1/eeg.age) 

colnames(fooofLong) <- c("luna","visitno", "Region", "age", "Condition", "Offset", "Exponent", "ageGroup", "inverseAge")


write.csv(fooofLong,"/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsDLPFCfooofMeasures_20230613.csv")


### 2.1.1 FOOOF error measures ---- 
fooofErrors <- read.csv(hera('/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/allSubjectsErrorMeasures_20230516.csv'))

fooof_channels <- merge(fooofErrors, chanLocs, by = c("Channel"))

fooof_regional <- rbind(
  fooof_channels %>% filter(Channel %in% c('F3', 'F5', 'F7')) %>% group_by(Subject, Condition) %>% 
    summarise(Error = mean(Error, na.rm=T), Rsqaured = mean(R.Squared, na.rm=T)) %>% mutate(Region = "LDLPFC"),
  
  fooof_channels %>% filter(Channel %in% c('F4', 'F6', 'F8')) %>% group_by(Subject, Condition) %>% 
    summarise(Error = mean(Error, na.rm=T), Rsqaured = mean(R.Squared, na.rm=T)) %>% mutate(Region = "RDLPFC")
)

avgError <- aggregate(.~Condition, data = fooof_regional[2:4], mean)



## 2.2  Prep MRS ----
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

### 2.2.1 Create MRSI derivatives ----

idx <- which(!is.na(MRSlong$Glu_gamadj) & !is.na(MRSlong$GABA_gamadj))
gabaglu.lm <- lm(Glu_gamadj ~ GABA_gamadj + Region, data = MRSlong[idx,])
MRSlong$GluGABAimbalance <- NA
MRSlong$GluGABAimbalanceABS <- NA
MRSlong[idx,]$GluGABAimbalanceABS <- abs(gabaglu.lm$residuals)
MRSlong[idx,]$GluGABAimbalance <- (gabaglu.lm$residuals)

MRSlong$Ratio_gamadj <- MRSlong$Glu_gamadj/MRSlong$GABA_gamadj 



## 2.3 Prep Behavior ----

behav <- merge7t[c("lunaid","visitno","eeg.age", "eeg.BestError_DelayAll", "eeg.BestError_sd_DelayAll","eeg.mgsLatency_DelayAll", "eeg.mgsLatency_sd_DelayAll","cantab.ssp_max.span","cantab.ssp_nerrors", "cantab.ssp_ntrials")]

colnames(behav) <- c("luna","visitno","age","absBestError","absBestError_sd","mgsLatency","mgsLatency_sd", "SSP_maxSpan", "SSP_nErrors", "SSP_nTrials")

### 2.3.1 Behavior vs age ----

#### Best Saccade ----
lunaize(ggplot(data = behav, 
               aes(y = absBestError, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='black', method="lm", formula = 'y~I(1/x)', alpha = 0.4, size = 1) + 
          scale_color_manual(values=c("gold3", "blue4")))

gam.model <- gam(absBestError ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, agevar = 'age', yvar = 'absBestError', idvar = "luna" ,draw_points = T, xplotname = "Age", yplotname = "Mean MGS Accuracy (degs)"))


gam.model <-  gamm(absBestError ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)


#### Best Saccade Var ----
gam.model <- gam(absBestError_sd ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, agevar = 'age', idvar = "luna", yvar = 'absBestError_sd', draw_points = T, xplotname = "Age", yplotname = "SD MGS Accuracy (degs)"))

gam.model <-  gamm(absBestError_sd ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

#### MGS Latency ----
gam.model <- gam(mgsLatency ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'mgsLatency', draw_points = T, xplotname = "Age", yplotname = "Mean MGS Latency (s)"))

gam.model <-  gamm(mgsLatency ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

#### MGS Latency Var ----
gam.model <- gam(mgsLatency_sd ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'mgsLatency_sd', draw_points = T, xplotname = "Age", yplotname = "SD MGS Latency (s)"))

gam.model <-  gamm(mgsLatency_sd ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

#### Spatial Span Max ----

lunaize(ggplot(data = behav, 
               aes(y = SSP_maxSpan, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='blue', method="lm", formula = 'y~I(1/x)', alpha = 0.4))

gam.model <-  gamm(SSP_maxSpan ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

lm.model <- lmerTest::lmer(SSP_maxSpan ~ age + (1|luna), data = behav)
summary(lm.model)


gam.model <- gam(SSP_maxSpan ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'SSP_maxSpan', draw_points = T, xplotname = "Age", yplotname = "Max Sequence Length"))

gam.model <-  gamm(SSP_maxSpan ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

AIC(lmer(Exponent ~ age + Condition + Region + (1|luna), data = fooofLong), 
    lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooofLong))


#### Spatial Span Errors ----

lunaize(ggplot(data = behav, 
               aes(y = SSP_nErrors, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='blue', method="lm", formula = 'y~I(1/x)', alpha = 0.4))

gam.model <-  gamm(SSP_nErrors ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

lm.model <- lmerTest::lmer(SSP_nErrors ~ age + (1|luna), data = behav)
summary(lm.model)


gam.model <- gam(SSP_nErrors ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'SSP_nErrors', draw_points = T, xplotname = "Age", yplotname = "SSP nErrors"))

#### Spatial Span nTrials ----

lunaize(ggplot(data = behav, 
               aes(y = SSP_nTrials, x = age, by = luna))+ 
          geom_line(aes(group=luna), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = 1), color='blue', method="lm", formula = 'y~I(1/x)', alpha = 0.4))

gam.model <-  gamm(SSP_nTrials ~ s(age, k = 3), data = behav, random=list(luna=~1))
summary(gam.model$gam)

lm.model <- lmerTest::lmer(SSP_nTrials ~ age + (1|luna), data = behav)
summary(lm.model)


gam.model <- gam(SSP_nTrials ~ s(age), data = behav, random=list(luna=~1))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
lunaize(gam_growthrate_plot(behav, gam.model, gam.growthrate, idvar = "luna", agevar = 'age', yvar = 'SSP_nTrials', draw_points = T, xplotname = "Age", yplotname = "SSP nTrials"))


# 3.0 Merge Data Frames ----

fooofLong$Region <- as.factor(fooofLong$Region)
MRSlong$Region <- as.factor(MRSlong$Region)

fooofLong$luna <- as.factor(fooofLong$luna)
MRSlong$luna <- as.factor(MRSlong$luna)

fooofMRS <- merge(fooofLong, MRSlong, by = c("luna", "visitno", "age", "Region"))

behav$luna <- as.factor(behav$luna)
fooofMRSbehavior <- merge(fooofMRS, behav, by = c("luna", "visitno", "age"))

# 4.0 FOOOF stats ----

## 4.1 Exponent vs Offset ----

lunaize(ggplot(data = fooofLong, 
               aes(y = Exponent, x = Offset, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(alpha=.4) + 
          stat_smooth(aes(group = Condition, linetype = Condition), color='black', method="gam", alpha = 0.4, size = 1) + 
          scale_color_manual(values=c("gold3", "blue4")))

corDF <- fooofLong %>% filter(Condition == "eyesOpen")
cor.test(corDF$Exponent , corDF$Offset, method = "pearson")

corDF <- fooofLong %>% filter(Condition == "eyesClosed")
cor.test(corDF$Exponent , corDF$Offset, method = "pearson")


lm.model <- lmerTest::lmer(Exponent ~ Offset + age + Condition + Region + (1|luna), data = fooofLong)
summary(lm.model)

lm.model.int <- lmerTest::lmer(Exponent ~ Offset*age + Condition + Region + (1|luna), data = fooofLong )
summary(lm.model.int)

AIC(lm.model, lm.model.int)


## 4.2 Exponent vs age ----

lunaize(ggplot(data = fooofLong %>% filter(Condition == 'eyesOpen'), 
               aes(x = age, y = Exponent)) + 
          geom_line(aes(group=interaction(luna,Region), shape =Region), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
  theme(text = element_text(size = 30)) + xlab("Age") + ylab("Exponent") + theme(legend.position = "none")

gam.model <-  gamm(Exponent ~ s(age, k = 3)  + Condition + Region, data = fooofLong, random=list(luna=~1))
summary(gam.model$gam)

lm.model <-  lmerTest::lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooofLong)
summary(lm.model)

AIC(lmer(Exponent ~ age + Condition + Region + (1|luna), data = fooofLong), 
    lmer(Exponent ~ inverseAge + Condition + Region + (1|luna), data = fooofLong))


### 4.2.1 condition interaction, controlling for region ----
fooofLong$oCondition <- ordered(fooofLong$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
model_formula <- as.formula("Exponent ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)


lm.model <-  lmerTest::lmer(Exponent ~ inverseAge*Condition + Region + (1|luna), data = fooofLong)
summary(lm.model)


### 4.2.2 region interaction, controlling for condition ----

fooofLong$oRegion <- ordered(fooofLong$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
model_formula <- as.formula("Exponent ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)



## 4.3 Offset vs age ----

lunaize(ggplot(data = fooofLong%>% filter(Condition == 'eyesClosed'), 
               aes(x = age, y = Offset)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + 
          geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = 1, alpha = 0.2), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=0.8,size=1)) + 
   theme(text = element_text(size = 30)) + xlab("Age") + ylab("Offset") + theme(legend.position = "none")

gam.model <-  gamm(Offset ~ s(age, k = 3)  + Condition + Region, data = fooofLong, random=list(luna=~1))
summary(gam.model$gam)

AIC(lmer(Offset ~ age + Condition + Region + (1|luna), data = fooofLong), 
    lmer(Offset ~ inverseAge + Condition + Region + (1|luna), data = fooofLong))


### 4.3.1 region interaction, controlling for condition ----

fooofLong$oRegion <- ordered(fooofLong$Region, levels = c("LDLPFC","RDLPFC")) # ldlpfc will be the reference
model_formula <- as.formula("Offset ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)


### 4.3.2 condition interaction, controlling for region ----

fooofLong$oCondition <- ordered(fooofLong$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group
model_formula <- as.formula("Offset ~ oCondition + s(age, k = 4, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = fooofLong  )
summary(model$gam)


# 5.0 MRSI stats ----

## 5.1 GLU VS AGE ----

lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = Glu_gamadj, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
          geom_smooth(aes(group = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1) + 
          scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glutamate")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

#ggsave("dlpfc_glu_age.png", height=6, width = 6, dpi = 300)
gam.model <-  gamm(Glu_gamadj ~ s(age, k = 3) + Region, data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

#test <-  lmer(Glu_gamadj ~ inverseAge * Region + (1|luna), data = MRSlong)
#summary(test)


MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLFPC will be the reference group
model_formula <- as.formula("Glu_gamadj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong )
summary(model$gam)


## 5.2 GABA VS AGE ----
lunaize(ggplot(data = MRSlong,
               aes(x = age, y = GABA_gamadj, by =luna, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
          geom_smooth(aes(group = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1)) + 
  xlab("Age") +ylab("GABA") + theme(text = element_text(size = 30)) + theme(legend.position='none')


ggsave("dlpfc_gaba_age.png", height=6, width = 6, dpi = 300)

gam.model <-  gamm(GABA_gamadj ~ s(age, k = 3) + Region, data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("gabaRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong )
summary(model$gam)


## 5.3 RATIO VS AGE ----
lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = Ratio_gamadj, by =luna, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = Region, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1)) + 
  xlab("Age") +ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


gam.model <-  gamm(Ratio_gamadj ~ s(age, k = 3) + Region, data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("Ratio_gamadj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong%>% filter(age < 26) )
summary(model$gam)



## 5.4 GLU GABA Imbalance VS AGE ----
lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = GluGABAimbalanceABS, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(aes(shape= Region),alpha=.5) + 
          geom_smooth(aes(group = Region, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1)) + 
  xlab("Age") +ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

gam.model <-  gamm(GluGABAimbalanceABS ~ s(age, k = 3) +Region , data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
model_formula <- as.formula("GluGABAimbalanceABS ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong )
summary(model$gam)


# 6.0 FOOOF vs MRS ----

## 6.1 GLU VS Exponent ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Exponent, y = Glu_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) +
          geom_point(alpha=.5) + geom_smooth(aes(group = Region, linetype = Region), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glutamate") + xlab("Exponent") + theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Glu_gamadj ~ Exponent + age + Condition + Region + (1|luna), data = fooofMRS )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(Glu_gamadj ~ Exponent * age + Condition + Region + (1|luna), data = fooofMRS )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)


summary(lmerTest::lmer(Exponent ~ Glu_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Glu_gamadj + inverseAge + Region + Condition + (1|luna), data = fooofMRS))

summary(lmerTest::lmer(Exponent ~ Glu_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Glu_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))


# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ Glu_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS, random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)


AIC(lmer(Exponent ~ Glu_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ Glu_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))



## 6.2 GLU VS Offset ----
lunaize(ggplot(data = fooofMRS, 
               aes(x = Offset, y = Glu_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + geom_smooth(aes(group = Region, linetype=Region), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glutamte") + xlab("Offset") + theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj ~ Offset + age + Condition + Region + (1|luna), data = fooofMRS )
car::Anova(lm.model)
summary(lm.model)


lm.model <- lmer(Glu_gamadj ~  Offset * age +Condition +Region +   (1|luna), data = fooofMRS )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)


summary(lmerTest::lmer(Offset ~ Glu_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
summary(lmerTest::lmer(Offset ~ Glu_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))

summary(lmerTest::lmer(Offset ~ Glu_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))
summary(lmerTest::lmer(Offset ~ Glu_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2)))


# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Glu_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Glu_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Offset ~ Glu_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ Glu_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))




## 6.3 GABA VS Exponent ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Exponent, y = GABA_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + geom_smooth(aes(group = Region, linetype=Region), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("GABA") + xlab("Exponent") + theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer(GABA_gamadj ~ Exponent + age + Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)


lm.model <- lmer(GABA_gamadj ~ Exponent *  age +Condition + Region+ (1|luna), data = fooofMRS )
summ(lm.model)


summary(lmerTest::lmer(Exponent ~ GABA_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ GABA_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

summary(lmerTest::lmer(Exponent ~ GABA_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ GABA_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ GABA_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zExp < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Exponent ~ GABA_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ GABA_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))

## 6.4 GABA VS Offset ----

lunaize(ggplot(data = fooofMRS , 
               aes(x = Offset, y = GABA_gamadj, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + geom_smooth(aes(group = Region, linetype=Region), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("GABA") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer(GABA_gamadj  ~ Offset + age + Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)


lm.model <- lmer(GABA_gamadj  ~ Offset * age + Condition +Region +   (1|luna), data = fooofMRS )
car::Anova(lm.model)
summ(lm.model)


summary(lmerTest::lmer(Offset ~ GABA_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GABA_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

summary(lmerTest::lmer(Offset ~ GABA_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GABA_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ GABA_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Offset ~ GABA_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ GABA_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))

## 6.5 Ratio VS Exponent ----

lunaize(ggplot(data = fooofMRS, 
               aes(x = Exponent, y = Ratio_gamadj, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = Region, linetype=Region, alpha = 0.1), method="lm", alpha = 0.8)  + 
          scale_color_manual(values=c("gold3", "blue4"))) +
  ylab("Glu/GABA Ratio") + xlab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Ratio_gamadj  ~ Exponent + age + Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)
car::Anova(lm.model)


lm.model <- lmer(Ratio_gamadj  ~ Exponent *age +Condition +Region+   (1|luna), data = fooofMRS )
summ(lm.model)
car::Anova(lm.model)


summary(lmerTest::lmer(Exponent ~ Ratio_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Ratio_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

summary(lmerTest::lmer(Exponent ~ Ratio_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ Ratio_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zExp < 2)))

# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ Ratio_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zExp < 2), random=list(luna=~1))

summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Exponent ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ Ratio_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))


## 6.6 Ratio VS Offset ----

lunaize(ggplot(data = fooofMRS, 
               aes(x = Offset, y = Ratio_gamadj, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = Region, linetype=Region, alpha = 0.1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  ylab("Glu/GABA Ratio") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')



lm.model <- lmer(Ratio_gamadj  ~  Offset + age + Condition  + Region + (1|luna), data = fooofMRS)
summ(lm.model)
car::Anova(lm.model)


lm.model <- lmer(Ratio_gamadj  ~ Offset *sex +age +Condition + Region + (1|luna), data = fooofMRS )
summ(lm.model)
car::Anova(lm.model)


summary(lmerTest::lmer(Offset ~ Ratio_gamadj + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ Ratio_gamadj + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

summary(lmerTest::lmer(Offset ~ Ratio_gamadj*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ Ratio_gamadj*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Ratio_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)


# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ Ratio_gamadj + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Offset ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ Ratio_gamadj + inverseAge + Condition + Region + (1|luna), data = fooofMRS))


## 6.7 Gaba glu imbalance VS offset ----

lunaize(ggplot(data = fooofMRS, 
               aes(x = Offset, y = GluGABAimbalanceABS, by =luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = Region, linetype=Region, alpha = 0.01), method="lm",alpha=.8,size=1)) + 
  scale_color_manual(values=c("gold3", "blue4")) + 
  ylab("Glu GABA Imbalance") + xlab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS + age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS + inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS*age + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))
summary(lmerTest::lmer(Offset ~ GluGABAimbalanceABS*inverseage + Region + Condition + (1|subjID), data = fooofMRS %>% filter(zOffset < 2) ))

# gam approach, with non-linear form of age
gam.model <- gamm(Offset ~ GluGABAimbalanceABS + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zOffset < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)


AIC(lmer(Offset ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Offset ~ GluGABAimbalanceABS + inverseAge + Condition + Region + (1|luna), data = fooofMRS))


## 6.8 Gaba glu imbalance VS exponent ----

lunaize(ggplot(data = fooofMRS %>% mutate(zExp = abs(scale(Exponent)[,1])) %>% filter(zExp < 2), 
               aes(x = GluGABAimbalance, y = Exponent, by =luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region, Condition)), alpha = 0.2) + geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1, linetype=Region, alpha = 0.01), method="lm",formula = y ~ poly(x,2),alpha=.8,size=1)) + 
  scale_color_manual(values=c("gold3", "blue4")) + 
  xlab("Glu GABA Imbalance") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')




# lmer approach, different forms of age, w/ and w/out interactions
summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + age + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2)))
summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + inverseage + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))

summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS*age + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))
summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS*inverseage + Region + Condition + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))

summary(lmerTest::lmer(Exponent ~ GluGABAimbalanceABS + inverseage + Condition + Region + (1|subjID), 
                       data = fooofMRS %>% filter(zExp < 2) ))


# gam approach, with non-linear form of age
gam.model <- gamm(Exponent ~ GluGABAimbalanceABS + s(age, k=3) + Region + Condition, 
                  data = fooofMRS %>% filter(zExp < 2), random=list(luna=~1))
summary(gam.model$gam)
print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)

AIC(lmer(Exponent ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = fooofMRS), 
    lmer(Exponent ~ GluGABAimbalanceABS + inverseAge + Condition + Region + (1|luna), data = fooofMRS))


# 7.0 effect size of imbalance-exponent on brain electrode map ----

fooofallChannels <- read.csv('~/scratch/fooof/allChannelsFOOOFMeasures_20230214.csv')
fooofallChannels <- fooofallChannels %>% separate(Subject,c("luna","vdate"), remove=F)

fooofallChannels_MRS <- merge(fooofallChannels, fooofMRSbehavior, by =c("luna", "visitno")) 

exponentImbalanceCoef <- data.frame()
for (chan in unique(fooofallChannels_MRS$Channel)) {
  ageRes <- lm(z~age + e + age + GMrat + Condition+ Region , data = fooofallChannels_MRS %>% filter(Channel == chan) %>% mutate(z = scale(GluGABAimbalanceABS)[,1],e = scale(Exponent)[,1] ,invage=1/age))   
  exponentImbalanceCoef <- rbind(exponentImbalanceCoef, data.frame(chan,ageRes$coefficients[2]))
}

exponentImbalanceCoef$labels <- exponentImbalanceCoef$chan
exponentImbalanceCoef <- merge(exponentImbalanceCoef, chanLocs, by = "labels")

lunaize(ggplot(exponentImbalanceCoef, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.03, low="white", mid="red", high="firebrick4") + ggtitle("Effect size of imbalance~exponent controlling for age, condition, region, and GMrat")  + theme(text = element_text(size = 30)))


offsetImbalanceCoef <- data.frame()
for (chan in unique(fooofallChannels_MRS$Channel)) {
  ageRes <- lm(z~age + e + age + GMrat + Condition+ Region , data = fooofallChannels_MRS %>% filter(Channel == chan) %>% mutate(z = scale(GluGABAimbalanceABS)[,1],e = scale(Offset)[,1] ,invage=1/age))   
  offsetImbalanceCoef <- rbind(offsetImbalanceCoef, data.frame(chan,ageRes$coefficients[2]))
}

offsetImbalanceCoef$labels <- offsetImbalanceCoef$chan
offsetImbalanceCoef <- merge(offsetImbalanceCoef, chanLocs, by = "labels")

lunaize(ggplot(offsetImbalanceCoef, aes(x = -Y, y = X, fill = -ageRes.coefficients.2., z = -ageRes.coefficients.2., label = labels)) + geom_topo(chan_markers = "text") + scale_fill_gradient2(midpoint=0.025, low="white", mid="red", high="firebrick4") + ggtitle("Effect size of imbalance~offset controlling for age, condition, region, and GMrat") + theme(text = element_text(size = 30))) 


# 8.0 Mediation ----

## 8.1 Glu gaba imbalance on Exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(GluGABAimbalanceABS, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on exponent (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on imbalnce  (a)
model.M <- lme4::lmer(GluGABAimbalanceABS ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS imbalance on exponent (b)
model.Y <- lme4::lmer(Exponent ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GluGABAimbalanceABS", boot = FALSE, sims = 1000)
(summary(results))


## 8.2 Glu gaba imbalance on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(GluGABAimbalanceABS, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on imbalance  (a)
model.M <- lme4::lmer(GluGABAimbalanceABS ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS imbalance on offset (b)
model.Y <- lme4::lmer(Offset ~ GluGABAimbalanceABS + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GluGABAimbalanceABS", boot = FALSE, sims = 1000)
(summary(results))


## 8.3 Glu gaba ratio on exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(Ratio_gamadj, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Ratio_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Exponent ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Ratio_gamadj", boot = FALSE, sims = 1000)
(summary(results))


## 8.4 Glu gaba ratio on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(Ratio_gamadj, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Ratio_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Offset ~ Ratio_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Ratio_gamadj", boot = FALSE, sims = 1000)
(summary(results))


## 8.5 Glu on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(Glu_gamadj, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on glu  (a)
model.M <- lme4::lmer(Glu_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS glu on offset (b)
model.Y <- lme4::lmer(Offset ~ Glu_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Glu_gamadj", boot = FALSE, sims = 1000)
(summary(results))

## 8.6 Glu on Exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(Glu_gamadj, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Glu_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Exponent ~ Glu_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Glu_gamadj", boot = FALSE, sims = 1000)
(summary(results))


## 8.7 GABA on Offset ----

mediationMatrix <- fooofMRS %>% dplyr::select(GABA_gamadj, Offset, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Offset ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(GABA_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Offset ~ GABA_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GABA_gamadj", boot = FALSE, sims = 1000)
(summary(results))

## 8.8 GABA on Exponent ----

mediationMatrix <- fooofMRS %>% dplyr::select(GABA_gamadj, Exponent, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(Exponent ~ age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(GABA_gamadj ~ age + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(Exponent ~ GABA_gamadj + age + Condition + Region + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "GABA_gamadj", boot = FALSE, sims = 1000)
(summary(results))

## 8.9 Glu on Latency ----

mediationMatrix <- fooofMRSbehavior %>% dplyr::select(Glu_gamadj, mgsLatency, Condition, Region, luna, visitno, age) 
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

#the effect of age on offset (c)
model.0 <- lme4::lmer(mgsLatency ~ age  + (1|luna), data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
summ(model.0)
#anova(model.00, model.0) # significant improvement in fit when you include age 

#the effect of age on Ratio  (a)
model.M <- lme4::lmer(Glu_gamadj ~ age  + Region + (1|luna) + Condition, data = mediationMatrix )
print(car::Anova(model.M))
print(summary(model.M))
summ(model.M)

#the effect of MRS ratio on offset (b)
model.Y <- lme4::lmer(mgsLatency  ~ Glu_gamadj +Region + age + (1|luna), data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))
summ(model.Y)


results <- mediation::mediate(model.M, model.Y, treat = "age", mediator = "Glu_gamadj", boot = FALSE, sims = 1000)
(summary(results))






# 9.0 FOOOF vs Beh ----


## 9.2 Compare FOOOF, MRS, & Behavior ----


### Accuracy vs Exponent ----

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

AIC(lmer(Exponent ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Accuracy vs Offset ----

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

AIC(lmer(Offset ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy vs Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = absBestError, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj ~ absBestError + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Glu_gamadj ~ absBestError* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy vs GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = absBestError, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab(" GABA ")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj ~ absBestError + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GABA_gamadj ~ absBestError* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





### Accuracy vs Imbalance ----

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

AIC(lmer(GluGABAimbalanceABS ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


###  Accuracy vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior, 
               aes(y = Ratio_gamadj, x = absBestError, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Ratio_gamadj~  absBestError+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  absBestError* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ absBestError + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ absBestError + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))




### Accuracy Var vs Exponent ----

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

AIC(lmer(Exponent ~ absBestError_sd+ age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ absBestError_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Accuracy Var vs Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = absBestError_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj ~ absBestError_sd + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Glu_gamadj ~ absBestError_sd* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy Var vs GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = absBestError_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab(" GABA ")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj ~ absBestError_sd + inverseAge + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GABA_gamadj ~ absBestError_sd* inverseAge  + Region+GMrat + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





### Accuracy Var vs Offset ----

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

AIC(lmer(Offset ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Accuracy Var  vs Imbalance ----

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

AIC(lmer(GluGABAimbalanceABS ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Accuracy Var vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = absBestError_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Accuracy Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Ratio_gamadj~  absBestError_sd+ ageage + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  absBestError_sd* age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ absBestError_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ absBestError_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Latency vs Exponent ----

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

AIC(lmer(Exponent ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


###  Latency  vs Offset ----

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

AIC(lmer(Offset ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Latency VS Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj~  mgsLatency+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Glu_gamadj~  mgsLatency* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


# visualizing age interactions 
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ 
  theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))

###  Latency VS GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("GABA")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj~  mgsLatency+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( GABA_gamadj~  mgsLatency* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


# visualizing age interactions 
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ 
  theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))



###  Latency vs Imbalance ----

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

AIC(lmer(GluGABAimbalanceABS ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

# visualizing age interactions
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
               aes(y = GluGABAimbalanceABS, x = mgsLatency, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ 
  theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))


lm.model <- lmer(GluGABAimbalanceABS ~  mgsLatency + Region + (1|luna), data = fooofMRSbehavior%>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30")))  %>% filter(ageGroup == "17-22"))
car::Anova(lm.model)
summ(lm.model)



###  Latency vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Ratio_gamadj~  mgsLatency+ age + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  mgsLatency* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ mgsLatency + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ mgsLatency + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


# visualizing age interactions 
lunaize(ggplot(data = fooofMRSbehavior %>% mutate(ageGroup = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30"))), 
               aes(y = Ratio_gamadj, x = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2)) + 
  geom_point(alpha=.5) + 
  geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.4) + 
  xlab("Latency") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ 
  theme(legend.position='none') + scale_color_manual(values=c("gold3", "blue4", "red4"))





### Latency Var  vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior, 
               aes(y = Exponent, x = mgsLatency_sd, by = luna, color = ageGroup)) + 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group =1), method="lm", alpha = 0.8)) + 
  scale_color_manual(values=c("gold3", "blue4")) +
  xlab("Latency Var") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Exponent~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Exponent~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Exponent ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





###  Latency Var vs Offset ----
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

AIC(lmer(Offset ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

###  Latency Var vs Glu ----
lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group =1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency Var") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( Glu_gamadj~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(Glu_gamadj ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





###  Latency Var vs GABA ----
lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group =1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Latency Var") + ylab("GABA")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GABA_gamadj~  mgsLatency_sd+ inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer(GABA_gamadj ~  mgsLatency_sd* inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))





### Latency Var vs Imbalance ----

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

AIC(lmer(GluGABAimbalanceABS ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Latency Var vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = mgsLatency_sd, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) +
  xlab("Latency Var") + ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer(Ratio_gamadj ~  mgsLatency_sd+ age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)


lm.model <- lmer( Ratio_gamadj~  mgsLatency_sd* age  + Region + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)


AIC(lmer(Ratio_gamadj ~ mgsLatency_sd + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ mgsLatency_sd + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))

### Spatial Span Max vs Exponent ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Exponent, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Exponent")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmerTest::lmer( Exponent ~  SSP_maxSpan + inverseAge  + Region + Condition + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summ(lm.model)

lm.model <- lmer( Exponent ~  SSP_maxSpan*age  + Region + Condition + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Exponent ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Exponent ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))




### Spatial Span Max vs Offset ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Offset, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Offset")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmerTest::lmer( Offset ~  SSP_maxSpan + inverseAge  + Region + Condition + (1|luna), data = fooofMRSbehavior)
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmerTest::lmer( Offset ~  SSP_maxSpan*inverseAge  + Region + Condition + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Offset ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Offset ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs Glu ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Glu_gamadj, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Glu")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer( Glu_gamadj ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior %>% filter(ageGroup == "Adult") )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( Glu_gamadj ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Glu_gamadj ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Glu_gamadj ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs GABA ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GABA_gamadj, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("GABA")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer( GABA_gamadj ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior  %>% filter(ageGroup == "Adol"))
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( GABA_gamadj ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GABA_gamadj ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GABA_gamadj ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs Ratio ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = Ratio_gamadj, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

lm.model <- lmer( Ratio_gamadj ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior  %>% filter(ageGroup == "Adol"))
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( Ratio_gamadj ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(Ratio_gamadj ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(Ratio_gamadj ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


### Spatial Span Max vs Imbalance ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(y = GluGABAimbalanceABS, x = SSP_maxSpan, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) + 
  xlab("Span Max") + ylab("Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


lm.model <- lmer( GluGABAimbalanceABS ~  SSP_maxSpan + age  + Region + (1|luna), data = fooofMRSbehavior  %>% filter(ageGroup == "Adol"))
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

lm.model <- lmer( GluGABAimbalanceABS ~  SSP_maxSpan*age  + Region  + (1|luna), data = fooofMRSbehavior )
car::Anova(lm.model)
summary(lm.model)
summ(lm.model)

AIC(lmer(GluGABAimbalanceABS ~ SSP_maxSpan + age + Condition + Region + (1|luna), data = fooofMRSbehavior), 
    lmer(GluGABAimbalanceABS ~ SSP_maxSpan + inverseAge + Condition + Region + (1|luna), data = fooofMRSbehavior))


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

## 9.4 Fooof loop ----

### Fooof vs Age ----

fooofVars <- c('Offset','Exponent')

output <- c()  
for (fooofVar in fooofVars) {
  
  model <- paste0(fooofVar, ' ~ ', 's(age, k = 3) + Region + Condition')
  model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofLong))
  
  F <- summary(model.out$gam)$s.table[1,3]
  p <- summary(model.out$gam)$s.table[1,4]
  RegionF <- summary(model.out$gam)$pTerms.table[1,2]
  Regionp <- summary(model.out$gam)$pTerms.table[1,3]
  ConditionF <- summary(model.out$gam)$pTerms.table[2,2]
  Conditionp <- summary(model.out$gam)$pTerms.table[2,3]
  
  output <- rbind(output, data.frame(fooofVar, F, p, RegionF, Regionp, ConditionF, Conditionp))
}
output  


## looking at interactions 
# Interaction of condition 
fooofLong$oCondition <- ordered(fooofLong$Condition, levels = c("eyesOpen","eyesClosed")) # eyes open will be the reference group

output <- c()  
for (fooofVar in fooofVars) {
  
  model_formula <- as.formula(paste0(fooofVar, ' ~ ', ' oCondition + s(age, k = 3, fx = T) + s(age, by = oCondition, k = 4, fx = T) + Region')) 
  
  model.out <- (mgcv::gamm(model_formula,
                           random = list(luna=~1),
                           data = fooofLong))
  
  ConditionF <- summary(model.out$gam)$s.table[2,3]
  Conditionp <- summary(model.out$gam)$s.table[2,4]
  
  output <- rbind(output, data.frame(fooofVar, ConditionF, Conditionp))
}
output  


# region interaction
fooofLong$oRegion <- ordered(fooofLong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC reference group

output <- c()  
for (fooofVar in fooofVars) {
  
  model_formula <- as.formula(paste0(fooofVar, ' ~ ', ' oRegion + s(age, k = 3, fx = T) + s(age, by = oRegion, k = 4, fx = T) + Condition')) 
  
  model.out <- (mgcv::gamm(model_formula,
                           random = list(luna=~1),
                           data = fooofLong))
  
  RegionF <- summary(model.out$gam)$s.table[2,3]
  Regionp <- summary(model.out$gam)$s.table[2,4]
  
  output <- rbind(output, data.frame(fooofVar, RegionF, Regionp))
}
output  



### Fooof vs Behavior ----

fooofVars <- c('Offset','Exponent')
behVars <- c('absBestError','absBestError_sd','mgsLatency', 'mgsLatency_sd', 'SSP_maxSpan')


output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    model <- paste0(fooofVar, ' ~ ', behVar, '+ inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))
    b <- model.out$coefficients[2,1]
    t <- model.out$coefficients[2,4]
    p <- model.out$coefficients[2,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(fooofVar, behVar, b, t, pcor, p))
  }
}  
output  

# Looking at age interaction
output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    model <- paste0(fooofVar, ' ~ ', behVar, '* inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))
    b <- model.out$coefficients[6,1]
    t <- model.out$coefficients[6,4]
    p <- model.out$coefficients[6,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(fooofVar, behVar, b, t, pcor, p))
  }
}  
output  

# looks at GAMs
output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    
    model <- paste0(fooofVar, ' ~ ', behVar, '+ s(age, k = 3) + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    b <- summary(model.out$gam)$p.table[2,1]
    t <- summary(model.out$gam)$p.table[2,3]
    p <- summary(model.out$gam)$p.table[2,4]
    regionb <- summary(model.out$gam)$p.table[3,1]
    regiont <- summary(model.out$gam)$p.table[3,3]
    regionp <- summary(model.out$gam)$p.table[3,4]
    condb <- summary(model.out$gam)$p.table[4,1]
    condt <- summary(model.out$gam)$p.table[4,3]
    condp <- summary(model.out$gam)$p.table[4,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    regionpcor <- p.adjust((regionp), method = "bonferroni", n = 4)
    condpcor <- p.adjust((condp), method = "bonferroni", n = 4)
    
    
    
    output <- rbind(output, data.frame(fooofVar, behVar, b, t, p, pcor, regionb, regiont, regionp, regionpcor, condb, condt, condp, condpcor))
  }
}
output  


# looks at age interactions
output <- c()  
for (fooofVar in fooofVars) {
  for (behVar in behVars) {
    
    model <- paste0(fooofVar, ' ~ ','s(age, k = 3, by = ', behVar,') + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    f <- summary(model.out$gam)$s.table[1,3]
    p <- summary(model.out$gam)$s.table[1,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(fooofVar, behVar, f, p, pcor))
  }
}
output  


## 9.5 mrsi loop ----

### MRSI vs Age ----
mrsiVars <- c('Ratio_gamadj','GluGABAimbalanceABS','Glu_gamadj','GABA_gamadj')

dim(MRSlong$Ratio_gamadj) <- NULL
dim(MRSlong$Glu_gamadj) <- NULL
dim(MRSlong$GABA_gamadj) <- NULL

output <- c()  
for (mrsiVar in mrsiVars) {
  
  model <- paste0(mrsiVar, ' ~ ', 's(age, k = 3) + Region')
  model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = MRSlong ))
  F <- summary(model.out$gam)$s.table[1,3]
  p <- summary(model.out$gam)$s.table[1,4]
  RegionF <- summary(model.out$gam)$pTerms.table[1,2]
  Regionp <- summary(model.out$gam)$pTerms.table[1,3]
  
  output <- rbind(output, data.frame(mrsiVar, F, p, RegionF, Regionp))
}
output  

# region interaction
MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC reference group

output <- c()  
for (mrsiVar in mrsiVars) {
  
  model_formula <- as.formula(paste0(mrsiVar, ' ~ ', ' oRegion + s(age, k = 3, fx = T) + s(age, by = oRegion, k = 3, fx = T)')) 
  
  model.out <- (mgcv::gamm(model_formula,
                           random = list(luna=~1),
                           data = MRSlong))
  
  RegionF <- summary(model.out$gam)$s.table[2,3]
  Regionp <- summary(model.out$gam)$s.table[2,4]
  
  output <- rbind(output, data.frame(mrsiVar, RegionF, Regionp))
}
output  


### MRS vs Behavior ----
#'absBestError','absBestError_sd','mgsLatency', 'mgsLatency_sd', 'Glu_gamadj', 'GABA_gamadj', 'GluGABAimbalanceABS'
mrsiVars <- c('Ratio_gamadj')
behVars <- c('SSP_maxSpan')

output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    model <- paste0(mrsiVar, ' ~ ', behVar, '+ inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))  
    b <- model.out$coefficients[2,1]
    t <- model.out$coefficients[2,4]
    p <- model.out$coefficients[2,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, b, t, pcor, p))
  }
}  
output  

# Looking at age interaction
output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    model <- paste0(mrsiVar, ' ~ ', behVar, '* age + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRSbehavior))
    b <- model.out$coefficients[6,1]
    t <- model.out$coefficients[6,4]
    p <- model.out$coefficients[6,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, b, t, pcor, p))
  }
}  
output  



# looks at GAMs
output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    
    model <- paste0(mrsiVar, ' ~ ', behVar, '+ s(age, k = 3) + Region ')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    b <- summary(model.out$gam)$p.table[2,1]
    t <- summary(model.out$gam)$p.table[2,3]
    p <- summary(model.out$gam)$p.table[2,4]
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, b, t, p, pcor))
  }
}
output  

# looks at age interactions
output <- c()  
for (mrsiVar in mrsiVars) {
  for (behVar in behVars) {
    
    model <- paste0(mrsiVar, ' ~ ','s(age, k = 3, by = ', behVar,') + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = fooofMRSbehavior))
    f <- summary(model.out$gam)$s.table[1,3]
    p <- summary(model.out$gam)$s.table[1,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 4)
    
    output <- rbind(output, data.frame(mrsiVar, behVar, f, p, pcor))
  }
}
output  


## 9.6 FOOOF vs MRS Loop ----
### Linear Models ----

mrsiVars <- c('Ratio_gamadj','GluGABAimbalanceABS','Glu_gamadj','GABA_gamadj')
fooofVars <- c('Exponent', 'Offset')

# main effect
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ', mrsiVar, '+ inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRS))  
    b <- model.out$coefficients[2,1]
    t <- model.out$coefficients[2,4]
    p <- model.out$coefficients[2,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    
    regb <- model.out$coefficients[4,1]
    regt <- model.out$coefficients[4,4]
    regp <- model.out$coefficients[4,5]
    regpcor <- p.adjust((regp), method = "bonferroni", n = 2)
    
    
    conb <- model.out$coefficients[5,1]
    cont <- model.out$coefficients[5,4]
    conp <- model.out$coefficients[5,5]
    conpcor <- p.adjust((conp), method = "bonferroni", n = 2)
    
    
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor,regb, regt, regp,regpcor, conb, cont, conp,conpcor))   
  }
}  
output  

# age interaction 
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ', mrsiVar, '* inverseAge + Region + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRS))  
    b <- model.out$coefficients[6,1]
    t <- model.out$coefficients[6,4]
    p <- model.out$coefficients[6,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor))
  }
}  
output 



# region interaction 
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ', mrsiVar, '* Region + inverseAge + Condition + (1|luna)')
    model.out <- summary(lmerTest::lmer(model, data = fooofMRS))  
    b <- model.out$coefficients[6,1]
    t <- model.out$coefficients[6,4]
    p <- model.out$coefficients[6,5]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor))
  }
}  
output 


### GAMs ----
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    model <- paste0(fooofVar, ' ~ ',mrsiVar, '+ s(age, k = 3) + Region + Condition')
    model.out <- (mgcv::gamm(as.formula(model), data = fooofMRS, random=list(luna=~1)))
    b <- summary(model.out$gam)$p.table[2,1]
    t <- summary(model.out$gam)$p.table[2,3]
    p <- summary(model.out$gam)$p.table[2,4]
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    
    regb <- summary(model.out$gam)$p.table[3,1]
    regt <- summary(model.out$gam)$p.table[3,3]
    regp <- summary(model.out$gam)$p.table[3,4]
    
    condb <- summary(model.out$gam)$p.table[4,1]
    condt <- summary(model.out$gam)$p.table[4,3]
    condp <- summary(model.out$gam)$p.table[4,4]
    
    regpcor <- p.adjust((regp), method = "bonferroni", n = 2)
    condpcor <- p.adjust((condp), method = "bonferroni", n = 2)
    
    output <- rbind(output, data.frame(fooofVar, mrsiVar, b, t, p, pcor,regb, regt, regp,regpcor, condb, condt, condp,condpcor))    }
}
output 

# look at interactions 
output <- c()  
for (mrsiVar in mrsiVars) {
  for (fooofVar in fooofVars) {
    
    model_formula <- as.formula(paste0(fooofVar, ' ~ ', ' s(age, k = 3, fx = T) + s(age, by = ',mrsiVar, ',k = 4, fx = T) + Region + Condition')) 
    
    model.out <- (mgcv::gamm(model_formula,
                             random=list(luna=~1),
                             data = fooofMRS))
    
    F <- summary(model.out$gam)$s.table[2,3]
    p <- summary(model.out$gam)$s.table[2,4]
    
    pcor <- p.adjust((p), method = "bonferroni", n = 2)
    
    output <- rbind(output, data.frame(fooofVar, mrsiVar, F, p, pcor))
  }
}
output  

## 9.7 imbalance-accuracy ----

lunaize(ggplot(data = fooofMRSbehavior  , 
               aes(x = Ratio_gamadj, y = mgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = ageGroup), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) +
  theme(text = element_text(size = 14))#+ theme(legend.position='none')



lunaize(ggplot(data = fooofMRSbehavior, 
               aes(x = Ratio_gamadj, y = vgsLatency, by = luna, color = ageGroup))+ 
          geom_line(aes(group=interaction(luna,Region,Condition)), alpha = 0.2) + 
          geom_point(alpha=.5) + 
          geom_smooth(aes(group = 1), method="lm", alpha = 0.8) + 
          scale_color_manual(values=c("gold3", "blue4"))) +
  theme(text = element_text(size = 14))#+ theme(legend.position='none')


## 10.0 Residualized correlations

behVars <- c('absBestError','mgsLatency', 'vgsLatency', 'absBestError_sd', 'mgsLatency_sd','latencyDiff')
mrsiVars <- c('Ratio_gamadj','GluGABAimbalance','Glu_gamadj','gabaRes.residuals')
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
