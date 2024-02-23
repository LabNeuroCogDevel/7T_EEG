
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
library(eegUtils)

# Load data sets ----

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/ChannelLocs.csv')%>% select(-type)
merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")

SNRallChans40 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_allChans_40Hz.csv')
SNRallChans40 <- SNRallChans40 %>% separate(Subject, c('lunaID','visitDate')) %>% merge(.,ageValues, by = c("lunaID", "visitDate"))
SNRallChans40$hertz <- 40
write.csv(SNRallChans40, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_allChans_40hz_withages.csv')

SNRallChans30 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_allChans_30Hz.csv')
SNRallChans30 <- SNRallChans30 %>% separate(Subject, c('lunaID','visitDate'))%>% merge(.,ageValues, by = c("lunaID", "visitDate"))
SNRallChans30$hertz <- 30

SNRallChans20 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_allChans_20Hz.csv')
SNRallChans20 <- SNRallChans20 %>% separate(Subject, c('lunaID','visitDate'))%>% merge(.,ageValues, by = c("lunaID", "visitDate"))
SNRallChans20$hertz <- 20

SNRallChans <- rbind(SNRallChans40, SNRallChans30) %>% rbind(., SNRallChans20)

# Define outlier functions ----
outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


fit_inverse_age<-function(thisdf, xvar, covars=''){
  
  modsum <- summary(lmerTest::lmer(as.formula(paste0(xvar, ' ~ invage ', covars, ' + (1|lunaID)')), 
                                   data = thisdf))
  int <- modsum$coefficients[1,1]
  t <- modsum$coefficients[2,4]
  p <- modsum$coefficients[2,5]
  
  return(data.frame(y0=int, p.value=p, t=-1*t))
}

fit_gam<-function(thisdf, xvar, covars=''){
  
  model <- paste0(xvar, ' ~ ',covars, 's(age, k = 3)')
  model.out <- (mgcv::gamm(as.formula(model), data = thisdf, random=list(lunaID=~1)))
  F <- summary(model.out$gam)$s.table[1,3]
  p <- summary(model.out$gam)$s.table[1,4]
  
  return(data.frame(p.value=p, F=F))
}

fit_gam_fooof<-function(thisdf, xvar, fooofvar){
 
  model <- paste0(fooofvar, ' ~ ',xvar, ' + s(age, k = 3) + Condition')
  model.out <- (mgcv::gamm(as.formula(model), data = thisdf, random=list(lunaID=~1)))
  Fvalue <- summary(model.out$gam)$p.table[2,1]
  p <- summary(model.out$gam)$p.table[2,4]
  
  return(data.frame(p.value=p, Fvalue=Fvalue))
}

fit_gam_mrs<-function(thisdf, xvar, mrsvar){
  
  model <- paste0(mrsvar, ' ~ ',xvar, ' + s(age, k = 3)')
  model.out <- (mgcv::gamm(as.formula(model), data = thisdf, random=list(lunaID=~1)))
  Fvalue <- summary(model.out$gam)$p.table[2,1]
  p <- summary(model.out$gam)$p.table[2,4]
  
  return(data.frame(p.value=p, Fvalue=Fvalue))
}


# Merge SNR with ages and Channel Info ----
SNRallChans <- merge(SNRallChans, chanLocs, by = "urchan")

# Create age groups and remove DC comp
SNRallChans <- SNRallChans %>% 
  mutate(ageGroup = as.factor(cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30'))))  %>% filter(freqs > 0) %>% filter(Total < 1000)

# SNR Outlier Detection ----
SNRallChans_outlier <- SNRallChans %>% group_by(urchan, freqs, hertz) %>%
  mutate(across(c("Evoked", "Induced", "Total"), naoutlier)) %>% ungroup()

SNRallChans_outlier$SNR <- log(SNRallChans_outlier$Evoked/SNRallChans_outlier$Induced)

SNRallChans_outlier <- SNRallChans_outlier %>% group_by(urchan, freqs) %>% mutate(across(c("SNR"), naoutlier)) %>% ungroup() %>% mutate(SNR = exp(SNR))
write.csv(SNRallChans_outlier, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_allChans_allfreqs.csv')


# SNR across freqs ----

# pivot to long format
SNRallChans_outlier_long <- pivot_longer(SNRallChans_outlier, cols = c("Total","Evoked", "SNR", "Induced"), names_to = 'measure', values_to = 'value')

# summarize by age group & measurement
SNRallChans_summary_byAge <- SNRallChans_outlier_long %>% 
  mutate(ageGroup = as.factor(cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))) %>%
  group_by(freqs, measure, ageGroup, hertz) %>%
  summarize(y = mean(value, na.rm=T),
            sd = sd(value, na.rm=T),
            n = n(),
            se = sd / sqrt(n))

# plot
colors_condition1 <- c("lightblue", "blue", "darkblue")
colors_condition2 <- c("lightcoral", "red", "darkred")
colors_condition3 <- c("lightgreen", "green", "darkgreen")

# evoked vs induced
ggplot(data = SNRallChans_summary_byAge %>% filter(freqs <= 50, measure %in% c('Evoked','Induced')), 
       aes(x = freqs, y = y, color = interaction(ageGroup, measure), 
           fill = interaction(ageGroup, measure), linetype=ageGroup)) +
  geom_line() +
  geom_ribbon(aes(ymin = y - se, ymax = y + se), alpha = 0.4, linewidth = 0) +
  scale_fill_manual(values = c(colors_condition1, colors_condition2)) +
  scale_color_manual(values = c(colors_condition1, colors_condition2)) +
  scale_y_log10() + 
  theme_bw() + facet_wrap(~ hertz)

# evoked vs total
ggplot(data = SNRallChans_summary_byAge %>% filter(freqs <= 50, measure %in% c('Evoked','Total')), 
       aes(x = freqs, y = y, color = interaction(ageGroup, measure), 
           fill = interaction(ageGroup, measure), linetype=ageGroup)) +
  geom_line() +
  geom_ribbon(aes(ymin = y - se, ymax = y + se), alpha = 0.4, linewidth = 0) +
  scale_fill_manual(values = c(colors_condition1, colors_condition2)) +
  scale_color_manual(values = c(colors_condition1, colors_condition2)) +
  scale_y_log10() + 
  theme_bw() + facet_wrap(~ hertz)

# SNR
ggplot(data = SNRallChans_summary_byAge %>% filter(freqs <= 50, measure %in% c('SNR')), 
       aes(x = freqs, y = y, color = interaction(ageGroup, measure), 
           fill = interaction(ageGroup, measure), linetype=ageGroup)) +
  geom_line() +
  geom_ribbon(aes(ymin = y - se, ymax = y + se), alpha = 0.4, linewidth = 0) +
  scale_fill_manual(values = c(colors_condition1)) +
  scale_color_manual(values = c(colors_condition1)) +
  scale_y_log10() + 
  theme_bw() + facet_wrap(~ hertz)

# Average SNR across all subjects ---- 
SNRaggregate <- aggregate(SNR~ labels+urchan +X + Y, data = SNRallChans_outlier %>% filter(freqs >= 38 & freqs <= 42), mean)

## SNR across topography ----
lunaize(ggplot(SNRaggregate, aes(x = -Y, y = X, fill = SNR, z = SNR, label = labels)) + geom_topo(chan_markers = "text") + 
          scale_fill_gradient2(low="blue", mid="white", high="red") + ggtitle("SNR")  + theme(text = element_text(size = 30)))

# SNR in age groups ----
lunaize(ggplot(SNRallChans_outlier %>% filter(freqs == 20 & hertz == 20), aes(x = -Y, y = X, fill = SNR, z = SNR, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          ggtitle("SNR, 40 hz")  + theme(text = element_text(size = 30))) + facet_wrap(~ageGroup)

# Evoked in age groups ----
lunaize(ggplot(SNRallChans_outlier %>% filter(freqs == 40), aes(x = -Y, y = X, fill = Evoked, z = Evoked, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          ggtitle("Evoked, 40 hz")  + theme(text = element_text(size = 30))) + facet_wrap(~ageGroup)

# Induced in age groups ----
lunaize(ggplot(SNRallChans_outlier %>% filter(freqs == 40), aes(x = -Y, y = X, fill = Induced, z = Induced, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(0,0.2)) + 
          ggtitle("Induced, 40 hz")  + theme(text = element_text(size = 30))) + facet_wrap(~ageGroup)

# SNR across age topo ----
## linear ----
SNRageEffects <- SNRallChans_outlier %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs, hertz) %>% 
  dplyr::do(fit_inverse_age(thisdf=., 'SNR'))

SNRageEffects <- merge(SNRageEffects, chanLocs, by = "urchan") 
SNRageEffects <- SNRageEffects %>% mutate(logp = -1*log10(p.value))

lunaize(ggplot(SNRageEffects %>% filter(freqs == 40), aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNRageEffects$logp)+.1)) + 
          ggtitle("SNR, 40 hz")  + theme(text = element_text(size = 30))) 

## gam ----
### 40-40 hz ----
SNRageEffects4040 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))


SNRageEffects4040 <- merge(SNRageEffects4040, chanLocs, by = "urchan") 
SNRageEffects4040 <- SNRageEffects4040 %>% mutate(logp = -1*log10(p.value))
SNRageEffects4040$Fvalue <- SNRageEffects4040$F
SNRageEffects4040$hertz <- '40 hz'


### 40-30 hz ----
SNRageEffects4030 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects4030 <- merge(SNRageEffects4030, chanLocs, by = "urchan") 
SNRageEffects4030 <- SNRageEffects4030 %>% mutate(logp = -1*log10(p.value))
SNRageEffects4030$Fvalue <- SNRageEffects4030$F
SNRageEffects4030$hertz <- '30 hz'


### 40-20 hz ----
SNRageEffects4020 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects4020 <- merge(SNRageEffects4020, chanLocs, by = "urchan") 
SNRageEffects4020 <- SNRageEffects4020 %>% mutate(logp = -1*log10(p.value))
SNRageEffects4020$Fvalue <- SNRageEffects4020$F
SNRageEffects4020$hertz <- '20 hz'

### 30-40hz ----
SNRageEffects3040 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects3040 <- merge(SNRageEffects3040, chanLocs, by = "urchan") 
SNRageEffects3040 <- SNRageEffects3040 %>% mutate(logp = -1*log10(p.value))
SNRageEffects3040$Fvalue <- SNRageEffects3040$F
SNRageEffects3040$hertz <- '40 hz'


### 30-30 hz ----
SNRageEffects3030 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects3030 <- merge(SNRageEffects3030, chanLocs, by = "urchan") 
SNRageEffects3030 <- SNRageEffects3030 %>% mutate(logp = -1*log10(p.value))
SNRageEffects3030$Fvalue <- SNRageEffects3030$F
SNRageEffects3030$hertz <- '30 hz'



### 30-20 hz ----
SNRageEffects3020 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects3020 <- merge(SNRageEffects3020, chanLocs, by = "urchan") 
SNRageEffects3020 <- SNRageEffects3020 %>% mutate(logp = -1*log10(p.value))
SNRageEffects3020$Fvalue <- SNRageEffects3020$F
SNRageEffects3020$hertz <- '20 hz'


### 20-40 hz ----
SNRageEffects2040 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects2040 <- merge(SNRageEffects2040, chanLocs, by = "urchan") 
SNRageEffects2040 <- SNRageEffects2040 %>% mutate(logp = -1*log10(p.value))
SNRageEffects2040$Fvalue <- SNRageEffects2040$F
SNRageEffects2040$hertz <- '40 hz'


### 20-30 hz ----
SNRageEffects2030 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects2030 <- merge(SNRageEffects2030, chanLocs, by = "urchan") 
SNRageEffects2030 <- SNRageEffects2030 %>% mutate(logp = -1*log10(p.value))
SNRageEffects2030$Fvalue <- SNRageEffects2030$F
SNRageEffects2030$hertz <- '30 hz'



### 20-20 hz ----
 SNRageEffects2020 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'SNR'))

SNRageEffects2020 <- merge(SNRageEffects2020, chanLocs, by = "urchan") 
SNRageEffects2020 <- SNRageEffects2020 %>% mutate(logp = -1*log10(p.value))
SNRageEffects2020$Fvalue <- SNRageEffects2020$F
SNRageEffects2020$hertz <- '20 hz'


### Combine all combinations 

SNRageEffectsAll <- rbind(SNRageEffects4020, SNRageEffects4030) %>% rbind(., SNRageEffects4040) %>% rbind(., SNRageEffects3040) %>% rbind(., SNRageEffects3030) %>% rbind(., SNRageEffects3020) %>% 
  rbind(., SNRageEffects2040) %>% rbind(., SNRageEffects2030) %>% rbind(., SNRageEffects2020)

lunaize(ggplot(SNRageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text", interpolate = F, interp_limit = "head") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR")

lunaize(ggplot(SNRageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text", interpolate = F, interp_limit = "head") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNRageEffectsAll$logp)+.5)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR")

# Evoked across age topo ----
## linear ----
EvokedageEffects <- SNRallChans_outlier %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_inverse_age(thisdf=., 'Evoked'))

EvokedageEffects <- merge(EvokedageEffects, chanLocs, by = "urchan") 
EvokedageEffects <- EvokedageEffects %>% mutate(logp = -1*log10(p.value))

lunaize(ggplot(EvokedageEffects %>% filter(freqs == 40), aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          ggtitle("Evoked, 40 hz")  + theme(text = element_text(size = 30))) 

## gam ----
### 40-40 hz ----
EvokedageEffects4040 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))


EvokedageEffects4040 <- merge(EvokedageEffects40, chanLocs, by = "urchan") 
EvokedageEffects4040 <- EvokedageEffects40 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects4040$Fvalue <- EvokedageEffects4040$F
EvokedageEffects4040$hertz <- '40 hz'


### 40-30 hz ----
EvokedageEffects4030 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects4030 <- merge(EvokedageEffects4030, chanLocs, by = "urchan") 
EvokedageEffects4030 <- EvokedageEffects4030 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects4030$Fvalue <- EvokedageEffects4030$F
EvokedageEffects4030$hertz <- '30 hz'


### 40-20 hz ----
EvokedageEffects4020 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects4020 <- merge(EvokedageEffects4020, chanLocs, by = "urchan") 
EvokedageEffects4020 <- EvokedageEffects4020 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects4020$Fvalue <- EvokedageEffects4020$F
EvokedageEffects4020$hertz <- '20 hz'

### 30-40hz ----
EvokedageEffects3040 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects3040 <- merge(EvokedageEffects3040, chanLocs, by = "urchan") 
EvokedageEffects3040 <- EvokedageEffects3040 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects3040$Fvalue <- EvokedageEffects3040$F
EvokedageEffects3040$hertz <- '40 hz'


### 30-30 hz ----
EvokedageEffects3030 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects3030 <- merge(EvokedageEffects3030, chanLocs, by = "urchan") 
EvokedageEffects3030 <- EvokedageEffects3030 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects3030$Fvalue <- EvokedageEffects3030$F
EvokedageEffects3030$hertz <- '30 hz'



### 30-20 hz ----
EvokedageEffects3020 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects3020 <- merge(EvokedageEffects3020, chanLocs, by = "urchan") 
EvokedageEffects3020 <- EvokedageEffects3020 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects3020$Fvalue <- EvokedageEffects3020$F
EvokedageEffects3020$hertz <- '20 hz'


### 20-40 hz ----
EvokedageEffects2040 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects2040 <- merge(EvokedageEffects2040, chanLocs, by = "urchan") 
EvokedageEffects2040 <- EvokedageEffects2040 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects2040$Fvalue <- EvokedageEffects2040$F
EvokedageEffects2040$hertz <- '40 hz'


### 20-30 hz ----
EvokedageEffects2030 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects2030 <- merge(EvokedageEffects2030, chanLocs, by = "urchan") 
EvokedageEffects2030 <- EvokedageEffects2030 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects2030$Fvalue <- EvokedageEffects2030$F
EvokedageEffects2030$hertz <- '30 hz'



### 20-20 hz ----
EvokedageEffects2020 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Evoked'))

EvokedageEffects2020 <- merge(EvokedageEffects2020, chanLocs, by = "urchan") 
EvokedageEffects2020 <- EvokedageEffects2020 %>% mutate(logp = -1*log10(p.value))
EvokedageEffects2020$Fvalue <- EvokedageEffects2020$F
EvokedageEffects2020$hertz <- '20 hz'


### Combine all combinations 

EvokedageEffectsAll <- rbind(EvokedageEffects4020, EvokedageEffects4030) %>% rbind(., EvokedageEffects4040) %>% rbind(., EvokedageEffects3040) %>% rbind(., EvokedageEffects3030) %>% rbind(., EvokedageEffects3020) %>% 
  rbind(., EvokedageEffects2040) %>% rbind(., EvokedageEffects2030) %>% rbind(., EvokedageEffects2020)

lunaize(ggplot(EvokedageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked Activity")

lunaize(ggplot(EvokedageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(EvokedageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz)+ ggtitle("Evoked Activity")


# Induced across age topo ----
## linear ----
InducedageEffects <- SNRallChans_outlier %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_inverse_age(thisdf=., 'Induced'))

InducedageEffects <- merge(InducedageEffects, chanLocs, by = "urchan") 
InducedageEffects <- InducedageEffects %>% mutate(logp = -1*log10(p.value))

lunaize(ggplot(InducedageEffects %>% filter(freqs == 40), aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          ggtitle("Induced, 40 hz")  + theme(text = element_text(size = 30))) 

## gam ----
### 40-40 hz ----
InducedageEffects4040 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))


InducedageEffects4040 <- merge(InducedageEffects40, chanLocs, by = "urchan") 
InducedageEffects4040 <- InducedageEffects40 %>% mutate(logp = -1*log10(p.value))
InducedageEffects4040$Fvalue <- InducedageEffects4040$F
InducedageEffects4040$hertz <- '40 hz'


### 40-30 hz ----
InducedageEffects4030 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects4030 <- merge(InducedageEffects4030, chanLocs, by = "urchan") 
InducedageEffects4030 <- InducedageEffects4030 %>% mutate(logp = -1*log10(p.value))
InducedageEffects4030$Fvalue <- InducedageEffects4030$F
InducedageEffects4030$hertz <- '30 hz'


### 40-20 hz ----
InducedageEffects4020 <- SNRallChans_outlier %>% filter(freqs == 40 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects4020 <- merge(InducedageEffects4020, chanLocs, by = "urchan") 
InducedageEffects4020 <- InducedageEffects4020 %>% mutate(logp = -1*log10(p.value))
InducedageEffects4020$Fvalue <- InducedageEffects4020$F
InducedageEffects4020$hertz <- '20 hz'

### 30-40hz ----
InducedageEffects3040 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects3040 <- merge(InducedageEffects3040, chanLocs, by = "urchan") 
InducedageEffects3040 <- InducedageEffects3040 %>% filter(p.value > 0) %>% mutate(logp = -1*log10(p.value))
InducedageEffects3040$Fvalue <- InducedageEffects3040$F
InducedageEffects3040$hertz <- '40 hz'


### 30-30 hz ----
InducedageEffects3030 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects3030 <- merge(InducedageEffects3030, chanLocs, by = "urchan") 
InducedageEffects3030 <- InducedageEffects3030 %>% mutate(logp = -1*log10(p.value))
InducedageEffects3030$Fvalue <- InducedageEffects3030$F
InducedageEffects3030$hertz <- '30 hz'



### 30-20 hz ----
InducedageEffects3020 <- SNRallChans_outlier %>% filter(freqs == 30 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects3020 <- merge(InducedageEffects3020, chanLocs, by = "urchan") 
InducedageEffects3020 <- InducedageEffects3020 %>% mutate(logp = -1*log10(p.value))
InducedageEffects3020$Fvalue <- InducedageEffects3020$F
InducedageEffects3020$hertz <- '20 hz'


### 20-40 hz ----
InducedageEffects2040 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 40) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects2040 <- merge(InducedageEffects2040, chanLocs, by = "urchan") 
InducedageEffects2040 <- InducedageEffects2040 %>% mutate(logp = -1*log10(p.value))
InducedageEffects2040$Fvalue <- InducedageEffects2040$F
InducedageEffects2040$hertz <- '40 hz'


### 20-30 hz ----
InducedageEffects2030 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 30) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects2030 <- merge(InducedageEffects2030, chanLocs, by = "urchan") 
InducedageEffects2030 <- InducedageEffects2030 %>% mutate(logp = -1*log10(p.value))
InducedageEffects2030$Fvalue <- InducedageEffects2030$F
InducedageEffects2030$hertz <- '30 hz'



### 20-20 hz ----
InducedageEffects2020 <- SNRallChans_outlier %>% filter(freqs == 20 & hertz == 20) %>% 
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam(thisdf=., 'Induced'))

InducedageEffects2020 <- merge(InducedageEffects2020, chanLocs, by = "urchan") 
InducedageEffects2020 <- InducedageEffects2020 %>% mutate(logp = -1*log10(p.value))
InducedageEffects2020$Fvalue <- InducedageEffects2020$F
InducedageEffects2020$hertz <- '20 hz'


### Combine all combinations 

InducedageEffectsAll <- rbind(InducedageEffects4020, InducedageEffects4030) %>% rbind(., InducedageEffects4040) %>% rbind(., InducedageEffects3040) %>% rbind(., InducedageEffects3030) %>% rbind(., InducedageEffects3020) %>% 
  rbind(., InducedageEffects2040) %>% rbind(., InducedageEffects2030) %>% rbind(., InducedageEffects2020)

lunaize(ggplot(InducedageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced Activity")

lunaize(ggplot(InducedageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(InducedageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz)+ ggtitle("Induced Activity")


# SNR across age ----

ggplot(data = SNRallChans_outlier %>% filter(freqs == 40 & hertz == 40), aes(x = age, y = SNR, color = labels)) +
  stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha = 0.02) + 
  xlab("Age") + ylab("SNR") + ggtitle('SNR, 40 hz response at 40 hz clicks')

ggplot(data = SNRallChans_outlier %>% filter(freqs == 30 & hertz == 30), aes(x = age, y = SNR, color = labels)) +
  stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha = 0.02) + 
  xlab("Age") + ylab("SNR") + ggtitle('SNR, 30 hz response at 30 hz clicks')

ggplot(data = SNRallChans_outlier %>% filter(freqs == 20 & hertz == 20), aes(x = age, y = SNR, color = labels)) +
  stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha = 0.02) + 
  xlab("Age") + ylab("SNR") + ggtitle('SNR, 20 hz response at 20 hz clicks')

## CPz ----

ggplot(data = SNRallChans_outlier %>% filter(freqs == 40 | freqs == 30 | freqs == 20 & labels == 'CPz'), aes(x = age, y = SNR)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("SNR") + ggtitle('SNR at CPz')+ facet_wrap(~freqs+hertz)

gam.model <-  gamm(SNR ~ s(age, k = 3), data = SNRallChans_outlier  %>% filter(freqs == 40 & labels == 'CPz'), random=list(lunaID=~1))
summary(gam.model$gam)

## Cz ----

ggplot(data = SNRallChans_outlier %>% filter(freqs == hertz & labels == 'Cz'), aes(x = age, y = SNR)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("SNR") + ggtitle('SNR at Cz')+ facet_wrap(~freqs+hertz)

gam.model <-  gamm(SNR ~ s(age, k = 3), data = SNRallChans_outlier  %>%
                     filter(freqs == 40 & labels == 'Cz'), random=list(lunaID=~1))

summary(gam.model$gam)

## FCz ----

SNRallChans_outlier_FCz <- SNRallChans_outlier %>% filter(freqs == 40 & labels == 'FCz') %>% 
  mutate(across(c("SNR"), naoutlier)) %>% ungroup() %>% mutate(SNR = exp(SNR))


ggplot(data = SNRallChans_outlier_FCz %>% filter(freqs == 40 & labels == 'FCz'), aes(x = age, y = SNR)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("SNR") + ggtitle('SNR at FCz')

gam.model <-  gamm(SNR ~ s(age, k = 3), data = SNRallChans_outlier_FCz  %>% filter(freqs == 40 & labels == 'FCz'), random=list(lunaID=~1))
summary(gam.model$gam)

# Evoked across age ----
## CPz ----

ggplot(data = SNRallChans_outlier %>% filter(freqs == hertz & labels == 'CPz'), aes(x = age, y = Evoked)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("Evoked") + ggtitle('Evoked at CPz') + facet_wrap(~hertz)

gam.model <-  gamm(Evoked ~ s(age, k = 3), data = SNRallChans_outlier  %>% filter(freqs == 40 & labels == 'CPz'), random=list(lunaID=~1))
summary(gam.model$gam)


## Cz ----

ggplot(data = SNRallChans_outlier %>% filter(freqs == hertz & labels == 'Cz'), aes(x = age, y = Evoked)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("Evoked") + ggtitle('Evoked at Cz')+ facet_wrap(~hertz)

gam.model <-  gamm(Evoked ~ s(age, k = 3), data = SNRallChans_outlier  %>%
                     filter(freqs == 40 & labels == 'Cz'), random=list(lunaID=~1))

summary(gam.model$gam)

## FCz ----

SNRallChans_outlier_FCz <- SNRallChans_outlier %>% filter(freqs == 40 & labels == 'FCz') %>% 
  mutate(across(c("Evoked"), naoutlier)) %>% ungroup()


ggplot(data = SNRallChans_outlier_FCz %>% filter(freqs == 40 & labels == 'FCz'), aes(x = age, y = Evoked)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("Evoked") + ggtitle('Evoked at FCz')

gam.model <-  gamm(Evoked ~ s(age, k = 3), data = SNRallChans_outlier_FCz  %>% filter(freqs == 40 & labels == 'FCz'), random=list(lunaID=~1))
summary(gam.model$gam)

# Induced across age ----
## CPz ----

ggplot(data = SNRallChans_outlier %>% filter(freqs == hertz & labels == 'CPz'), aes(x = age, y = Induced)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("Induced") + ggtitle('Induced at CPz') + facet_wrap(~hertz)

gam.model <-  gamm(Induced ~ s(age, k = 3), data = SNRallChans_outlier  %>% filter(freqs == 40 & labels == 'CPz'), random=list(lunaID=~1))
summary(gam.model$gam)

## Cz ----

ggplot(data = SNRallChans_outlier %>% filter(freqs == hertz & labels == 'Cz'), aes(x = age, y = Induced)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("Induced") + ggtitle('Induced at Cz') + facet_wrap(~hertz)

gam.model <-  gamm(Induced ~ s(age, k = 3), data = SNRallChans_outlier  %>%
                     filter(freqs == 40 & labels == 'Cz'), random=list(lunaID=~1))

summary(gam.model$gam)

## FCz ----

SNRallChans_outlier_FCz <- SNRallChans_outlier %>% filter(freqs == hertz & labels == 'FCz') %>% 
  mutate(across(c("Induced"), naoutlier)) %>% ungroup()

ggplot(data = SNRallChans_outlier %>% filter(freqs == hertz & labels == 'FCz'), aes(x = age, y = Induced)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Age") + ylab("Induced") + ggtitle('Induced at FCz')+ facet_wrap(~hertz)

gam.model <-  gamm(Induced ~ s(age, k = 3), data = SNRallChans_outlier_FCz %>% 
                     filter(freqs == 40 & labels == 'FCz'), random=list(lunaID=~1))
summary(gam.model$gam)


# Load in FOOOF all channels ----
fooof <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/FOOOF/Results/allSubjectsFooofMeasures_20230516.csv') %>% select(-X) %>% separate(Subject, c('lunaID','visitDate'))
names(fooof)[names(fooof) == "Channel"] <- "labels"

# FOOOF Outlier Detection ----
fooof_outlier <- fooof %>% group_by(labels, Condition) %>%
  mutate(across(c("Exponent", "Offset"), naoutlier)) %>% ungroup()

# Merge SNR and FOOOF ----

SNR_fooof <- merge(SNRallChans_outlier, fooof_outlier, by = c("lunaID", "visitDate", "labels"))
write.csv(SNR_fooof, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_FOOOF.csv')


# SNR vs Exponent  ----
## Across Electrodes ----

### 40 hz ----
SNR_fooofEffects40 <- SNR_fooof %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'SNR','Exponent'))

SNR_fooofEffects40 <- merge(SNR_fooofEffects40, chanLocs, by = "urchan") 
SNR_fooofEffects40 <- SNR_fooofEffects40 %>% mutate(logp = -1*log10(p.value))

SNRfooof40 <- lunaize(ggplot(SNR_fooofEffects40, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                  limits = c(-1*log10(.05),NA)) + 
          ggtitle("SNR vs Exponent, 40 hz")  + theme(text = element_text(size = 30))) 

### 30 hz ----
SNR_fooofEffects30 <- SNR_fooof %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'SNR','Exponent'))

SNR_fooofEffects30 <- merge(SNR_fooofEffects30, chanLocs, by = "urchan") 
SNR_fooofEffects30 <- SNR_fooofEffects30 %>% mutate(logp = -1*log10(p.value))

SNRfooof30 <-lunaize(ggplot(SNR_fooofEffects30, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                  limits = c(-1*log10(.05),NA)) + 
          ggtitle("SNR vs Exponent, 30 hz")  + theme(text = element_text(size = 30))) 


### 20 hz ----
SNR_fooofEffects20 <- SNR_fooof %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'SNR','Exponent'))

SNR_fooofEffects20 <- merge(SNR_fooofEffects20, chanLocs, by = "urchan") 
SNR_fooofEffects20 <- SNR_fooofEffects20 %>% mutate(logp = -1*log10(p.value))

SNRfooof20 <-lunaize(ggplot(SNR_fooofEffects20, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                       geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                               limits = c(-1*log10(.05),NA)) + 
                       ggtitle("SNR vs Exponent, 20 hz")  + theme(text = element_text(size = 30))) 

cowplot::plot_grid(SNRfooof40, SNRfooof30, SNRfooof20)

## Cz ----
### gam ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'Cz'), aes(x = SNR, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("SNR") + ylab("Exponent")+ ggtitle('SNR vs Exponent at Cz')

gam.model <-  gamm(Exponent ~ SNR + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'Cz'), random=list(lunaID=~1))
summary(gam.model$gam)

### linear ----
#### 40 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'Cz'), aes(x = SNR, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("SNR") + ylab("Exponent")+ ggtitle('SNR vs Exponent at Cz')

summary(lmerTest::lmer(Exponent ~ SNR + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 40 & labels == 'Cz')))

#### 30 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 30 & labels == 'Cz'), aes(x = SNR, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("SNR") + ylab("Exponent")+ ggtitle('SNR vs Exponent at Cz')

summary(lmerTest::lmer(Exponent ~ SNR + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 30 & labels == 'Cz')))

#### 20 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 20 & labels == 'Cz'), aes(x = SNR, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("SNR") + ylab("Exponent")+ ggtitle('SNR vs Exponent at Cz')

summary(lmerTest::lmer(Exponent ~ SNR + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 20 & labels == 'Cz')))


## FC2 ----

ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'FC2'), aes(x = SNR, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("SNR") + ylab("Exponent")+ ggtitle('SNR vs Exponent at FC2')

gam.model <-  gamm(Exponent ~ SNR + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'FC2'), random=list(lunaID=~1))
summary(gam.model$gam)



# Evoked vs Exponent  ----
## Across Electrodes ----

### 40 hz ----
Exp_fooofEffects40 <- SNR_fooof %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Evoked','Exponent'))

Exp_fooofEffects40 <- merge(Exp_fooofEffects40, chanLocs, by = "urchan") 
Exp_fooofEffects40 <- Exp_fooofEffects40 %>% mutate(logp = -1*log10(p.value))

evkExp40 <- lunaize(ggplot(Exp_fooofEffects40, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                  limits = c(-1*log10(.05),NA)) + 
          ggtitle("Evoked vs Exponent, 40 hz")  + theme(text = element_text(size = 30))) 

### 30 hz ----
Exp_fooofEffects30 <- SNR_fooof %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Evoked','Exponent'))

Exp_fooofEffects30 <- merge(Exp_fooofEffects30, chanLocs, by = "urchan") 
Exp_fooofEffects30 <- Exp_fooofEffects30 %>% mutate(logp = -1*log10(p.value))

evkExp30 <- lunaize(ggplot(Exp_fooofEffects30, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                              limits = c(-1*log10(.05),NA)) + 
                      ggtitle("Evoked vs Exponent, 30 hz")  + theme(text = element_text(size = 30))) 


### 20 hz ----
Exp_fooofEffects20 <- SNR_fooof %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Evoked','Exponent'))

Exp_fooofEffects20 <- merge(Exp_fooofEffects20, chanLocs, by = "urchan") 
Exp_fooofEffects20 <- Exp_fooofEffects20 %>% mutate(logp = -1*log10(p.value))

evkExp20 <- lunaize(ggplot(Exp_fooofEffects20, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                              limits = c(-1*log10(.05),NA)) + 
                      ggtitle("Evoked vs Exponent, 20 hz")  + theme(text = element_text(size = 30))) 

cowplot::plot_grid(evkExp40, evkExp30, evkExp20)


## CPz ----
### gam ----

SNRFOOOFfiltered <- SNR_fooof %>% filter(freqs == 40 & labels == 'CPz' & !is.na(Exponent) & !is.na(age))

gamExp <- gamm(Exponent ~ s(age, k=3) + Condition, data = SNRFOOOFfiltered)
SNRFOOOFfiltered$ExpRes <- residuals(gamExp$gam)


ggplot(data = SNRFOOOFfiltered, aes(y = Evoked, x = ExpRes, color = Condition)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  ylab("Evoked") + xlab("Exponent Res")+ ggtitle('Evoked vs Exponent Res at CPz')

gam.model <-  gamm(Exponent ~ Evoked + s(age, k = 3) + Condition, data = SNRFOOOFfiltered, random=list(lunaID=~1))
summary(gam.model$gam)

summary(lmerTest::lmer(Exponent ~ Evoked + age + Condition + (1|lunaID), data = SNRFOOOFfiltered))
summary(lmerTest::lmer(ExpRes ~ Evoked + age + Condition + (1|lunaID), data = SNRFOOOFfiltered))


print(plot(getViz(gam.model$gam), allTerms = T), pages = 1)


ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'CPz'), aes(x = Evoked, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Evoked") + ylab("Exponent")+ ggtitle('Evoked vs Exponent at CPz')

gam.model <-  gamm(Exponent ~ Evoked + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'CPz'), random=list(lunaID=~1))
summary(gam.model$gam)

### linear ----
#### 40 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & hertz == 40 & labels == 'CPz'), aes(x = Evoked, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Evoked") + ylab("Exponent")+ ggtitle('Evoked vs Exponent at CPz')

summary(lmerTest::lmer(Exponent ~ Evoked + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 40 & hertz == 40 & labels == 'CPz')))

#### 30 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 30 & hertz == 30 & labels == 'CPz'), aes(x = Evoked, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Evoked") + ylab("Exponent")+ ggtitle('Evoked vs Exponent at CPz')

summary(lmerTest::lmer(Exponent ~ Evoked + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 30 & hertz == 30 & labels == 'CPz')))

#### 20 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 20 & hertz == 20 & labels == 'CPz'), aes(x = Evoked, y = Exponent)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Evoked") + ylab("Exponent")+ ggtitle('Evoked vs Exponent at CPz')

summary(lmerTest::lmer(Exponent ~ Evoked + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 20 & hertz == 20 & labels == 'CPz')))



# Induced vs Exponent ----
## Across Electrodes ----

### 40 hz ----
indExp_Effects40 <- SNR_fooof %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Induced','Exponent'))

indExp_Effects40 <- merge(indExp_Effects40, chanLocs, by = "urchan") 
indExp_Effects40 <- indExp_Effects40 %>% mutate(logp = -1*log10(p.value))

indExp40 <- lunaize(ggplot(indExp_Effects40, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                  limits = c(-1*log10(.05),NA)) + 
          ggtitle("Induced vs Exponent, 40 hz")  + theme(text = element_text(size = 30))) 

### 30 hz ----
indExp_Effects30 <- SNR_fooof %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Induced','Exponent'))

indExp_Effects30 <- merge(indExp_Effects30, chanLocs, by = "urchan") 
indExp_Effects30 <- indExp_Effects30 %>% mutate(logp = -1*log10(p.value))

indExp30 <- lunaize(ggplot(indExp_Effects30, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                              limits = c(-1*log10(.05),NA)) + 
                      ggtitle("Induced vs Exponent, 30 hz")  + theme(text = element_text(size = 30))) 

### 20 hz ----
  indExp_Effects20 <- SNR_fooof %>% filter(freqs == 20 & hertz == 20) %>%
    mutate(invage = 1/age) %>% 
    group_by(urchan, freqs) %>% 
    dplyr::do(fit_gam_fooof(thisdf=., 'Induced','Exponent'))
  
  indExp_Effects20 <- merge(indExp_Effects20, chanLocs, by = "urchan") 
  indExp_Effects20 <- indExp_Effects20 %>% mutate(logp = -1*log10(p.value))
  
  indExp20 <- lunaize(ggplot(indExp_Effects20, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                        geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", 
                                                                                limits = c(-1*log10(.05),NA)) + 
                        ggtitle("Induced vs Exponent, 20 hz")  + theme(text = element_text(size = 30))) 
  
  cowplot::plot_grid(indExp20, indExp30, indExp40)
  
  ## F4 ----
  ### gam ----
  ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'F4'), aes(y = Exponent, x = Induced)) +
    geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
    geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
    ylab("Exponent") + xlab("Induced")+ ggtitle('Exponent vs Induced at F4')
  
  gam.model <-  gamm(Exponent ~ Induced + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                       filter(freqs == 40 & labels == 'F4'), random=list(lunaID=~1))
  summary(gam.model$gam)
  
  ### linear ----
  ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'F4'), aes(y = Exponent, x = Induced)) +
    geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
    geom_point() + stat_smooth(method=lm) + 
    ylab("Exponent") + xlab("Induced")+ ggtitle('Exponent vs Induced at F4')
  
  summary(lmerTest::lmer(Exponent ~ Induced + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 40 & labels == 'F4')))
  
  
  ## FC4 ----
  ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'FC4'), aes(y = Exponent, x = Induced)) +
    geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
    geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
    ylab("Exponent") + xlab("Induced")+ ggtitle('Exponent vs Induced at FC4')
  
  gam.model <-  gamm(Exponent ~ Induced + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                       filter(freqs == 40 & labels == 'FC4'), random=list(lunaID=~1))
  summary(gam.model$gam)
  
# SNR vs Offset ----
## Across Electrodes ----
  
### 40 hz ----
SNR_offEffects40 <- SNR_fooof %>% filter(freqs == 40 & hertz ==40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'SNR','Offset'))

SNR_offEffects40 <- merge(SNR_offEffects40, chanLocs, by = "urchan") 
SNR_offEffects40 <- SNR_offEffects40 %>% mutate(logp = -1*log10(p.value))

SNRoff40 <- lunaize(ggplot(SNR_offEffects40, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
          ggtitle("SNR vs Offset, 40 hz")  + theme(text = element_text(size = 30))) 

### 30 hz ----
SNR_offEffects30 <- SNR_fooof %>% filter(freqs == 30 & hertz ==30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'SNR','Offset'))

SNR_offEffects30 <- merge(SNR_offEffects30, chanLocs, by = "urchan") 
SNR_offEffects30 <- SNR_offEffects30 %>% mutate(logp = -1*log10(p.value))

SNRoff30 <- lunaize(ggplot(SNR_offEffects30, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
                      ggtitle("SNR vs Offset, 30 hz")  + theme(text = element_text(size = 30))) 

### 20 hz ----
SNR_offEffects20 <- SNR_fooof %>% filter(freqs == 20 & hertz ==20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'SNR','Offset'))

SNR_offEffects20 <- merge(SNR_offEffects20, chanLocs, by = "urchan") 
SNR_offEffects20 <- SNR_offEffects20 %>% mutate(logp = -1*log10(p.value))

SNRoff20 <- lunaize(ggplot(SNR_offEffects20, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
                      ggtitle("SNR vs Offset, 20 hz")  + theme(text = element_text(size = 30))) 

cowplot::plot_grid(SNRoff20, SNRoff30, SNRoff40)

## CPz ----
### gam ----

ggplot(data =  SNR_fooof %>% filter(freqs == 40 & labels == 'CPz'), aes(x = SNR, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("SNR") + ylab("Offset")+ ggtitle('SNR vs Offset at CPz')

gam.model <-  gamm(Offset ~ SNR + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'CPz'), random=list(lunaID=~1))
summary(gam.model$gam)

### linear ----

ggplot(data =  SNR_fooof %>% filter(freqs == 40 & labels == 'CPz'), aes(x = SNR, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("SNR") + ylab("Offset")+ ggtitle('SNR vs Offset at CPz')

summary(lmerTest::lmer(Offset ~ SNR + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 40 & labels == 'CPz')))


## C1 ----

ggplot(data =  SNR_fooof %>% filter(freqs == 40 & labels == 'C1'), aes(x = SNR, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = , fx = T)) + 
  xlab("SNR") + ylab("Offset")+ ggtitle('SNR vs Offset at C1')

gam.model <-  gamm(Offset ~ SNR + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'C1'), random=list(lunaID=~1))
summary(gam.model$gam)


# Evoked vs Offset  ----
## Across Electrodes ----

### 40 hz ----
evkOffEffects40 <- SNR_fooof %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Evoked','Offset'))

evkOffEffects40 <- merge(evkOffEffects40, chanLocs, by = "urchan") 
evkOffEffects40 <- evkOffEffects40 %>% mutate(logp = -1*log10(p.value))

evkOff40 <- lunaize(ggplot(evkOffEffects40, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
          ggtitle("Evoked vs Offset, 40 hz")  + theme(text = element_text(size = 30))) 


### 30 hz ----
evkOffEffects30 <- SNR_fooof %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Evoked','Offset'))

evkOffEffects30 <- merge(evkOffEffects30, chanLocs, by = "urchan") 
evkOffEffects30 <- evkOffEffects30 %>% mutate(logp = -1*log10(p.value))

evkOff30 <- lunaize(ggplot(evkOffEffects30, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
                      ggtitle("Evoked vs Offset, 30 hz")  + theme(text = element_text(size = 30))) 

### 20 hz ----
evkOffEffects20 <- SNR_fooof %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Evoked','Offset'))

evkOffEffects20 <- merge(evkOffEffects20, chanLocs, by = "urchan") 
evkOffEffects20 <- evkOffEffects20 %>% mutate(logp = -1*log10(p.value))

evkOff20 <- lunaize(ggplot(evkOffEffects20, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
                      ggtitle("Evoked vs Offset, 20 hz")  + theme(text = element_text(size = 30))) 

cowplot::plot_grid(evkOff20, evkOff30, evkOff40)

## FC4 ----
### gam ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'FC4'), aes(x = Evoked, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Evoked") + ylab("Offset")+ ggtitle('Evoked vs Offset at FC4')

gam.model <-  gamm(Offset ~ Evoked + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'FC4'), random=list(lunaID=~1))
summary(gam.model$gam)

### linear ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'FC4'), aes(x = Evoked, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Evoked") + ylab("Offset")+ ggtitle('Evoked vs Offset at FC4')

summary(lmerTest::lmer(Offset ~ Evoked + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 40 & labels == 'FC4')))


# Induced vs Offset  ----
## Across Electrodes ----

### 40 hz ----
indOffEffects40 <- SNR_fooof %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Induced','Offset'))

indOffEffects40 <- merge(indOffEffects40, chanLocs, by = "urchan") 
indOffEffects40 <- indOffEffects40 %>% mutate(logp = -1*log10(p.value))

indOff40 <- lunaize(ggplot(indOffEffects40, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
          ggtitle("Induced vs Offset, 40 hz")  + theme(text = element_text(size = 30))) 

### 30 hz ----
indOffEffects30 <- SNR_fooof %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Induced','Offset'))

indOffEffects30 <- merge(indOffEffects30, chanLocs, by = "urchan") 
indOffEffects30 <- indOffEffects30 %>% mutate(logp = -1*log10(p.value))

indOff30 <- lunaize(ggplot(indOffEffects30, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
                      ggtitle("Induced vs Offset, 30 hz")  + theme(text = element_text(size = 30))) 

### 20 hz ----
indOffEffects20 <- SNR_fooof %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_fooof(thisdf=., 'Induced','Offset'))

indOffEffects20 <- merge(indOffEffects20, chanLocs, by = "urchan") 
indOffEffects20 <- indOffEffects20 %>% mutate(logp = -1*log10(p.value))

indOff20 <- lunaize(ggplot(indOffEffects20, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
                      geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),NA)) + 
                      ggtitle("Induced vs Offset, 20 hz")  + theme(text = element_text(size = 30))) 

cowplot::plot_grid(indOff20, indOff30, indOff40)


## CPz ----
### gam ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'CPz'), aes(x = Induced, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Induced") + ylab("Offset")+ ggtitle('Induced vs Offset at CPz')

gam.model <-  gamm(Offset ~ Induced + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'CPz'), random=list(lunaID=~1))
summary(gam.model$gam)

### linear ----
#### 40 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & hertz == 40 & labels == 'CPz'), aes(x = Induced, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Induced") + ylab("Offset")+ ggtitle('Induced vs Offset at CPz')

summary(lmerTest::lmer(Offset ~ Induced + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 40 & hertz == 40 & labels == 'CPz')))

#### 30 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 30 & hertz == 30 & labels == 'CPz'), aes(x = Induced, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Induced") + ylab("Offset")+ ggtitle('Induced vs Offset at CPz')

summary(lmerTest::lmer(Offset ~ Induced + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 30 & hertz == 30 & labels == 'CPz')))

#### 20 hz ----
ggplot(data = SNR_fooof %>% filter(freqs == 20 & hertz == 20 & labels == 'CPz'), aes(x = Induced, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Induced") + ylab("Offset")+ ggtitle('Induced vs Offset at CPz')

summary(lmerTest::lmer(Offset ~ Induced + age + Condition + (1|lunaID), data = SNR_fooof %>% filter(freqs == 20 & hertz == 20 & labels == 'CPz')))

#### All Freqs ----
ggplot(data = SNR_fooof %>% filter(freqs == hertz & labels == 'CPz'), aes(x = Induced, y = Offset, color = as.factor(hertz))) +
  stat_smooth(method=lm) + 
  xlab("Induced") + ylab("Offset")+ ggtitle('Induced vs Offset at CPz')

## P2 ----
ggplot(data = SNR_fooof %>% filter(freqs == 40 & labels == 'P2'), aes(x = Induced, y = Offset)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T)) + 
  xlab("Induced") + ylab("Offset")+ ggtitle('Induced vs Offset at P2')

gam.model <-  gamm(Offset ~ Induced + s(age, k = 3) + Condition, data = SNR_fooof %>% 
                     filter(freqs == 40 & labels == 'P2'), random=list(lunaID=~1))
summary(gam.model$gam)



# Load in MRSI from Merge 7T ----

MRSvalues <- merge7t[c("lunaid","sipfc.date","visitno","sipfc.Thalamus_GABA_gamadj", "sipfc.Thalamus_Glu_gamadj", "sipfc.RDLPFC_Glu_gamadj", "sipfc.LDLPFC_Glu_gamadj", 
                       "sipfc.RDLPFC_GABA_gamadj", "sipfc.LDLPFC_GABA_gamadj")]

MRSlong <- MRSvalues %>% 
  select(matches('lunaid|visitno|sipfc.date|(sipfc).*[LR]DLPFC.*(gamadj)')) %>%
  pivot_longer(cols=matches('DLPFC'),
               names_pattern='(.*).([LR]DLPFC)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','sipfc.date','Region'),
              names_from=c('data','measure'))


colnames(MRSlong) <- c("lunaID","visitno","visitDate", "Region", "Glu", "GABA")


### Create MRSI derivatives ----

idx <- which(!is.na(MRSlong$Glu) & !is.na(MRSlong$GABA))
gabaglu.lm <- lm(Glu ~ GABA + Region, data = MRSlong[idx,])
MRSlong$GluGABAimbalance <- NA
MRSlong$GluGABAimbalanceABS <- NA
MRSlong[idx,]$GluGABAimbalanceABS <- abs(gabaglu.lm$residuals)
MRSlong[idx,]$GluGABAimbalance <- (gabaglu.lm$residuals)

MRSlong$Ratio <- MRSlong$Glu/MRSlong$GABA

MRSlong<- MRSlong %>% pivot_wider(names_from = Region,
                                  values_from = c("Glu", "GABA", "Ratio", "GluGABAimbalance", "GluGABAimbalanceABS"),
                                  names_glue = "{Region}_{.value}")


# calculate imbalance with thalamus 

MRSlongThalamus <- MRSvalues %>% 
  select(matches('lunaid|visitno|sipfc.date|(sipfc).*Thalamus.*(gamadj)')) %>%
  pivot_longer(cols=matches('Thalamus.'),
               names_pattern='(.*).(Thalamus)_(.*)',
               names_to=c("data", "Region","measure"))  %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols=c('lunaid','visitno','sipfc.date','Region'),
              names_from=c('data','measure'))


colnames(MRSlongThalamus) <- c("lunaID","visitno","visitDate", "Region", "Glu", "GABA")


### Create MRSI derivatives ----

idx <- which(!is.na(MRSlongThalamus$Glu) & !is.na(MRSlongThalamus$GABA))
gabaglu.lm <- lm(Glu ~ GABA, data = MRSlongThalamus[idx,])
MRSlongThalamus$GluGABAimbalance <- NA
MRSlongThalamus$GluGABAimbalanceABS <- NA
MRSlongThalamus[idx,]$GluGABAimbalanceABS <- abs(gabaglu.lm$residuals)
MRSlongThalamus[idx,]$GluGABAimbalance <- (gabaglu.lm$residuals)

MRSlongThalamus$Ratio <- MRSlongThalamus$Glu/MRSlongThalamus$GABA

MRSlongThalamus<- MRSlongThalamus %>% pivot_wider(names_from = Region,
                                  values_from = c("Glu", "GABA", "Ratio", "GluGABAimbalance", "GluGABAimbalanceABS"),
                                  names_glue = "{Region}_{.value}")

MRSregions <- merge(MRSlongThalamus, MRSlong, by = c("lunaID", "visitno", "visitDate"))

# Merge MRS and SNR ----
SNR_MRS <- merge(SNRallChans_outlier, MRSregions, by = c("lunaID", "visitno"))
write.csv(SNR_MRS, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR_MRS.csv')

write.csv(SNR_MRS, '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Lab/Projects/SNR/rMarkdown/allSubjectsSNR_MRS.csv')


# SNR vs RDLPFC GABA ----
## Across Electrodes ----
### 40-40 hz ----
SNR_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects4040 <- merge(SNR_MRSEffects4040, chanLocs, by = "urchan") 
SNR_MRSEffects4040 <- SNR_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
SNR_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects4030 <- merge(SNR_MRSEffects4030, chanLocs, by = "urchan") 
SNR_MRSEffects4030 <- SNR_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
SNR_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects4020 <- merge(SNR_MRSEffects4020, chanLocs, by = "urchan") 
SNR_MRSEffects4020 <- SNR_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
SNR_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects3040 <- merge(SNR_MRSEffects3040, chanLocs, by = "urchan") 
SNR_MRSEffects3040 <- SNR_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
SNR_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects3030 <- merge(SNR_MRSEffects3030, chanLocs, by = "urchan") 
SNR_MRSEffects3030 <- SNR_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
SNR_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects3020 <- merge(SNR_MRSEffects3020, chanLocs, by = "urchan") 
SNR_MRSEffects3020 <- SNR_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
SNR_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects2040 <- merge(SNR_MRSEffects2040, chanLocs, by = "urchan") 
SNR_MRSEffects2040 <- SNR_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
SNR_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects2030 <- merge(SNR_MRSEffects2030, chanLocs, by = "urchan") 
SNR_MRSEffects2030 <- SNR_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
SNR_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_GABA'))

SNR_MRSEffects2020 <- merge(SNR_MRSEffects2020, chanLocs, by = "urchan") 
SNR_MRSEffects2020 <- SNR_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(SNR_MRSEffects4020, SNR_MRSEffects4030) %>% rbind(., SNR_MRSEffects4040) %>% rbind(., SNR_MRSEffects3040) %>% 
  rbind(., SNR_MRSEffects3030) %>% rbind(., SNR_MRSEffects3020) %>% rbind(., SNR_MRSEffects2040) %>% rbind(., SNR_MRSEffects2030) %>% rbind(., SNR_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs RDLPFC GABA")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNRageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs RDLPFC GABA")


# SNR vs LDLPFC GABA ----
## Across Electrodes ----
### 40-40 hz ----
SNR_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects4040 <- merge(SNR_MRSEffects4040, chanLocs, by = "urchan") 
SNR_MRSEffects4040 <- SNR_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
SNR_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects4030 <- merge(SNR_MRSEffects4030, chanLocs, by = "urchan") 
SNR_MRSEffects4030 <- SNR_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
SNR_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects4020 <- merge(SNR_MRSEffects4020, chanLocs, by = "urchan") 
SNR_MRSEffects4020 <- SNR_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
SNR_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects3040 <- merge(SNR_MRSEffects3040, chanLocs, by = "urchan") 
SNR_MRSEffects3040 <- SNR_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
SNR_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects3030 <- merge(SNR_MRSEffects3030, chanLocs, by = "urchan") 
SNR_MRSEffects3030 <- SNR_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
SNR_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects3020 <- merge(SNR_MRSEffects3020, chanLocs, by = "urchan") 
SNR_MRSEffects3020 <- SNR_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
SNR_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects2040 <- merge(SNR_MRSEffects2040, chanLocs, by = "urchan") 
SNR_MRSEffects2040 <- SNR_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
SNR_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects2030 <- merge(SNR_MRSEffects2030, chanLocs, by = "urchan") 
SNR_MRSEffects2030 <- SNR_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
SNR_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_GABA'))

SNR_MRSEffects2020 <- merge(SNR_MRSEffects2020, chanLocs, by = "urchan") 
SNR_MRSEffects2020 <- SNR_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(SNR_MRSEffects4020, SNR_MRSEffects4030) %>% rbind(., SNR_MRSEffects4040) %>% rbind(., SNR_MRSEffects3040) %>% 
  rbind(., SNR_MRSEffects3030) %>% rbind(., SNR_MRSEffects3020) %>% rbind(., SNR_MRSEffects2040) %>% rbind(., SNR_MRSEffects2030) %>% rbind(., SNR_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs LDLPFC GABA")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs LDLPFC GABA")



# SNR vs RDLPFC Glu ----
## Across Electrodes ----
### 40-40 hz ----
SNR_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects4040 <- merge(SNR_MRSEffects4040, chanLocs, by = "urchan") 
SNR_MRSEffects4040 <- SNR_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
SNR_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects4030 <- merge(SNR_MRSEffects4030, chanLocs, by = "urchan") 
SNR_MRSEffects4030 <- SNR_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
SNR_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects4020 <- merge(SNR_MRSEffects4020, chanLocs, by = "urchan") 
SNR_MRSEffects4020 <- SNR_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
SNR_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects3040 <- merge(SNR_MRSEffects3040, chanLocs, by = "urchan") 
SNR_MRSEffects3040 <- SNR_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
SNR_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects3030 <- merge(SNR_MRSEffects3030, chanLocs, by = "urchan") 
SNR_MRSEffects3030 <- SNR_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
SNR_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects3020 <- merge(SNR_MRSEffects3020, chanLocs, by = "urchan") 
SNR_MRSEffects3020 <- SNR_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
SNR_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects2040 <- merge(SNR_MRSEffects2040, chanLocs, by = "urchan") 
SNR_MRSEffects2040 <- SNR_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
SNR_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects2030 <- merge(SNR_MRSEffects2030, chanLocs, by = "urchan") 
SNR_MRSEffects2030 <- SNR_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
SNR_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','RDLPFC_Glu'))

SNR_MRSEffects2020 <- merge(SNR_MRSEffects2020, chanLocs, by = "urchan") 
SNR_MRSEffects2020 <- SNR_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(SNR_MRSEffects4020, SNR_MRSEffects4030) %>% rbind(., SNR_MRSEffects4040) %>% rbind(., SNR_MRSEffects3040) %>% 
  rbind(., SNR_MRSEffects3030) %>% rbind(., SNR_MRSEffects3020) %>% rbind(., SNR_MRSEffects2040) %>% rbind(., SNR_MRSEffects2030) %>% rbind(., SNR_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs RDLPFC Glu")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs RDLPFC Glu")


# SNR vs LDLPFC Glu ----
## Across Electrodes ----
### 40-40 hz ----
SNR_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects4040 <- merge(SNR_MRSEffects4040, chanLocs, by = "urchan") 
SNR_MRSEffects4040 <- SNR_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
SNR_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects4030 <- merge(SNR_MRSEffects4030, chanLocs, by = "urchan") 
SNR_MRSEffects4030 <- SNR_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
SNR_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects4020 <- merge(SNR_MRSEffects4020, chanLocs, by = "urchan") 
SNR_MRSEffects4020 <- SNR_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
SNR_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects3040 <- merge(SNR_MRSEffects3040, chanLocs, by = "urchan") 
SNR_MRSEffects3040 <- SNR_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
SNR_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects3030 <- merge(SNR_MRSEffects3030, chanLocs, by = "urchan") 
SNR_MRSEffects3030 <- SNR_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
SNR_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects3020 <- merge(SNR_MRSEffects3020, chanLocs, by = "urchan") 
SNR_MRSEffects3020 <- SNR_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
SNR_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects2040 <- merge(SNR_MRSEffects2040, chanLocs, by = "urchan") 
SNR_MRSEffects2040 <- SNR_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
SNR_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects2030 <- merge(SNR_MRSEffects2030, chanLocs, by = "urchan") 
SNR_MRSEffects2030 <- SNR_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
SNR_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'SNR','LDLPFC_Glu'))

SNR_MRSEffects2020 <- merge(SNR_MRSEffects2020, chanLocs, by = "urchan") 
SNR_MRSEffects2020 <- SNR_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
SNR_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(SNR_MRSEffects4020, SNR_MRSEffects4030) %>% rbind(., SNR_MRSEffects4040) %>% rbind(., SNR_MRSEffects3040) %>% 
  rbind(., SNR_MRSEffects3030) %>% rbind(., SNR_MRSEffects3020) %>% rbind(., SNR_MRSEffects2040) %>% rbind(., SNR_MRSEffects2030) %>% rbind(., SNR_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs LDLPFC Glu")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("SNR vs LDLPFC Glu")


# Evoked vs RDLPFC Glu ----
## Across Electrodes ----
### 40-40 hz ----
Evoked_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects4040 <- merge(Evoked_MRSEffects4040, chanLocs, by = "urchan") 
Evoked_MRSEffects4040 <- Evoked_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Evoked_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects4030 <- merge(Evoked_MRSEffects4030, chanLocs, by = "urchan") 
Evoked_MRSEffects4030 <- Evoked_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Evoked_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects4020 <- merge(Evoked_MRSEffects4020, chanLocs, by = "urchan") 
Evoked_MRSEffects4020 <- Evoked_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Evoked_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects3040 <- merge(Evoked_MRSEffects3040, chanLocs, by = "urchan") 
Evoked_MRSEffects3040 <- Evoked_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Evoked_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects3030 <- merge(Evoked_MRSEffects3030, chanLocs, by = "urchan") 
Evoked_MRSEffects3030 <- Evoked_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Evoked_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects3020 <- merge(Evoked_MRSEffects3020, chanLocs, by = "urchan") 
Evoked_MRSEffects3020 <- Evoked_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Evoked_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects2040 <- merge(Evoked_MRSEffects2040, chanLocs, by = "urchan") 
Evoked_MRSEffects2040 <- Evoked_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Evoked_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects2030 <- merge(Evoked_MRSEffects2030, chanLocs, by = "urchan") 
Evoked_MRSEffects2030 <- Evoked_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Evoked_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_Glu'))

Evoked_MRSEffects2020 <- merge(Evoked_MRSEffects2020, chanLocs, by = "urchan") 
Evoked_MRSEffects2020 <- Evoked_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Evoked_MRSEffects4020, Evoked_MRSEffects4030) %>% rbind(., Evoked_MRSEffects4040) %>% rbind(., Evoked_MRSEffects3040) %>% 
  rbind(., Evoked_MRSEffects3030) %>% rbind(., Evoked_MRSEffects3020) %>% rbind(., Evoked_MRSEffects2040) %>% rbind(., Evoked_MRSEffects2030) %>% rbind(., Evoked_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs RDLPFC Glu")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs RDLPFC Glu")


# Evoked vs LDLPFC Glu ----
## Across Electrodes ----
### 40-40 hz ----
Evoked_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects4040 <- merge(Evoked_MRSEffects4040, chanLocs, by = "urchan") 
Evoked_MRSEffects4040 <- Evoked_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Evoked_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects4030 <- merge(Evoked_MRSEffects4030, chanLocs, by = "urchan") 
Evoked_MRSEffects4030 <- Evoked_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Evoked_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects4020 <- merge(Evoked_MRSEffects4020, chanLocs, by = "urchan") 
Evoked_MRSEffects4020 <- Evoked_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Evoked_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects3040 <- merge(Evoked_MRSEffects3040, chanLocs, by = "urchan") 
Evoked_MRSEffects3040 <- Evoked_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Evoked_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects3030 <- merge(Evoked_MRSEffects3030, chanLocs, by = "urchan") 
Evoked_MRSEffects3030 <- Evoked_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Evoked_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects3020 <- merge(Evoked_MRSEffects3020, chanLocs, by = "urchan") 
Evoked_MRSEffects3020 <- Evoked_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Evoked_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects2040 <- merge(Evoked_MRSEffects2040, chanLocs, by = "urchan") 
Evoked_MRSEffects2040 <- Evoked_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Evoked_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects2030 <- merge(Evoked_MRSEffects2030, chanLocs, by = "urchan") 
Evoked_MRSEffects2030 <- Evoked_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Evoked_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_Glu'))

Evoked_MRSEffects2020 <- merge(Evoked_MRSEffects2020, chanLocs, by = "urchan") 
Evoked_MRSEffects2020 <- Evoked_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Evoked_MRSEffects4020, Evoked_MRSEffects4030) %>% rbind(., Evoked_MRSEffects4040) %>% rbind(., Evoked_MRSEffects3040) %>% 
  rbind(., Evoked_MRSEffects3030) %>% rbind(., Evoked_MRSEffects3020) %>% rbind(., Evoked_MRSEffects2040) %>% rbind(., Evoked_MRSEffects2030) %>% rbind(., Evoked_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs LDLPFC Glu")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs LDLPFC Glu")


# Evoked vs RDLPFC GABA ----
## Across Electrodes ----
### 40-40 hz ----
Evoked_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects4040 <- merge(Evoked_MRSEffects4040, chanLocs, by = "urchan") 
Evoked_MRSEffects4040 <- Evoked_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Evoked_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects4030 <- merge(Evoked_MRSEffects4030, chanLocs, by = "urchan") 
Evoked_MRSEffects4030 <- Evoked_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Evoked_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects4020 <- merge(Evoked_MRSEffects4020, chanLocs, by = "urchan") 
Evoked_MRSEffects4020 <- Evoked_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Evoked_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects3040 <- merge(Evoked_MRSEffects3040, chanLocs, by = "urchan") 
Evoked_MRSEffects3040 <- Evoked_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Evoked_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects3030 <- merge(Evoked_MRSEffects3030, chanLocs, by = "urchan") 
Evoked_MRSEffects3030 <- Evoked_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Evoked_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects3020 <- merge(Evoked_MRSEffects3020, chanLocs, by = "urchan") 
Evoked_MRSEffects3020 <- Evoked_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Evoked_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects2040 <- merge(Evoked_MRSEffects2040, chanLocs, by = "urchan") 
Evoked_MRSEffects2040 <- Evoked_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Evoked_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects2030 <- merge(Evoked_MRSEffects2030, chanLocs, by = "urchan") 
Evoked_MRSEffects2030 <- Evoked_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Evoked_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','RDLPFC_GABA'))

Evoked_MRSEffects2020 <- merge(Evoked_MRSEffects2020, chanLocs, by = "urchan") 
Evoked_MRSEffects2020 <- Evoked_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Evoked_MRSEffects4020, Evoked_MRSEffects4030) %>% rbind(., Evoked_MRSEffects4040) %>% rbind(., Evoked_MRSEffects3040) %>% 
  rbind(., Evoked_MRSEffects3030) %>% rbind(., Evoked_MRSEffects3020) %>% rbind(., Evoked_MRSEffects2040) %>% rbind(., Evoked_MRSEffects2030) %>% rbind(., Evoked_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs RDLPFC GABA")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs RDLPFC GABA")


# Evoked vs LDLPFC GABA ----
## Across Electrodes ----
### 40-40 hz ----
Evoked_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects4040 <- merge(Evoked_MRSEffects4040, chanLocs, by = "urchan") 
Evoked_MRSEffects4040 <- Evoked_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Evoked_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects4030 <- merge(Evoked_MRSEffects4030, chanLocs, by = "urchan") 
Evoked_MRSEffects4030 <- Evoked_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Evoked_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects4020 <- merge(Evoked_MRSEffects4020, chanLocs, by = "urchan") 
Evoked_MRSEffects4020 <- Evoked_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Evoked_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects3040 <- merge(Evoked_MRSEffects3040, chanLocs, by = "urchan") 
Evoked_MRSEffects3040 <- Evoked_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Evoked_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects3030 <- merge(Evoked_MRSEffects3030, chanLocs, by = "urchan") 
Evoked_MRSEffects3030 <- Evoked_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Evoked_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects3020 <- merge(Evoked_MRSEffects3020, chanLocs, by = "urchan") 
Evoked_MRSEffects3020 <- Evoked_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Evoked_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects2040 <- merge(Evoked_MRSEffects2040, chanLocs, by = "urchan") 
Evoked_MRSEffects2040 <- Evoked_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Evoked_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects2030 <- merge(Evoked_MRSEffects2030, chanLocs, by = "urchan") 
Evoked_MRSEffects2030 <- Evoked_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Evoked_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Evoked','LDLPFC_GABA'))

Evoked_MRSEffects2020 <- merge(Evoked_MRSEffects2020, chanLocs, by = "urchan") 
Evoked_MRSEffects2020 <- Evoked_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Evoked_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Evoked_MRSEffects4020, Evoked_MRSEffects4030) %>% rbind(., Evoked_MRSEffects4040) %>% rbind(., Evoked_MRSEffects3040) %>% 
  rbind(., Evoked_MRSEffects3030) %>% rbind(., Evoked_MRSEffects3020) %>% rbind(., Evoked_MRSEffects2040) %>% rbind(., Evoked_MRSEffects2030) %>% rbind(., Evoked_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs LDLPFC GABA")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Evoked vs LDLPFC GABA")



# Induced vs RDLPFC Glu ----
## Across Electrodes ----
### 40-40 hz ----
Induced_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects4040 <- merge(Induced_MRSEffects4040, chanLocs, by = "urchan") 
Induced_MRSEffects4040 <- Induced_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Induced_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects4030 <- merge(Induced_MRSEffects4030, chanLocs, by = "urchan") 
Induced_MRSEffects4030 <- Induced_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Induced_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects4020 <- merge(Induced_MRSEffects4020, chanLocs, by = "urchan") 
Induced_MRSEffects4020 <- Induced_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Induced_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects3040 <- merge(Induced_MRSEffects3040, chanLocs, by = "urchan") 
Induced_MRSEffects3040 <- Induced_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Induced_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects3030 <- merge(Induced_MRSEffects3030, chanLocs, by = "urchan") 
Induced_MRSEffects3030 <- Induced_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Induced_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects3020 <- merge(Induced_MRSEffects3020, chanLocs, by = "urchan") 
Induced_MRSEffects3020 <- Induced_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Induced_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects2040 <- merge(Induced_MRSEffects2040, chanLocs, by = "urchan") 
Induced_MRSEffects2040 <- Induced_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Induced_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects2030 <- merge(Induced_MRSEffects2030, chanLocs, by = "urchan") 
Induced_MRSEffects2030 <- Induced_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Induced_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_Glu'))

Induced_MRSEffects2020 <- merge(Induced_MRSEffects2020, chanLocs, by = "urchan") 
Induced_MRSEffects2020 <- Induced_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Induced_MRSEffects4020, Induced_MRSEffects4030) %>% rbind(., Induced_MRSEffects4040) %>% rbind(., Induced_MRSEffects3040) %>% 
  rbind(., Induced_MRSEffects3030) %>% rbind(., Induced_MRSEffects3020) %>% rbind(., Induced_MRSEffects2040) %>% rbind(., Induced_MRSEffects2030) %>% rbind(., Induced_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs RDLPFC Glu")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs RDLPFC Glu")


# Induced vs LDLPFC Glu ----
## Across Electrodes ----
### 40-40 hz ----
Induced_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects4040 <- merge(Induced_MRSEffects4040, chanLocs, by = "urchan") 
Induced_MRSEffects4040 <- Induced_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Induced_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects4030 <- merge(Induced_MRSEffects4030, chanLocs, by = "urchan") 
Induced_MRSEffects4030 <- Induced_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Induced_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects4020 <- merge(Induced_MRSEffects4020, chanLocs, by = "urchan") 
Induced_MRSEffects4020 <- Induced_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Induced_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects3040 <- merge(Induced_MRSEffects3040, chanLocs, by = "urchan") 
Induced_MRSEffects3040 <- Induced_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Induced_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects3030 <- merge(Induced_MRSEffects3030, chanLocs, by = "urchan") 
Induced_MRSEffects3030 <- Induced_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Induced_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects3020 <- merge(Induced_MRSEffects3020, chanLocs, by = "urchan") 
Induced_MRSEffects3020 <- Induced_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Induced_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects2040 <- merge(Induced_MRSEffects2040, chanLocs, by = "urchan") 
Induced_MRSEffects2040 <- Induced_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Induced_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects2030 <- merge(Induced_MRSEffects2030, chanLocs, by = "urchan") 
Induced_MRSEffects2030 <- Induced_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Induced_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_Glu'))

Induced_MRSEffects2020 <- merge(Induced_MRSEffects2020, chanLocs, by = "urchan") 
Induced_MRSEffects2020 <- Induced_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Induced_MRSEffects4020, Induced_MRSEffects4030) %>% rbind(., Induced_MRSEffects4040) %>% rbind(., Induced_MRSEffects3040) %>% 
  rbind(., Induced_MRSEffects3030) %>% rbind(., Induced_MRSEffects3020) %>% rbind(., Induced_MRSEffects2040) %>% rbind(., Induced_MRSEffects2030) %>% rbind(., Induced_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs LDLPFC Glu")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs LDLPFC Glu")


# Induced vs RDLPFC GABA ----
## Across Electrodes ----
### 40-40 hz ----
Induced_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects4040 <- merge(Induced_MRSEffects4040, chanLocs, by = "urchan") 
Induced_MRSEffects4040 <- Induced_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Induced_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects4030 <- merge(Induced_MRSEffects4030, chanLocs, by = "urchan") 
Induced_MRSEffects4030 <- Induced_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Induced_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects4020 <- merge(Induced_MRSEffects4020, chanLocs, by = "urchan") 
Induced_MRSEffects4020 <- Induced_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Induced_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects3040 <- merge(Induced_MRSEffects3040, chanLocs, by = "urchan") 
Induced_MRSEffects3040 <- Induced_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Induced_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects3030 <- merge(Induced_MRSEffects3030, chanLocs, by = "urchan") 
Induced_MRSEffects3030 <- Induced_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Induced_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects3020 <- merge(Induced_MRSEffects3020, chanLocs, by = "urchan") 
Induced_MRSEffects3020 <- Induced_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Induced_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects2040 <- merge(Induced_MRSEffects2040, chanLocs, by = "urchan") 
Induced_MRSEffects2040 <- Induced_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Induced_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects2030 <- merge(Induced_MRSEffects2030, chanLocs, by = "urchan") 
Induced_MRSEffects2030 <- Induced_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Induced_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','RDLPFC_GABA'))

Induced_MRSEffects2020 <- merge(Induced_MRSEffects2020, chanLocs, by = "urchan") 
Induced_MRSEffects2020 <- Induced_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Induced_MRSEffects4020, Induced_MRSEffects4030) %>% rbind(., Induced_MRSEffects4040) %>% rbind(., Induced_MRSEffects3040) %>% 
  rbind(., Induced_MRSEffects3030) %>% rbind(., Induced_MRSEffects3020) %>% rbind(., Induced_MRSEffects2040) %>% rbind(., Induced_MRSEffects2030) %>% rbind(., Induced_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs RDLPFC GABA")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs RDLPFC GABA")


# Induced vs LDLPFC GABA ----
## Across Electrodes ----
### 40-40 hz ----
Induced_MRSEffects4040 <- SNR_MRS %>% filter(freqs == 40 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects4040 <- merge(Induced_MRSEffects4040, chanLocs, by = "urchan") 
Induced_MRSEffects4040 <- Induced_MRSEffects4040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4040$hertz <- '40 hz'

### 40-30 hz ----
Induced_MRSEffects4030 <- SNR_MRS %>% filter(freqs == 40 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects4030 <- merge(Induced_MRSEffects4030, chanLocs, by = "urchan") 
Induced_MRSEffects4030 <- Induced_MRSEffects4030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4030$hertz <- '30 hz'

### 40-20 hz ----
Induced_MRSEffects4020 <- SNR_MRS %>% filter(freqs == 40 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects4020 <- merge(Induced_MRSEffects4020, chanLocs, by = "urchan") 
Induced_MRSEffects4020 <- Induced_MRSEffects4020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects4020$hertz <- '20 hz'

### 30-40 hz ----
Induced_MRSEffects3040 <- SNR_MRS %>% filter(freqs == 30 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects3040 <- merge(Induced_MRSEffects3040, chanLocs, by = "urchan") 
Induced_MRSEffects3040 <- Induced_MRSEffects3040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3040$hertz <- '40 hz'

### 30-30 hz ----
Induced_MRSEffects3030 <- SNR_MRS %>% filter(freqs == 30 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects3030 <- merge(Induced_MRSEffects3030, chanLocs, by = "urchan") 
Induced_MRSEffects3030 <- Induced_MRSEffects3030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3030$hertz <- '30 hz'

### 30-20 hz ----
Induced_MRSEffects3020 <- SNR_MRS %>% filter(freqs == 30 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects3020 <- merge(Induced_MRSEffects3020, chanLocs, by = "urchan") 
Induced_MRSEffects3020 <- Induced_MRSEffects3020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects3020$hertz <- '20 hz'

### 20-40 hz ----
Induced_MRSEffects2040 <- SNR_MRS %>% filter(freqs == 20 & hertz == 40) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects2040 <- merge(Induced_MRSEffects2040, chanLocs, by = "urchan") 
Induced_MRSEffects2040 <- Induced_MRSEffects2040 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2040$hertz <- '40 hz'

### 20-30 hz ----
Induced_MRSEffects2030 <- SNR_MRS %>% filter(freqs == 20 & hertz == 30) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects2030 <- merge(Induced_MRSEffects2030, chanLocs, by = "urchan") 
Induced_MRSEffects2030 <- Induced_MRSEffects2030 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2030$hertz <- '30 hz'

### 20-20 hz ----
Induced_MRSEffects2020 <- SNR_MRS %>% filter(freqs == 20 & hertz == 20) %>%
  mutate(invage = 1/age) %>% 
  group_by(urchan, freqs) %>% 
  dplyr::do(fit_gam_mrs(thisdf=., 'Induced','LDLPFC_GABA'))

Induced_MRSEffects2020 <- merge(Induced_MRSEffects2020, chanLocs, by = "urchan") 
Induced_MRSEffects2020 <- Induced_MRSEffects2020 %>% mutate(logp = -1*log10(p.value))
Induced_MRSEffects2020$hertz <- '20 hz'

### Combine all combinations 

SNR_MRSageEffectsAll <- rbind(Induced_MRSEffects4020, Induced_MRSEffects4030) %>% rbind(., Induced_MRSEffects4040) %>% rbind(., Induced_MRSEffects3040) %>% 
  rbind(., Induced_MRSEffects3030) %>% rbind(., Induced_MRSEffects3020) %>% rbind(., Induced_MRSEffects2040) %>% rbind(., Induced_MRSEffects2030) %>% rbind(., Induced_MRSEffects2020)

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = Fvalue, z = Fvalue, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red") + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs LDLPFC GABA")

lunaize(ggplot(SNR_MRSageEffectsAll, aes(x = -Y, y = X, fill = logp, z = logp, label = labels)) + 
          geom_topo(chan_markers = "text") + scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1*log10(.05),max(SNR_MRSageEffectsAll$logp)+.1)) + 
          theme(text = element_text(size = 30))) + facet_wrap(~freqs+hertz) + ggtitle("Induced vs LDLPFC GABA")





## Per Electrode ----
### Cz ----
ggplot(data = SNR_MRS %>% filter(freqs == 30 & hertz == 30 & labels == 'Cz'), aes(x = Evoked, y = Thalamus_Glu)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Evoked") + ylab("Glu")+ ggtitle('Evoked vs Thalamus Glu at Cz, 30 Hz')

summary(lmerTest::lmer(Thalamus_Glu ~ Evoked + age  + (1|lunaID), data = SNR_MRS %>% filter(freqs == 30 & hertz == 30 & labels == 'Cz')))


### C2 ----
ggplot(data = SNR_MRS %>% filter(freqs == 30 & hertz == 30 & labels == 'C2'), aes(x = Evoked, y = Thalamus_Glu)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Evoked") + ylab("Glu")+ ggtitle('Evoked vs Thalamus Glu at C2, 30 Hz')

summary(lmerTest::lmer(Thalamus_Glu ~ Evoked + age  + (1|lunaID), data = SNR_MRS %>% filter(freqs == 30 & hertz == 30 & labels == 'C2')))


### FC3 ----
ggplot(data = SNR_MRS %>% filter(freqs == 20 & hertz == 20 & labels == 'FC3'), aes(x = Evoked, y = Thalamus_Glu)) +
  geom_line(aes(group=interaction(lunaID)), alpha = 0.2) + 
  geom_point() + stat_smooth(method=lm) + 
  xlab("Evoked") + ylab("Glu")+ ggtitle('Evoked vs Thalamus Glu at FC3, 20 Hz')

summary(lmerTest::lmer(Thalamus_Glu ~ Evoked + age  + (1|lunaID), data = SNR_MRS %>% filter(freqs == 20 & hertz == 20 & labels == 'FC3')))


