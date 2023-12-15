
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

# Load in dataframes ----

itc40 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubsITC_40Hz.csv')
itc40$hertz <- '40'
itc30 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubsITC_30Hz.csv')
itc30$hertz <- '30'
itc20 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubsITC_20Hz.csv')
itc20$hertz <- '20'

itc_allfreq <- rbind(itc40, itc30) %>% rbind(., itc20)

# Find rows where itc is 0
rows_to_delete <- which(itc_allfreq$ITC == 0)

# Delete rows with itc equal to 0
itc_allfreq <- itc_allfreq[-rows_to_delete, ]

str(itc_allfreq)

itc_byGroup <- itc_allfreq %>% 
  mutate(ageGroup = as.factor(cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))) %>%
  group_by(ageGroup, time, freqs, hertz) %>%
  summarize(meanITC = mean(ITC, na.rm=T),
            sdITC = sd(ITC, na.rm=T),
            n = n()) %>%
  mutate(seITC = sdITC / sqrt(n))

unique(itc_allfreq$freqs)
unique(itc_allfreq$time)

cowplot::plot_grid(
  ggplot(data = itc_byGroup %>% filter(abs(freqs-20.2)<.2) %>% filter(hertz =='20') %>% arrange(ageGroup, time), aes(x = time, y = meanITC, color = ageGroup, fill = ageGroup)) +
    geom_ribbon(aes(ymin = meanITC - seITC, ymax = meanITC + seITC), alpha=0.2, linewidth=NA) +
    geom_line() +
    labs(x = 'Time (sec)', y = 'ITC', title = '20 Hz'),
  ggplot(data = itc_byGroup %>% filter(abs(freqs-30.4)<.2) %>% filter(hertz =='30') %>% arrange(ageGroup, time), aes(x = time, y = meanITC, color = ageGroup, fill = ageGroup)) +
    geom_ribbon(aes(ymin = meanITC - seITC, ymax = meanITC + seITC), alpha=0.2, linewidth=NA) +
    geom_line() +
    labs(x = 'Time (sec)', y = 'ITC', title = '30 Hz'),
  ggplot(data = itc_byGroup %>% filter(abs(freqs-40.6)<.2)%>% filter(hertz =='40') %>% arrange(ageGroup, time), aes(x = time, y = meanITC, color = ageGroup, fill = ageGroup)) +
    geom_ribbon(aes(ymin = meanITC - seITC, ymax = meanITC + seITC), alpha=0.2, linewidth=NA) +
    geom_line() +
    labs(x = 'Time (sec)', y = 'ITC', title = '40 Hz')
  
)

ggplot(data = itc_byGroup %>% mutate(f = round(freqs, 1)) %>% filter(f %in% c(20.2, 30.4, 40.6)) %>% arrange(ageGroup, time), aes(x = time, y = meanITC, color = ageGroup, fill = ageGroup)) +
  geom_ribbon(aes(ymin = meanITC - seITC, ymax = meanITC + seITC), alpha=0.2, linewidth=NA) +
  geom_line() + facet_wrap(~hertz+f)

ggplot(data = itc_allfreq %>% mutate(f = round(freqs, 1)) %>% filter(f %in% c(20.2, 30.4, 40.6) & time>200 & time < 500) %>% group_by(lunaID, visitDate, age, hertz,f) %>% summarize(ITC = mean(ITC, na.rm=T)), 
       aes(x = age, y = ITC)) + 
  geom_point() + stat_smooth(method='loess')+ facet_wrap(~hertz+f, scales='free')


cowplot::plot_grid(
  ggplot(data = itc_byGroup %>% filter(abs(freqs-10)<.2) %>% arrange(ageGroup, time), aes(x = time, y = meanITC, color = ageGroup, fill = ageGroup)) +
    geom_ribbon(aes(ymin = meanITC - seITC, ymax = meanITC + seITC), alpha=0.2, linewidth=NA) +
    geom_line() +
    labs(x = 'Time (sec)', y = 'ITC', title = '10 Hz'),
  ggplot(data = itc_byGroup %>% filter(abs(freqs-40.6)<.2) %>% arrange(ageGroup, time), aes(x = time, y = meanITC, color = ageGroup, fill = ageGroup)) +
    geom_ribbon(aes(ymin = meanITC - seITC, ymax = meanITC + seITC), alpha=0.2, linewidth=NA) +
    geom_line() +
    labs(x = 'Time (sec)', y = 'ITC', title = '40 Hz')
  
)

ggplot(data = itc_allfreq %>% filter(abs(freqs-40.6)<.2 & abs(time-252.15)<.1 & hertz =='40'), aes(x = age, y = ITC)) + 
  geom_point() + stat_smooth(method='loess') +
  labs(title = '40 Hz, t=250ms') 

ggplot(data = itc_allfreq %>% filter(abs(freqs-40.6)<.2 & time>200 & time < 300) %>% group_by(lunaID, visitDate, age) %>% summarize(ITC = mean(ITC, na.rm=T)), 
       aes(x = age, y = ITC)) + 
  geom_point() + stat_smooth(method='loess') +
  labs(title = '40 Hz, t=200-300ms')


cowplot::plot_grid(
  ggplot(data = itc_allfreq %>% filter(freqs < 15 & time>50 & time < 150) %>% group_by(lunaID, visitDate, age) %>% summarize(ITC = mean(ITC, na.rm=T)), 
         aes(x = age, y = ITC)) + 
    geom_point() + stat_smooth(method='gam', formula=y~s(x, k=7)) +
    labs(title = '10-15 Hz, t=50-150ms'),
  
  ggplot(data = itc_allfreq %>% filter(freqs > 35 & freqs < 45 & time>200 & time < 300) %>% group_by(lunaID, visitDate, age) %>% summarize(ITC = mean(ITC, na.rm=T)), 
         aes(x = age, y = ITC)) + 
    geom_point() + stat_smooth(method='gam', formula=y~s(x, k=7)) +
    labs(title = '35-45 Hz, t=200-300ms')
)
