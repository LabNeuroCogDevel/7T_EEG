
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

# Load in dataframe ----

itc <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubsTFvalues.csv')

str(itc)

itc_byGroup <- itc %>% 
  mutate(ageGroup = as.factor(cut(age, c(0,17,23,Inf), labels = c('10-16','17-22','23-30')))) %>%
  group_by(ageGroup, time, freqs) %>%
  summarize(meanITC = mean(ITC, na.rm=T),
            sdITC = sd(ITC, na.rm=T),
            n = n()) %>%
  mutate(seITC = sdITC / sqrt(n))

unique(itc$freqs)
unique(itc$time)

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

ggplot(data = itc %>% filter(abs(freqs-40.6)<.2 & abs(time-252.15)<.1), aes(x = age, y = ITC)) + 
  geom_point() + stat_smooth(method='loess') +
  labs(title = '40 Hz, t=250ms')

ggplot(data = itc %>% filter(abs(freqs-40.6)<.2 & time>200 & time < 300) %>% group_by(lunaID, visitDate, age) %>% summarize(ITC = mean(ITC, na.rm=T)), 
       aes(x = age, y = ITC)) + 
  geom_point() + stat_smooth(method='loess') +
  labs(title = '40 Hz, t=200-300ms')


cowplot::plot_grid(
  ggplot(data = itc %>% filter(freqs < 15 & time>50 & time < 150) %>% group_by(lunaID, visitDate, age) %>% summarize(ITC = mean(ITC, na.rm=T)), 
         aes(x = age, y = ITC)) + 
    geom_point() + stat_smooth(method='gam', formula=y~s(x, k=7)) +
    labs(title = '10-15 Hz, t=50-150ms'),
  
  ggplot(data = itc %>% filter(freqs > 35 & freqs < 45 & time>200 & time < 300) %>% group_by(lunaID, visitDate, age) %>% summarize(ITC = mean(ITC, na.rm=T)), 
         aes(x = age, y = ITC)) + 
    geom_point() + stat_smooth(method='gam', formula=y~s(x, k=7)) +
    labs(title = '35-45 Hz, t=200-300ms')
)
