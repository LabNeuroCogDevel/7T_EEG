
# Trying to find periods of significant growth 

library(LNCDR)
library(mgcv)
library(plotly)
library(graphics)

SigGrowth <- function() {

alldf <- Only_Take_One_Delay_Bin_TrialLevel()
alldata_TrialLevel <- alldf$alldata_TrialLevel

Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')

#behavioral data contains each trial 
Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
Behavior$absPositionError <- abs(Behavior$PositionError)

alldata_TrialLevel_behavior <- merge(alldata_TrialLevel, Behavior, by = c("Subject", "Trial"))

# First look at PE and gamma events 
# clean the data, removing outliers above 2 SDs
cleanData_PE_events <- alldata_TrialLevel_behavior %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number)) < (sd(Gamma.Event_Number) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))

# Gam model, accounting for age, length of delay period, and the fact that not all data points are unique, ie grouped by Subject 
# these models are only to really look at one variable against age 
gam.model <- gam(absPositionError ~ s(age) + Delay + s(Trial) + s(Subject, bs="re"), data = cleanData_PE_events)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
PEplot <- gam_growthrate_plot(cleanData_PE_events,gam.model,gam.growthrate, agevar = 'age',yvar = 'absPositionError', idvar = 'Subject', draw_points = F)


# Gam model, accounting for age, length of delay period, and the fact that not all data points are unique, ie grouped by Subject 
# these models are only to really look at one variable against age 
gam.model <- gam(Gamma.Event_Number ~ s(age) + Delay + s(Trial) + s(Subject, bs="re"), data = cleanData_PE_events)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
eventsPlot <- gam_growthrate_plot(cleanData_PE_events,gam.model,gam.growthrate, agevar = 'age',yvar = 'Gamma.Event_Number', idvar = 'Subject', draw_points = F)

par(mfrow = c(2,1))
plot(PEplot$both)
plot(eventsPlot$both)


# compare to loess, gam, and lm using inverse age formula  
lunaize(ggplot(data = cleanData_PE_events, aes(x = age, y = Gamma.Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)'))


# clean data for PE and power
cleanData_PE_power <- alldata_TrialLevel_behavior %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power)) < (sd(Gamma.log1p_Trial_Power) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))

# gam model for power
gam.model <- gam(Gamma.log1p_Trial_Power ~ s(age) + Delay + s(Trial) + s(Subject, bs="re"), data = cleanData_PE_power)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
powerPlot <- gam_growthrate_plot(cleanData_PE_events,gam.model,gam.growthrate, agevar = 'age',yvar = 'Gamma.log1p_Trial_Power', idvar = 'Subject', draw_points = F)


# clean data for latency
cleanData_latency <- alldata_TrialLevel_behavior %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))

# gam model for power
gam.model <- gam(mgsLatency ~ s(age) + Delay + s(Trial) + s(Subject, bs="re"), data = cleanData_latency)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
latencyPlot <- gam_growthrate_plot(cleanData_latency, gam.model, gam.growthrate, agevar = 'age',yvar = 'mgsLatency', idvar = 'Subject', draw_points = F)

}

sigGrowth_GammaBetaRatio <- function () {
  
  alldf <- Only_Take_One_Delay_Bin_TrialLevel()
  alldata_TrialLevel <- alldf$alldata_TrialLevel
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  alldata_TrialLevel_behavior <- merge(alldata_TrialLevel, Behavior, by = c("Subject", "Trial"))
  
  alldata_TrialLevel_Behavior$RatioNumberofEvents <- alldata_TrialLevel_Behavior$Gamma.Event_Number/alldata_TrialLevel_Behavior$Beta.Event_Number
  alldata_TrialLevel_Behavior$RatioNumberofEvents[which(is.infinite(alldata_TrialLevel_Behavior$RatioNumberofEvents))] <- 0
  
  # First look at PE and gamma: beta events ratio
  # clean the data, removing outliers above 2 SDs
  cleandata <-  alldata_TrialLevel_Behavior %>% filter(abs(RatioNumberofEvents - mean(RatioNumberofEvents, na.rm= T)) < (sd(RatioNumberofEvents, na.rm= T) * 2)) 
 
  # Gam model, accounting for age, length of delay period, and the fact that not all data points are unique, ie grouped by Subject 
  # these models are only to really look at one variable against age 
  gam.model <- gam(RatioNumberofEvents ~ s(age) + Delay + s(Trial) + s(Subject, bs="re"), data = cleandata)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
  Ratioplot <- gam_growthrate_plot(cleandata, gam.model, gam.growthrate, agevar = 'age',yvar = 'RatioNumberofEvents', idvar = 'Subject', draw_points = F)
  
  
  # gamma:beta ratio variability
  sd_cleanData <- aggregate(. ~ Subject, cleandata, sd)
  sd_cleanData <- merge(sd_cleanData, agefile, by = 'Subject')
  sd_cleanData$age <- sd_cleanData$age.y
  
  gam.model <- gam(RatioNumberofEvents ~ s(age), data = sd_cleanData)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  Ratioplot <- gam_growthrate_plot(sd_cleanData, gam.model, gam.growthrate, agevar = 'age',yvar = 'RatioNumberofEvents', draw_points = F)
  
}

sigGrowth_SubLevel <- function() {
  
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R")
  
  alldf <- DelayOnly_Sublevel()
  alldata_delayOnly <- alldf$alldata_delayOnly
  
  avgBehavior <- Behavior_Sublevel_Maria()
  
  allData <- merge(alldata_delayOnly, avgBehavior, by = c("Subject", "age"))
  
  # Gam model, accounting for age, looking a PE
  # these models are only to really look at one variable against age 
  gam.model <- gam(absPositionError ~ s(age), data = allData)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  gam_growthrate_plot(allData, gam.model, gam.growthrate, agevar = 'age', yvar = 'absPositionError', draw_points = T)
  
  
  # Gam model, accounting for age, looking at gamma events
  # these models are only to really look at one variable against age 
  gam.model <- gam(Gamma.Gamma_Event_Number ~ s(age), data = allData)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  gam_growthrate_plot(allData, gam.model, gam.growthrate, agevar = 'age', yvar = 'Gamma.Gamma_Event_Number', draw_points = F)
  
  
  # Gam model, accounting for age, looking at PE variability 
  # these models are only to really look at one variable against age 
  gam.model <- gam(absPositionError ~ s(age), data = sd_cleanData_PE_events)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  gam_growthrate_plot(sd_cleanData_PE_events, gam.model, gam.growthrate, agevar = 'age', yvar = 'absPositionError', draw_points = F)
  
  
  # Gam model, accounting for age, looking at Gamma Event variability 
  # these models are only to really look at one variable against age 
  gam.model <- gam(Gamma.Gamma_Event_Number ~ s(age), data = sd_cleanData_PE_events)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  gam_growthrate_plot(sd_cleanData_PE_events, gam.model, gam.growthrate, agevar = 'age', yvar = 'Gamma.Gamma_Event_Number', draw_points = F)
  
  # clean data for latency and power
  cleanData_latency_power <- alldata_TrialLevel_behavior %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2)) %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm = T) * 2))
  
  #aggregate Latency and power data
  avg_cleanData_latency_power <- merge(aggregate(. ~ Subject, cleanData_latency_power, mean), agefile, by = "Subject")
  sd_cleanData_latency_power <- merge(aggregate(. ~ Subject, cleanData_latency_power, sd), agefile, by = "Subject")
  
  # Gam model, accounting for age, looking a latency
  # these models are only to really look at one variable against age 
  gam.model <- gam(mgsLatency ~ s(age), data = avg_cleanData_latency_power)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  Latplot <- gam_growthrate_plot(avg_cleanData_latency_power, gam.model, gam.growthrate, agevar = 'age', yvar = 'mgsLatency', draw_points = F)
  
  # Gam model, accounting for age, looking a gamma power
  # these models are only to really look at one variable against age 
  gam.model <- gam(Gamma.log1p_Trial_Power ~ s(age), data = avg_cleanData_latency_power)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  powerPlot <- gam_growthrate_plot(avg_cleanData_latency_power, gam.model, gam.growthrate, agevar = 'age', yvar = 'Gamma.log1p_Trial_Power', draw_points = F)
  
  # Gam model, accounting for age, looking a latency variability
  # these models are only to really look at one variable against age 
  gam.model <- gam(mgsLatency ~ s(age), data = sd_cleanData_latency_power)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  Latplot <- gam_growthrate_plot(sd_cleanData_latency_power, gam.model, gam.growthrate, agevar = 'age', yvar = 'mgsLatency', draw_points = F)
  
  # Gam model, accounting for age, looking a gamma power variability
  # these models are only to really look at one variable against age 
  gam.model <- gam(Gamma.log1p_Trial_Power ~ s(age), data = sd_cleanData_latency_power)
  gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  powerPlot <- gam_growthrate_plot(sd_cleanData_latency_power, gam.model, gam.growthrate, agevar = 'age', yvar = 'Gamma.log1p_Trial_Power', draw_points = F)
  
  
}
