# Trying to find periods of significant growth 
library(LNCDR)
library(mgcv)
library(plotly)
library(graphics)

## Read in data
## merge trial level eeg and behavioral

# windows vs linux
analysis_root <- "H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis"
if(!dir.exists(analysis_root))
    analysis_root <- gsub('^H:/','/Volumes/Hera/', analysis_root)

# trial level data
# see 01.PrepData.R:Only_Take_One_Delay_Bin_TrialLevel()
eeg_trial_csv <- file.path(analysis_root, "alldata_TrialLevel_20210204.csv")
alldata_TrialLevel <- read.csv(eeg_trial_csv) #15008x23

#behavioral data contains each trial 
behave_csv <- file.path(analysis_root, 'Behavior_20200921.csv')
Behavior <- read.csv(behave_csv) # 15354x9
Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
Behavior$absPositionError <- abs(Behavior$PositionError)


alldata_TrialLevel_behavior <- merge(alldata_TrialLevel, Behavior,
                                     by = c("Subject", "Trial"))
# 
# > group_by(alldata_TrialLevel_behavior, Subject,visitno) %>% tally
# # A tibble: 144 x 3
# # Groups:   Subject [144]
#    Subject        visitno     n
#    <chr>            <int> <int>
#  1 10129_20180919       1    95
#  2 10173_20180801       1    96
#  3 10195_20180201       1    96
#  4 10202_20191001       1    92
#  ....


# First look at PE and gamma events 
# clean the data, removing outliers above 2 SDs
below_2sd <- function(x) abs(x - mean(x,na.rm=T)) < sd(x, na.rm=T)*2
# N.B. loosing 10 trials b/c PositionError calc after Gama_Event_Number is removed
cleanData_PE_events <- alldata_TrialLevel_behavior %>%
    filter(below_2sd(Gamma.Gamma_Event_Number)) %>%
    filter(below_2sd(PositionError)) 

# Gam model, accounting for age, length of delay period, and the fact that not all data points are unique, ie grouped by Subject 
# these models are only to really look at one variable against age 
gam.model <- cleanData_PE_events %>%
    mutate(Subject=as.factor(gsub('_','',Subject))) %>%
    gam(data=., formula=absPositionError ~ s(age) + Delay + s(Trial) + s(Subject, re="bs"))
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
PEplot <- gam_growthrate_plot(cleanData_PE_events,gam.model,gam.growthrate, agevar = 'age',yvar = 'absPositionError', idvar = 'Subject', draw_points = F)


# Gam model, accounting for age, length of delay period, and the fact that not all data points are unique, ie grouped by Subject 
# these models are only to really look at one variable against age 
gam.model <- gam(Gamma.Gamma_Event_Number ~ s(age) + Delay + s(Trial) + s(Subject, bs="re"), data = cleanData_PE_events)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
eventsPlot <- gam_growthrate_plot(cleanData_PE_events,gam.model,gam.growthrate, agevar = 'age',yvar = 'Gamma.Gamma_Event_Number', idvar = 'Subject', draw_points = F)

par(mfrow = c(2,1))
plot(PEplot$both)
plot(eventsPlot$both)


# compare to loess 
lunaize(ggplot(data = cleanData_PE_events, aes(x = age, y = Gamma.Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)'))

# clean data for PE and power
cleanData_PE_power <- alldata_TrialLevel_behavior %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power)) < (sd(Gamma.log1p_Trial_Power) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))

# gam model for power
gam.model <- gam(Gamma.log1p_Trial_Power ~ s(age) + Delay + s(Trial) + s(Subject, bs="re"), data = cleanData_PE_power)
gam.growthrate <- gam_growthrate(gam.model, agevar = 'age', idvar = 'Subject')
powerPlot <- gam_growthrate_plot(cleanData_PE_events,gam.model,gam.growthrate, agevar = 'age',yvar = 'Gamma.log1p_Trial_Power', idvar = 'Subject', draw_points = F)



