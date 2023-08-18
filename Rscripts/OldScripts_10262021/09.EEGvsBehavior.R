library(dplyr)
library(ggplot2)
library(lme4)
library(car)
library(LNCDR)

EventNumbervsBehavior_Ageplots <- function () {
setwd('~/Downloads/')
gamma <- readRDS(file = "gammadata.rds")
alpha <- readRDS(file = "alphaData.rds")
alpha_summary <-
  alpha %>% group_by(Alpha_Event_Number) %>% dplyr::summarize(
    meanLat = mean(mgsLatency),
    n = n(),
    sd = sd(mgsLatency),
    se = sd(mgsLatency) / sqrt(n())
  )
ggplot(data = alpha_summary[alpha_summary$Alpha_Event_Number <= 4,], aes(x =
                                                                           Alpha_Event_Number, y = meanLat)) +
  geom_bar(stat = "identity", position = "dodge") + geom_errorbar(aes(ymin =
                                                                        meanLat - se, ymax = meanLat + se)) + coord_cartesian(ylim = c(.4, .72))

ggplot(data = alpha_summary[alpha_summary$Alpha_Event_Number <= 3,], aes(x =
                                                                           Alpha_Event_Number, y = meanLat)) +
  geom_bar(stat = "identity", position = "dodge") + geom_errorbar(aes(ymin =
                                                                        meanLat - se, ymax = meanLat + se)) + coord_cartesian(ylim = c(.4, .46))
alpha$ageGroup <-
  ifelse(alpha$age < 16, '<16', ifelse(alpha$age >= 23, '>23', '16-23'))
alpha_summary_byage <-
  alpha %>% group_by(Alpha_Event_Number, ageGroup) %>% dplyr::summarize(
    meanLat = mean(mgsLatency),
    n = n(),
    sd = sd(mgsLatency),
    se = sd(mgsLatency) / sqrt(n())
  )
alpha_summary_byage$Alpha_Event_Number_fac <-
  as.factor(alpha_summary_byage$Alpha_Event_Number)
alpha_summary_byage$ageGroup <-
  factor(alpha_summary_byage$ageGroup, levels = c('<16', '16-23', '>23'))
ggplot(data = alpha_summary_byage[alpha_summary_byage$Alpha_Event_Number <= 3,],
       aes(x = ageGroup, y = meanLat, fill = Alpha_Event_Number_fac)) +
  geom_bar(stat = "identity", position = "dodge") + coord_cartesian(ylim =
                                                                      c(.35, .6))
lunaize(ggplot(
  data = alpha[alpha$Alpha_Event_Number <= 3,],
  aes(x = Alpha_Event_Number, y = mgsLatency, color = ageGroup)) +
  stat_smooth(method = 'lm')) + ggtitle("Increases in the number of alpha bursts are associated with increases in latency") + xlab("Alpha Event Number") + ylab("Response Latency (s)") + theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    axis.title.y = element_text(margin = margin(
      t = 0,
      r = 10,
      b = 0,
      l = 0
    )),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 20)
  ) +  scale_color_discrete("Age Group")


lunaize(ggplot(data = alpha[alpha$Alpha_Event_Number <= 3, ], aes(x = Alpha_Event_Number, y = mgsLatency)) + stat_smooth(method = 'lm'))
lmer.model <- lmer(mgsLatency ~ Alpha_Event_Number * age + (1 | Subject), data = alpha[alpha$Alpha_Event_Number <= 3, ])
car::Anova(lmer.model)

gamma_summary <-
  gamma %>% group_by(Gamma_Event_Number) %>% dplyr::summarize(
    meanAcc = mean(absPositionError),
    n = n(),
    sd = sd(absPositionError),
    se = sd(absPositionError) / sqrt(n())
  )

ggplot(data = gamma_summary[gamma_summary$Gamma_Event_Number <= 100, ], aes(x = Gamma_Event_Number, y = meanAcc)) +
  geom_bar(stat = "identity", position = "dodge") + geom_errorbar(aes(ymin = meanAcc - se, ymax = meanAcc + se)) + coord_cartesian(ylim = c(0, 10))

ggplot(data = gamma_summary[gamma_summary$Gamma_Event_Number <= 100, ], aes(x = Gamma_Event_Number, y = meanAcc)) +
  stat_smooth(method = 'lm') + coord_cartesian(ylim = c(0, 10))

gamma$ageGroup <-
  ifelse(gamma$age < 16, '<16', ifelse(gamma$age >= 23, '>23', '16-23'))
gamma_summary_byage <-
  gamma %>% group_by(Gamma_Event_Number, ageGroup) %>% dplyr::summarize(
    meanAcc = mean(absPositionError),
    n = n(),
    sd = sd(absPositionError),
    se = sd(absPositionError) / sqrt(n())
  )

gamma_summary_byage$Gamma_Event_Number_fac <-
  as.factor(gamma_summary_byage$Gamma_Event_Number)
gamma_summary_byage$ageGroup <-
  factor(gamma_summary_byage$ageGroup, levels = c('<16', '16-23', '>23'))

ggplot(data = gamma_summary_byage[gamma_summary_byage$Gamma_Event_Number <= 100,],
       aes(x = ageGroup, y = meanAcc, fill = Gamma_Event_Number_fac)) +
  stat_smooth(method = 'lm') + coord_cartesian(ylim = c(.35, .6))

lunaize(ggplot(
  data = gamma[gamma$Gamma_Event_Number <= 100,],
  aes(x = Gamma_Event_Number, y = absPositionError, color = ageGroup)
) +
  stat_smooth(method = 'lm'))  + ggtitle("Increases in the number of gamma bursts are associated with increases in position error") + xlab("Gamma Event Number") + ylab("abs(Position Error) (degs)") + theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    axis.title.y = element_text(margin = margin(
      t = 0,
      r = 10,
      b = 0,
      l = 0
    )),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 20)
  ) +  scale_color_discrete("Age Group")



lmer.model <-
  lmer(absPositionError ~ Gamma_Event_Number * age + (1 |
                                                        Subject), data = gamma[gamma$Gamma_Event_Number <= 100,])
car::Anova(lmer.model)
}

#Subject Level
EEGvsBehavior <- function(selected_vars, alldata) {
  library(matrixStats)
  
  alldf <- DelayOnly_Sublevel()
  alldata_delayOnly <- alldf$alldata_delayOnly
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  Behavior_SubLevel$absPositionError <- abs(Behavior_SubLevel$PositionError)

  eeg_behavior_subLevel <- merge(alldata_delayOnly, Behavior_SubLevel, by = "Subject")
  eeg_behavior_subLevel$Group <- factor(eeg_behavior_subLevel$Group)
  
  #PE var vs events var
  cleanData_sub_var_PE_events <- eeg_behavior_subLevel %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability)) < (sd(Gamma.Event_Number_Variability) * 2)) %>% filter(abs(PositionError_Variability - mean(PositionError_Variability, na.rm= T)) < (sd(PositionError_Variability, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = cleanData_sub_var_PE_events[], aes(x=Gamma.Event_Number_Variability, y=absPositionError_Variability, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Event Number Variability (Delay Only) vs Accuracy Variability") + xlab("Number of Gamma Events (sd)") + ylab("Accuracy (degs)(sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_PE_events[cleanData_sub_var_PE_events$Group == 1,], absPositionError_Variability ~ Gamma.Event_Number_Variability)
  print(car::Anova(lm.model))
  
  #events var vs age 
  lunaize(ggplot(data = cleanData_sub_var_PE_events[], aes(x=age, y=Gamma.Event_Number_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Event Number Variability (Delay Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events (sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_PE_events[], Gamma.Event_Number_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE var vs age 
  lunaize(ggplot(data = cleanData_sub_var_PE_events[], aes(x=age, y=absPositionError_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy Variability (Delay Only) vs Age") + xlab("Age") + ylab("Accuracy (degs) (sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_PE_events[], absPositionError_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  
  #PE var vs power var
  cleanData_sub_var_PE_power <- eeg_behavior_subLevel %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability)) < (sd(Gamma.log1p_Trial_Power_Variability) * 2)) %>% filter(abs(PositionError_Variability - mean(PositionError_Variability, na.rm= T)) < (sd(PositionError_Variability, na.rm = T) * 2))
  
  lunaize(ggplot(data = cleanData_sub_var_PE_power[], aes(x=Gamma.log1p_Trial_Power_Variability, y=absPositionError_Variability, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Power Variability (Delay Only) vs Accuracy Variability") + xlab("log(Gamma Power)(sd)") + ylab("Accuracy (degs)(sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_PE_power[cleanData_sub_var_PE_power$Group == 4,], absPositionError_Variability ~ Gamma.log1p_Trial_Power_Variability)
  print(car::Anova(lm.model))
  
  #power var vs age 
  lunaize(ggplot(data = cleanData_sub_var_PE_power[], aes(x=age, y=Gamma.log1p_Trial_Power_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Power Variability (Delay Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)(sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_PE_power[], Gamma.log1p_Trial_Power_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE var vs age 
  lunaize(ggplot(data = cleanData_sub_var_PE_power[], aes(x=age, y=absPositionError_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy Variability (Delay Only) vs Age") + xlab("Age") + ylab("Accuracy (degs) (sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_PE_power[], absPositionError_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  
  
  #latency var vs events var
  cleanData_sub_var_Lat_events <- eeg_behavior_subLevel %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability)) < (sd(Gamma.Event_Number_Variability) * 2)) %>% filter(abs(mgsLatency_Variability - mean(mgsLatency_Variability, na.rm= T)) < (sd(mgsLatency_Variability, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = cleanData_sub_var_Lat_events[], aes(x=Gamma.Event_Number_Variability, y=mgsLatency_Variability, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Event Number Variability (Delay Only) vs Latency Variability") + xlab("Number of Gamma Events (sd)") + ylab("Latency (secs)(sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_Lat_events[cleanData_sub_var_Lat_events$Group == 4,], mgsLatency_Variability ~ Gamma.Event_Number_Variability)
  print(car::Anova(lm.model))
  
  #events var vs age 
  lunaize(ggplot(data = cleanData_sub_var_Lat_events[], aes(x=age, y=Gamma.Event_Number_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Event Number Variability (Delay Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events (sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_Lat_events[], Gamma.Event_Number_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE var vs age 
  lunaize(ggplot(data = cleanData_sub_var_Lat_events[], aes(x=age, y=mgsLatency_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy Variability (Delay Only) vs Age") + xlab("Age") + ylab("Accuracy (degs) (sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_Lat_events[], mgsLatency_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  
  #Lat var vs power var
  cleanData_sub_var_Lat_power <- eeg_behavior_subLevel %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability)) < (sd(Gamma.log1p_Trial_Power_Variability) * 2)) %>% filter(abs(mgsLatency_Variability - mean(mgsLatency_Variability, na.rm= T)) < (sd(mgsLatency_Variability, na.rm = T) * 2))
  
  lunaize(ggplot(data = cleanData_sub_var_Lat_power[], aes(x=Gamma.log1p_Trial_Power_Variability, y=mgsLatency_Variability, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Power Variability (Delay Only) vs Latency Variability") + xlab("log(Gamma Power)(sd)") + ylab("Latency (secs)(sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_Lat_power[cleanData_sub_var_Lat_power$Group == 1,], mgsLatency_Variability ~ Gamma.log1p_Trial_Power_Variability)
  print(car::Anova(lm.model))
  
  #power var vs age 
  lunaize(ggplot(data = cleanData_sub_var_Lat_power[], aes(x=age, y=Gamma.log1p_Trial_Power_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Power Variability (Delay Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)(sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_Lat_power[], Gamma.log1p_Trial_Power_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE var vs age 
  lunaize(ggplot(data = cleanData_sub_var_Lat_power[], aes(x=age, y=mgsLatency_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy Variability (Delay Only) vs Age") + xlab("Age") + ylab("Accuracy (degs) (sd)"))
  
  lm.model <- lm(data = cleanData_sub_var_Lat_power[], mgsLatency_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
  #Duration vs Latency 
  cleanData_sub_duration_lat <- eeg_behavior_subLevel %>% filter(abs(Gamma.Event_Duration - mean(Gamma.Event_Duration)) < (sd(Gamma.Event_Duration) * 2)) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))
  
  lunaize(ggplot(data = cleanData_sub_duration_lat[], aes(x=Gamma.Event_Duration, y=mgsLatency)) + geom_point() + stat_smooth(method='lm', se = T) + ggtitle("Gamma Duration (Delay Only) vs Latency") + xlab("Gamma Event Duration") + ylab("Latency (secs)"))
  
  lm.model <- lm(data = cleanData_sub_duration_lat, mgsLatency ~ Gamma.Event_Duration)
  print(car::Anova(lm.model))
  
  #Duration vs Accuracy 
  cleanData_sub_duration_lat <- eeg_behavior_subLevel %>% filter(abs(Gamma.Event_Duration - mean(Gamma.Event_Duration)) < (sd(Gamma.Event_Duration) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))
  
  lunaize(ggplot(data = cleanData_sub_duration_lat[], aes(x=Gamma.Event_Duration, y=absPositionError)) + geom_point() + stat_smooth(method='lm', se = T) + ggtitle("Gamma Duration (Delay Only) vs Accuracy") + xlab("Gamma Event Duration") + ylab("Accuracy (degs)"))
  
  lm.model <- lm(data = cleanData_sub_duration_lat, PositionError ~ Gamma.Event_Duration)
  print(car::Anova(lm.model))
}

EEGvsBehavior_Rest <- function(selected_vars, alldata) {
  library(matrixStats)
  
  Resting_State_Data_SubjectLevel()
  rest_SubLevel <- alldf_SubjectLevel_RS$rest_SubLevel
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  Resteeg_behavior_subLevel <- merge(rest_SubLevel, Behavior_SubLevel, by = "Subject")
  Resteeg_behavior_subLevel$Group <- factor(Resteeg_behavior_subLevel$Group)
  
  #PE vs events
  cleanData_sub_PE_events <- Resteeg_behavior_subLevel %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number)) < (sd(Gamma.Event_Number) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = cleanData_sub_PE_events[], aes(x=Gamma.Event_Number, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Event Number (Rest Only) vs Accuracy") + xlab("Number of Gamma Events") + ylab("Accuracy (degs)"))
  
  lm.model <- lm(data = cleanData_sub_PE_events[cleanData_sub_PE_events$Group == 4,], absPositionError ~ Gamma.Event_Number)
  print(car::Anova(lm.model))
  
  #events vs age 
  lunaize(ggplot(data = cleanData_sub_PE_events[], aes(x=age, y=Gamma.Event_Number)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Event Number (Rest Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events"))
  
  lm.model <- lm(data = cleanData_sub_PE_events[], Gamma.Event_Number ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE vs age 
  lunaize(ggplot(data = cleanData_sub_PE_events[], aes(x=age, y=absPositionError)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy (Rest Only) vs Age") + xlab("Age") + ylab("Accuracy (degs)"))
  
  lm.model <- lm(data = cleanData_sub_PE_events[], absPositionError ~ inverseAge)
  print(car::Anova(lm.model))
  
  
  #PE vs power 
  cleanData_sub_PE_power <- Resteeg_behavior_subLevel %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power)) < (sd(Gamma.log1p_Trial_Power) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))
  
  lunaize(ggplot(data = cleanData_sub_PE_power[], aes(x=Gamma.log1p_Trial_Power, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Power (Rest Only) vs Accuracy") + xlab("log(Gamma Power)") + ylab("Accuracy (degs)"))
  
  lm.model <- lm(data = cleanData_sub_PE_power[cleanData_sub_PE_power$Group == 1,], absPositionError ~ Gamma.log1p_Trial_Power)
  print(car::Anova(lm.model))
  
  #power vs age 
  lunaize(ggplot(data = cleanData_sub_PE_power[], aes(x=age, y=Gamma.log1p_Trial_Power)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Power (Rest Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)"))
  
  lm.model <- lm(data = cleanData_sub_PE_power[], Gamma.log1p_Trial_Power ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE vs age 
  lunaize(ggplot(data = cleanData_sub_PE_power[], aes(x=age, y=absPositionError)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy (Rest Only) vs Age") + xlab("Age") + ylab("Accuracy (degs)"))
  
  lm.model <- lm(data = cleanData_sub_PE_power[], absPositionError ~ inverseAge)
  print(car::Anova(lm.model))
  
  
  
  #latency vs events 
  cleanData_sub_Lat_events <- Resteeg_behavior_subLevel %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number)) < (sd(Gamma.Event_Number) * 2)) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = cleanData_sub_Lat_events[], aes(x=Gamma.Event_Number, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Event Number (Rest Only) vs Latency") + xlab("Number of Gamma Events") + ylab("Latency (secs)"))
  
  lm.model <- lm(data = cleanData_sub_Lat_events[cleanData_sub_Lat_events$Group == 1,], mgsLatency ~ Gamma.Event_Number)
  print(car::Anova(lm.model))
  
  #events vs age 
  lunaize(ggplot(data = cleanData_sub_Lat_events[], aes(x=age, y=Gamma.Event_Number)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Event Number (Rest Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events"))
  
  lm.model <- lm(data = cleanData_sub_Lat_events[], Gamma.Event_Number ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE vs age 
  lunaize(ggplot(data = cleanData_sub_Lat_events[], aes(x=age, y=mgsLatency)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy (Rest Only) vs Age") + xlab("Age") + ylab("Accuracy (degs)"))
  
  lm.model <- lm(data = cleanData_sub_Lat_events[], mgsLatency ~ inverseAge)
  print(car::Anova(lm.model))
  
  
  #Lat vs power 
  cleanData_sub_Lat_power <- Resteeg_behavior_subLevel %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power)) < (sd(Gamma.log1p_Trial_Power) * 2)) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))
  
  lunaize(ggplot(data = cleanData_sub_Lat_power[], aes(x=Gamma.log1p_Trial_Power, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Power (Rest Only) vs Latency") + xlab("log(Gamma Power)") + ylab("Latency (secs)"))
  
  lm.model <- lm(data = cleanData_sub_Lat_power[cleanData_sub_Lat_power$Group == 4,], mgsLatency ~ Gamma.log1p_Trial_Power)
  print(car::Anova(lm.model))
  
  #power vs age 
  lunaize(ggplot(data = cleanData_sub_Lat_power[], aes(x=age, y=Gamma.log1p_Trial_Power)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Gamma Power (Rest Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)"))
  
  lm.model <- lm(data = cleanData_sub_Lat_power[], Gamma.log1p_Trial_Power ~ inverseAge)
  print(car::Anova(lm.model))
  
  #PE  vs age 
  lunaize(ggplot(data = cleanData_sub_Lat_power[], aes(x=age, y=mgsLatency_Variability)) + stat_smooth(method='lm', formula = 'y~I(1/x)') + ggtitle("Accuracy Variability (Delay Only) vs Age") + xlab("Age") + ylab("Accuracy (degs) (sd)"))
  
  lm.model <- lm(data = cleanData_sub_Lat_power[], mgsLatency_Variability ~ inverseAge)
  print(car::Anova(lm.model))
  
}


DelayMinusFixvsBehavior_TrialPower <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  EEG_and_Behavior <- merge(alldata_DelayminusFix, Behavior_SubLevel, by = "Subject")
  
  # Trial Power
  
  avgGamma <- aggregate(Gamma.Trial_Power ~ visitno , EEG_and_Behavior, mean)
  sdGamma <- aggregate(Gamma.Trial_Power ~ visitno , EEG_and_Behavior, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Trial_Power ~ visitno , EEG_and_Behavior, mean)
  sdBeta <- aggregate(Beta.Trial_Power ~ visitno , EEG_and_Behavior, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Trial_Power ~ visitno , EEG_and_Behavior, mean)
  sdAlpha <- aggregate(Alpha.Trial_Power ~ visitno , EEG_and_Behavior, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Trial_Power ~ visitno , EEG_and_Behavior, mean)
  sdTheta <- aggregate(Theta.Trial_Power ~ visitno , EEG_and_Behavior, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2] 
  
  # Gamma 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Trial_Power < gammaCutoff,], aes(x = absPositionError, y = Gamma.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Gamma Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power (Delay-Fix)"))

  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Trial_Power < gammaCutoff,], aes(x = mgsLatency, y = Gamma.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Gamma Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power (Delay-Fix)"))
  
  # Beta 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Trial_Power < BetaCutoff,], aes(x = absPositionError, y = Beta.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Beta Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Trial_Power < BetaCutoff,], aes(x = mgsLatency, y = Beta.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Beta Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power (Delay-Fix)"))
  
  
  # Alpha
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Trial_Power < AlphaCutoff,], aes(x = absPositionError, y = Alpha.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Alpha Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Trial_Power < AlphaCutoff,], aes(x = mgsLatency, y = Alpha.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Alpha Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Trial_Power < AlphaCutoff,],  Alpha.Trial_Power ~ mgsLatency)
  print(anova(lm.model))
  
  #Theta
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Trial_Power < ThetaCutoff,], aes(x = absPositionError, y = Theta.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Theta Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Trial_Power < ThetaCutoff,], aes(x = mgsLatency, y = Theta.Trial_Power))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Theta Trial Power (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power (Delay-Fix)"))
  
  
  
}

DelayMinusFixvsBehavior_TrialPowerVar <- function() {
  
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  EEG_and_Behavior <- merge(alldata_DelayminusFix, Behavior_SubLevel, by = "Subject")
  
  
  # Trial Power Variability 
  
  avgGamma <- aggregate(Gamma.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdGamma <- aggregate(Gamma.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdBeta <- aggregate(Beta.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdAlpha <- aggregate(Alpha.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdTheta <- aggregate(Theta.Trial_Power_Variability ~ visitno , alldata_DelayminusFix, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  # Gamma 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Trial_Power_Variability < gammaCutoff,], aes(x = absPositionError, y = Gamma.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Gamma Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Trial_Power_Variability < gammaCutoff,], aes(x = mgsLatency, y = Gamma.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Gamma Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power Variability (Delay-Fix)"))
  
  # Beta 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Trial_Power_Variability < BetaCutoff,], aes(x = absPositionError, y = Beta.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Beta Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Trial_Power_Variability < BetaCutoff,], aes(x = mgsLatency, y = Beta.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Beta Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power Variability (Delay-Fix)"))
  
  
  # Alpha
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Trial_Power_Variability < AlphaCutoff,], aes(x = absPositionError, y = Alpha.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Alpha Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Trial_Power_Variability < AlphaCutoff,], aes(x = mgsLatency, y = Alpha.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Alpha Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power Variability (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Trial_Power_Variability < AlphaCutoff,],  Alpha.Trial_Power_Variability ~ mgsLatency)
  print(anova(lm.model))
  
  #Theta
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Trial_Power_Variability < ThetaCutoff,], aes(x = absPositionError, y = Theta.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Theta Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Trial Power Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Trial_Power_Variability < ThetaCutoff,], aes(x = mgsLatency, y = Theta.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Theta Trial Power Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Trial Power Variability (Delay-Fix)"))
  
  
}

DelayMinusFixvsBehavior_EventNumber <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  EEG_and_Behavior <- merge(alldata_DelayminusFix, Behavior_SubLevel, by = "Subject")
  
  
  avgGamma <- aggregate(Gamma.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdGamma <- aggregate(Gamma.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdBeta <- aggregate(Beta.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdAlpha <- aggregate(Alpha.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdTheta <- aggregate(Theta.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  # Gamma 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number < gammaCutoff,], aes(x = absPositionError, y = Gamma.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Gamma Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number < gammaCutoff,],  Gamma.Event_Number ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number < gammaCutoff,], aes(x = mgsLatency, y = Gamma.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Gamma Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number < gammaCutoff,],  Gamma.Event_Number ~ mgsLatency)
  print(anova(lm.model))
  
  # Beta 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Number < BetaCutoff,], aes(x = absPositionError, y = Beta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Beta Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Number < BetaCutoff,], aes(x = mgsLatency, y = Beta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Beta Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  
  
  # Alpha
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Number < AlphaCutoff,], aes(x = absPositionError, y = Alpha.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Alpha Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Number < AlphaCutoff,], aes(x = mgsLatency, y = Alpha.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Alpha Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  
  
  #Theta
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Number < ThetaCutoff,], aes(x = absPositionError, y = Theta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Theta Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Number < ThetaCutoff,], aes(x = mgsLatency, y = Theta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Theta Event Number (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  
  
  
}

DelayMinusFixvsBehavior_EventNumberVar <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  EEG_and_Behavior <- merge(alldata_DelayminusFix, Behavior_SubLevel, by = "Subject")
  
  
  avgGamma <- aggregate(Gamma.Event_Number_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdGamma <- aggregate(Gamma.Event_Number_Variability ~ visitno , alldata_DelayminusFix, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Number_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdBeta <- aggregate(Beta.Event_Number_Variability ~ visitno , alldata_DelayminusFix, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Number_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdAlpha <- aggregate(Alpha.Event_Number_Variability ~ visitno , alldata_DelayminusFix, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Number_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdTheta <- aggregate(Theta.Event_Number_Variability ~ visitno , alldata_DelayminusFix, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  # Gamma 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number_Variability < gammaCutoff,], aes(x = absPositionError, y = Gamma.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Gamma Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number Variability (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number_Variability < gammaCutoff,],  Gamma.Event_Number_Variability ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number_Variability < gammaCutoff,], aes(x = mgsLatency, y = Gamma.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Gamma Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number Variability (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number_Variability < gammaCutoff,],  Gamma.Event_Number_Variability ~ mgsLatency)
  print(anova(lm.model))
  
  # Beta 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Number_Variability < BetaCutoff,], aes(x = absPositionError, y = Beta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Beta Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Number_Variability < BetaCutoff,], aes(x = mgsLatency, y = Beta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Beta Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number Variability (Delay-Fix)"))
  
  
  # Alpha
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Number_Variability < AlphaCutoff,], aes(x = absPositionError, y = Alpha.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Alpha Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Number_Variability < AlphaCutoff,], aes(x = mgsLatency, y = Alpha.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Alpha Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number Variability (Delay-Fix)"))
  
  
  #Theta
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Number_Variability < ThetaCutoff,], aes(x = absPositionError, y = Theta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Theta Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Number_Variability < ThetaCutoff,], aes(x = mgsLatency, y = Theta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Theta Event Number Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number Variability (Delay-Fix)"))
  
  
  
}

DelayMinusFixvsBehavior_EventDuration <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  EEG_and_Behavior <- merge(alldata_DelayminusFix, Behavior_SubLevel, by = "Subject")
  
  
  avgGamma <- aggregate(Gamma.Event_Duration ~ visitno , alldata_DelayminusFix, mean)
  sdGamma <- aggregate(Gamma.Event_Duration ~ visitno , alldata_DelayminusFix, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Duration ~ visitno , alldata_DelayminusFix, mean)
  sdBeta <- aggregate(Beta.Event_Duration ~ visitno , alldata_DelayminusFix, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Duration ~ visitno , alldata_DelayminusFix, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration ~ visitno , alldata_DelayminusFix, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Duration ~ visitno , alldata_DelayminusFix, mean)
  sdTheta <- aggregate(Theta.Event_Duration ~ visitno , alldata_DelayminusFix, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  # Gamma 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration < gammaCutoff,], aes(x = absPositionError, y = Gamma.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Gamma Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration < gammaCutoff,],  Gamma.Event_Duration ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration < gammaCutoff,], aes(x = mgsLatency, y = Gamma.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Gamma Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration < gammaCutoff,],  Gamma.Event_Duration ~ mgsLatency)
  print(anova(lm.model))
  
  # Beta 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Duration < BetaCutoff,], aes(x = absPositionError, y = Beta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Beta Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Duration < BetaCutoff,], aes(x = mgsLatency, y = Beta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Beta Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration (Delay-Fix)"))
  
  
  # Alpha
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration < AlphaCutoff,], aes(x = absPositionError, y = Alpha.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Alpha Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration (Delay-Fix)"))
  
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration < AlphaCutoff,],  Alpha.Event_Duration ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration < AlphaCutoff,], aes(x = mgsLatency, y = Alpha.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Alpha Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration (Delay-Fix)"))
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration < AlphaCutoff,],  Alpha.Event_Duration ~ mgsLatency)
  print(anova(lm.model))
  
  #Theta
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Duration < ThetaCutoff,], aes(x = absPositionError, y = Theta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Theta Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Duration < ThetaCutoff,], aes(x = mgsLatency, y = Theta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Theta Event Duration (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration (Delay-Fix)"))
  
  
  
}

DelayMinusFixvsBehavior_EventDurationVar <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  EEG_and_Behavior <- merge(alldata_DelayminusFix, Behavior_SubLevel, by = "Subject")
  
  
  avgGamma <- aggregate(Gamma.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdGamma <- aggregate(Gamma.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdBeta <- aggregate(Beta.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, mean)
  sdTheta <- aggregate(Theta.Event_Duration_Variability ~ visitno , alldata_DelayminusFix, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  # Gamma 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration_Variability < gammaCutoff,], aes(x = absPositionError, y = Gamma.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Gamma Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration Variability (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration_Variability < gammaCutoff,],  Gamma.Event_Duration_Variability ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration_Variability < gammaCutoff,], aes(x = mgsLatency, y = Gamma.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Gamma Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration Variability (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Duration_Variability < gammaCutoff,],  Gamma.Event_Duration_Variability ~ mgsLatency)
  print(anova(lm.model))
  
  # Beta 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Duration_Variability < BetaCutoff,], aes(x = absPositionError, y = Beta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Beta Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Duration_Variability < BetaCutoff,], aes(x = mgsLatency, y = Beta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Beta Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration Variability (Delay-Fix)"))
  
  
  # Alpha
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration_Variability < AlphaCutoff,], aes(x = absPositionError, y = Alpha.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Alpha Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration Variability (Delay-Fix)"))
  
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration_Variability < AlphaCutoff,],  Alpha.Event_Duration_Variability ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration_Variability < AlphaCutoff,], aes(x = mgsLatency, y = Alpha.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Alpha Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration Variability (Delay-Fix)"))
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Duration_Variability < AlphaCutoff,],  Alpha.Event_Duration_Variability ~ mgsLatency)
  print(anova(lm.model))
  
  #Theta
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Duration_Variability < ThetaCutoff,], aes(x = absPositionError, y = Theta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Theta Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Duration Variability (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Duration_Variability < ThetaCutoff,], aes(x = mgsLatency, y = Theta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Theta Event Duration Variability (Delay-Fix)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Duration Variability (Delay-Fix)"))
  
  
  
}


DelayMinusFixvsBehavior_EventNumber_AdultsOnly <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  EEG_and_Behavior <- merge(alldata_DelayminusFix, Behavior_SubLevel, by = "Subject")
  
  
  avgGamma <- aggregate(Gamma.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdGamma <- aggregate(Gamma.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdBeta <- aggregate(Beta.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdAlpha <- aggregate(Alpha.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Number ~ visitno , alldata_DelayminusFix, mean)
  sdTheta <- aggregate(Theta.Event_Number ~ visitno , alldata_DelayminusFix, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  # Gamma 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number < gammaCutoff  & EEG_and_Behavior$age > 20.99,], aes(x = absPositionError, y = Gamma.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Gamma Event Number (Delay-Fix) Adults Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number < gammaCutoff & EEG_and_Behavior$age > 20.99,],  Gamma.Event_Number ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Gamma.Event_Number < gammaCutoff & EEG_and_Behavior$age > 20.99,], aes(x = mgsLatency, y = Gamma.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Gamma Event Number (Delay-Fix) Adults Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  

  # Beta 
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Number < BetaCutoff & EEG_and_Behavior$age > 20.99,], aes(x = absPositionError, y = Beta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Beta Event Number (Delay-Fix) Adults Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Number < BetaCutoff & EEG_and_Behavior$age > 20.99,],  Beta.Event_Number ~ absPositionError)
  print(anova(lm.model))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Beta.Event_Number < BetaCutoff & EEG_and_Behavior$age > 20.99,], aes(x = mgsLatency, y = Beta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Beta Event Number (Delay-Fix) Adults Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  
  
  # Alpha
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Number < AlphaCutoff & EEG_and_Behavior$age > 20.99,], aes(x = absPositionError, y = Alpha.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Alpha Event Number (Delay-Fix) Adult Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Alpha.Event_Number < AlphaCutoff& EEG_and_Behavior$age > 20.99,], aes(x = mgsLatency, y = Alpha.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Alpha Event Number (Delay-Fix) Adult Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  
  
  #Theta
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Number < ThetaCutoff & EEG_and_Behavior$age > 20.99,], aes(x = absPositionError, y = Theta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Absolute Position Error vs Theta Event Number (Delay-Fix) Adults Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Position Error") + ylab("Event Number (Delay-Fix)"))
  
  print(ggplot(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Number < ThetaCutoff & EEG_and_Behavior$age > 20.99,], aes(x = mgsLatency, y = Theta.Event_Number))  + geom_point() + stat_smooth(method = 'lm') + ggtitle("Latency vs Theta Event Number (Delay-Fix) Adults Only") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Latency") + ylab("Event Number (Delay-Fix)"))
  
  
  lm.model <- lm(data = EEG_and_Behavior[EEG_and_Behavior$visitno < 2 & EEG_and_Behavior$Theta.Event_Number < ThetaCutoff & EEG_and_Behavior$age > 20.99,],  Theta.Event_Number ~ mgsLatency)
  print(anova(lm.model))
}




#Trial Level
DelayvsBehavior_TrialLevel_zScored <- function() {
  
  alldf <- Only_Take_One_Delay_Bin_TrialLevel()
  alldata_TrialLevel <- alldf$alldata_TrialLevel
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv')

  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  alldata_TrialLevel <- merge(alldata_TrialLevel, Behavior, by = c("Subject", "Trial"))
  alldata_TrialLevel$Group <- factor(alldata_TrialLevel$Group)
  
  #OUTLIER DETECTION
  alldata_TrialLevel_outlierDetection  <- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2)) %>% filter(abs(Gamma.Gamma_Event_Number - mean(Gamma.Gamma_Event_Number, na.rm= T)) < (sd(Gamma.Gamma_Event_Number, na.rm = T) * 2))
  
  accuracyOutliers  <- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) > (sd(PositionError, na.rm = T) * 2)) 
  
  
  accuracyOutliers <- merge(accuracyOutliers, agefile, by = "Subject")
  accuracyOutliers %>% group_by(Group) %>% tally()
  
  
 gammaNumberOutliers <- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(Gamma.Gamma_Event_Number - mean(Gamma.Gamma_Event_Number, na.rm= T)) > (sd(Gamma.Gamma_Event_Number, na.rm = T) * 2))
  
  
 gammaNumberOutliers <- merge(gammaNumberOutliers, agefile, by = "Subject")
 gammaNumberOutliers %>% group_by(Group) %>% tally()
  
 
 
 # ZSCORING
  alldata_TrialLevel_notZscored <- merge(alldata_TrialLevel_outlierDetection, agefile, by = "Subject")
  alldata_TrialLevel_notZscored$Group <- factor(alldata_TrialLevel_notZscored$Group)
  
  alldata_TrialLevel_zScored <- alldata_TrialLevel_outlierDetection[] %>% group_by(Subject) %>%  mutate_if(is.numeric, ~ as.numeric(scale(., center=T)))
  
  alldata_TrialLevel_zScored <- merge(alldata_TrialLevel_zScored, agefile, by = "Subject")
  
  alldata_TrialLevel_zScored$Group <- factor(alldata_TrialLevel_zScored$Group)
  
  gammaNEvsPE <- lmer(absPositionError ~ Gamma.Gamma_Event_Number * age + (1|Subject), data = alldata_TrialLevel_zScored)
  summary(gammaNEvsPE)
  car::Anova(gammaNEvsPE)
  
  
  ggplot(data = alldata_TrialLevel_zScored[], aes(x=Gamma.Gamma_Event_Number, y=absPositionError)) + stat_smooth(method='lm') + ggtitle("Gamma Event Number (Delay Only) vs Position Error")
  lm.model <- lm(data = alldata_TrialLevel_zScored[], absPositionError ~ Gamma.Gamma_Event_Number)
  print(anova(lm.model))
  
  
  ggplot(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$age >= 25,], aes(x=Gamma.Gamma_Event_Number, y=absPositionError)) + stat_smooth(method='lm') + ggtitle("Gamma Event Number (Delay Only) vs Position Error: Adults Only")
  lm.model <- lm(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$age >= 25,], absPositionError ~ Gamma.Gamma_Event_Number)
  print(anova(lm.model))
  
  gammaTPvsPE <- lmer(absPositionError ~ Gamma.Gamma_Trial_Power * age + (1|Subject), data = alldata_TrialLevel_zScored)
  summary(gammaTPvsPE)
  car::Anova(gammaTPvsPE)
}

DelayvsBehavior_TrialLevel <- function() {
  
  alldf <- Only_Take_One_Delay_Bin_TrialLevel()
  alldata_TrialLevel <- alldf$alldata_TrialLevel
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv')
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  alldata_TrialLevel <- merge(alldata_TrialLevel, Behavior, by = c("Subject", "Trial"))
  alldata_TrialLevel$Group <- factor(alldata_TrialLevel$Group)
  
  # Gamma event number 
  alldata_TrialLevel_outlierDetection  <- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2)) %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm = T) * 2))
  
  # Gamma duration vs position error
  alldata_TrialLevel_outlierDetection  <- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2)) %>% filter(abs(Gamma.Event_Duration - mean(Gamma.Event_Duration, na.rm= T)) < (sd(Gamma.Event_Duration, na.rm = T) * 2))
  

  lunaize(ggplot(data = alldata_TrialLevel_outlierDetection, aes(x=Gamma.Event_Duration, y=absPositionError)) + stat_smooth(method='lm', se = T) + ggtitle("Gamma Event \nDuration (Delay Only) vs Accuracy") + xlab("Duration of Gamma Events") + ylab("Accuracy (degs)")) + theme(plot.title = element_text(hjust = 0.5))
  
  lm.model <- lmer(data = alldata_TrialLevel_outlierDetection, absPositionError ~ Gamma.Event_Duration + (1|Subject))
  print(car::Anova(lm.model))
  
  # Gamma duration vs latency
  alldata_TrialLevel_outlierDetection  <- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2)) %>% filter(abs(Gamma.Event_Duration - mean(Gamma.Event_Duration, na.rm= T)) < (sd(Gamma.Event_Duration, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = alldata_TrialLevel_outlierDetection, aes(x=Gamma.Event_Duration, y=mgsLatency)) + stat_smooth(method='lm', se = T) + ggtitle("Gamma Event \nDuration (Delay Only) vs Latency") + xlab("Duration of Gamma Events") + ylab("Latency (secs)")) + theme(plot.title = element_text(hjust = 0.5))
  
  lm.model <- lmer(data = alldata_TrialLevel_outlierDetection, mgsLatency ~ Gamma.Event_Duration + (1|Subject))
  print(car::Anova(lm.model))
  
}


DelayvsBehavior_TrialLevel_AgeGroups <- function() {
  
  # Look at individual age groups dealy vs behavior relationships 
   Only_Take_One_Delay_Bin_TrialLevel()
  alldata_TrialLevel <- alldf$alldata_TrialLevel
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv')
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  alldata_TrialLevel <- merge(alldata_TrialLevel, Behavior, by = c("Subject", "Trial"))
  alldata_TrialLevel$Group <- factor(alldata_TrialLevel$Group)
 
  #number of events vs position error
  lunaize(ggplot(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$visitno <2,], aes(x=Gamma.Gamma_Event_Number, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = T) + ggtitle("Gamma Event Number (Delay Only) vs Accuracy") + xlab("Number of Gamma Events") + ylab("Accuracy") + scale_fill_brewer(palette = "BuGn"))
  
  lm.model <- lmer(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$Group == 4,], absPositionError ~ Gamma.Gamma_Event_Number + (1|Subject))
  print(car::anova(lm.model))
  
  #number of events vs position error, NOT ZSCORED WITH OUTLIER DETECTION
  
  #OUTLIER DETECTION
  alldata_TrialLevel_CleanPE_events<- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2)) %>% filter(abs(Gamma.Gamma_Event_Number - mean(Gamma.Gamma_Event_Number, na.rm= T)) < (sd(Gamma.Gamma_Event_Number, na.rm = T) * 2))
  
  lunaize(ggplot(data = alldata_TrialLevel_CleanPE_events, aes(x=Gamma.Gamma_Event_Number, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = F) + ggtitle("Gamma Event Number (Delay Only) vs Accuracy") + xlab("Number of Gamma Events") + ylab("Accuracy (degs)"))
  
  
  lm.model <- lmer(data = alldata_TrialLevel_CleanPE_events[alldata_TrialLevel_CleanPE_events$Group == 4,], absPositionError ~ Gamma.Gamma_Event_Number + (1|Subject))
  print(car::Anova(lm.model))
  
  #plot gamma events by age
  
  lunaize(ggplot(data = alldata_TrialLevel_CleanPE_events, aes(x=age, y=Gamma.Gamma_Event_Number)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("Gamma Event Number (Delay Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events"))
  
  lm.model <- lmer(data = alldata_TrialLevel_CleanPE_events[], Gamma.Gamma_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  #trial power vs position error 
  lunaize(ggplot(data = alldata_TrialLevel_zScored[], aes(x=Gamma.Gamma_Trial_Power, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Event Power (Delay Only) vs Accuracy") + xlab("Power") + ylab("Accuracy") + scale_fill_brewer(palette = "BuGn"))
  
  
  lm.model <- lm(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$Group == 3,], absPositionError ~ Gamma.Gamma_Trial_Power)
  print(anova(lm.model))
  
  #trial power vs PE, NOT ZSCORED WITH OUTLIER DETECTION 
  alldata_TrialLevel_CleanPE_Power<- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2)) %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm = T) * 2))
  
  lunaize(ggplot(data = alldata_TrialLevel_CleanPE_Power, aes(x=Gamma.log1p_Trial_Power, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = F) + ggtitle("Gamma Power (Delay Only) vs Accuracy") + xlab("log(Gamma Power)") + ylab("Accuracy (degs)"))
  
  lm.model <- lmer(data = alldata_TrialLevel_CleanPE_Power[alldata_TrialLevel_CleanPE_Power$Group == 1,], absPositionError ~ Gamma.log1p_Trial_Power + (1|Subject))
  print(car::Anova(lm.model))
  
  #plot gamma events by age
  
  lunaize(ggplot(data = alldata_TrialLevel_CleanPE_Power, aes(x=age, y=Gamma.log1p_Trial_Power)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("Gamma Power (Delay Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)"))
  
  lm.model <- lmer(data = alldata_TrialLevel_CleanPE_Power[], Gamma.log1p_Trial_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  #mgs latency
  lunaize(ggplot(data = alldata_TrialLevel_zScored[], aes(x=Gamma.Gamma_Event_Number, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Event Number (Delay Only) vs mgs Latency") + xlab("Number of Gamma Events") + ylab("Latency"))
  
  lm.model <- lm(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$Group == 1,], mgsLatency ~ Gamma.Gamma_Event_Number)
  print(anova(lm.model))
  
  #trial power vs MGS Latency, NOT ZSCORED WITH OUTLIER DETECTION 
  alldata_TrialLevel_Clean_Lat_Power<- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2)) %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm = T) * 2))
  
  lunaize(ggplot(data = alldata_TrialLevel_Clean_Lat_Power, aes(x=Gamma.log1p_Trial_Power, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = F) + ggtitle("Gamma Power (Delay Only) vs MGS Latency") + xlab("log(Gamma Power)") + ylab("Latency (sec)"))
  
  lm.model <- lmer(data = alldata_TrialLevel_Clean_Lat_Power[alldata_TrialLevel_Clean_Lat_Power$Group == 4,], mgsLatency ~ Gamma.log1p_Trial_Power + (1|Subject))
  print(car::Anova(lm.model))
  
  #plot gamma events by age
  
  lunaize(ggplot(data = alldata_TrialLevel_Clean_Lat_Power, aes(x=age, y=Gamma.log1p_Trial_Power)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("Gamma Power (Delay Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)"))
  
  lm.model <- lmer(data = alldata_TrialLevel_Clean_Lat_Power[], Gamma.log1p_Trial_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  #latency vs power, zscored
  lunaize(ggplot(data = alldata_TrialLevel_zScored[], aes(x=Gamma.Gamma_Trial_Power, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("Gamma Power (Delay Only) vs mgs Latency") + xlab("Power") + ylab("Latency"))
  
  lm.model <- lm(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$Group == 4,], mgsLatency ~ Gamma.Gamma_Trial_Power)
  print(anova(lm.model))
  
  #gamma events vs MGS Latency, NOT ZSCORED WITH OUTLIER DETECTION 
  alldata_TrialLevel_Clean_Lat_Events<- alldata_TrialLevel %>% group_by(Subject) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2)) %>% filter(abs(Gamma.Gamma_Event_Number - mean(Gamma.Gamma_Event_Number, na.rm= T)) < (sd(Gamma.Gamma_Event_Number, na.rm = T) * 2))
  
  lunaize(ggplot(data = alldata_TrialLevel_Clean_Lat_Events, aes(x=Gamma.Gamma_Event_Number, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = F) + ggtitle("Number of Gamma Events (Delay Only) vs MGS Latency") + xlab("Number of Gamma Events") + ylab("Latency (sec)"))
  
  lm.model <- lmer(data = alldata_TrialLevel_Clean_Lat_Events[alldata_TrialLevel_Clean_Lat_Events$Group == 1,], mgsLatency ~ Gamma.Gamma_Event_Number + (1|Subject))
  print(car::Anova(lm.model))
  
  #plot gamma events by age
  
  lunaize(ggplot(data = alldata_TrialLevel_Clean_Lat_Events, aes(x=age, y=Gamma.Gamma_Event_Number)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("Number of Gamma Events (Delay Only) vs Age") + xlab("Age") + ylab("Number of Events"))
  
  lm.model <- lmer(data = alldata_TrialLevel_Clean_Lat_Events[], Gamma.Gamma_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  
  # Brendans method to see where the age effect actually shows up
  install.packages("tvem")
  library(tvem)
  
  alldata_TrialLevel <- subset(alldata_TrialLevel, select = -c(LunaID,ScanDate,calR2))
  alldata_TrialLevel$Gamma.Gamma_Event_Number <- as.vector(alldata_TrialLevel$Gamma.Gamma_Event_Number)
  alldata_TrialLevel <- alldata_TrialLevel[complete.cases(alldata_TrialLevel),]
  
  alldata_TrialLevel_zScored <- subset(alldata_TrialLevel_zScored, select = -c(LunaID,ScanDate,calR2))
  alldata_TrialLevel_zScored$Gamma.Gamma_Event_Number <- as.vector(alldata_TrialLevel_zScored$Gamma.Gamma_Event_Number)
  alldata_TrialLevel_zScored <- alldata_TrialLevel_zScored[complete.cases(alldata_TrialLevel_zScored),]
  
  test <- subset(alldata_TrialLevel, select = c(Subject, Gamma.Gamma_Event_Number, age, absPositionError))
  
  
  tvemmodelobj <- tvem(data = test,
                       formula = absPositionError ~ Gamma.Gamma_Event_Number,
                       id = Subject,
                       time = age)
  
  bammodel <- tvemmodelobj$back_end_model
  bammodelsum <- summary(bammodel)
  pvi <- grep("primaryvar",row.names(bammodelsum$s.table))
  bammodelpval <- bammodelsum$s.table[pvi,4]
  
  tvem_grid_fitted_coefficients<-tvemmodelobj$grid_fitted_coefficients$primaryvar
  tvem_grid_fitted_coefficients$time_grid<-tvemmodelobj$time_grid
  
}


FixvsBehavior_TrialLevel <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  alldata_Fix_TrialLevel <- merge(alldata_Fix_TrialLevel[1:14], Behavior, by = c("Subject", "Trial"), all.x = TRUE, all.y = TRUE)
  
  alldata_Fix_TrialLevel_zScored <- alldata_Fix_TrialLevel %>% group_by(Subject) %>% mutate_if(is.numeric, scale, center = T)
  
  alldata_Fix_TrialLevel_zScored <- merge(alldata_Fix_TrialLevel_zScored, agefile, by = "Subject")
  
  gammaNEvsPE <- lmer(absPositionError ~ Gamma_Event_Number * age + (1|Subject), data = alldata_Fix_TrialLevel_zScored)
  summary(gammaNEvsPE)
  car::Anova(gammaNEvsPE)
  
  
  ggplot(data = alldata_Fix_TrialLevel_zScored[], aes(x=Gamma_Event_Number, y=absPositionError)) + stat_smooth(method='lm') + ggtitle("Gamma Event Number (Fix Only) vs Position Error")
  lm.model <- lm(data = alldata_Fix_TrialLevel_zScored[], absPositionError ~ Gamma_Event_Number)
  print(anova(lm.model))
  
  
  ggplot(data = alldata_Fix_TrialLevel_zScored[alldata_Fix_TrialLevel_zScored$age >= 25,], aes(x=Gamma_Event_Number, y=absPositionError)) + stat_smooth(method='lm')+ ggtitle("Gamma Event Number (Fix Only) vs Position Error: Adults Only")
  lm.model <- lm(data = alldata_Fix_TrialLevel_zScored[alldata_Fix_TrialLevel_zScored$age >= 25,], absPositionError ~ Gamma_Event_Number)
  print(anova(lm.model))
  
  gammaTPvsPE <- lmer(absPositionError ~ Gamma_Trial_Power * age + (1|Subject), data = alldata_Fix_TrialLevel_zScored)
  summary(gammaTPvsPE)
  car::Anova(gammaTPvsPE)
  
  
  
}

DelayMinusFixvsBehavior_TrialLevel <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  
  DelayAndFix_TrialLevel <- merge(alldata_Fix_TrialLevel[1:14], alldata_TrialLevel[1:14], by = c("Subject", "Trial"))
 
  DelayMinusFix_TrialLevel <-DelayAndFix_TrialLevel[1:2]
  DelayMinusFix_TrialLevel$Gamma.Event_Number <- DelayAndFix_TrialLevel$Gamma.Gamma_Event_Number - DelayAndFix_TrialLevel$Gamma_Event_Number

  
  DelayMinusFix_TrialLevel <- merge(DelayMinusFix_TrialLevel, Behavior, by = c("Subject", "Trial"), all.x = TRUE, all.y = TRUE)

 DelayMinusFix_TrialLevel_zScored <- DelayMinusFix_TrialLevel %>% group_by(Subject) %>% mutate_if(is.numeric, scale, center = T)
 
 DelayMinusFix_TrialLevel_zScored <- merge(DelayMinusFix_TrialLevel_zScored, agefile, by = "Subject")
  
 gammaNEvsPE <- lmer(absPositionError ~ Gamma.Event_Number * age + (1|Subject), data = DelayMinusFix_TrialLevel_zScored)
 summary(gammaNEvsPE)
 car::Anova(gammaNEvsPE)

 
 ggplot(data = DelayMinusFix_TrialLevel_zScored[], aes(x=Gamma.Event_Number, y=absPositionError)) + stat_smooth(method='lm') + ggtitle("Gamma Event Numbers (Delay-Fix) vs Position Error")
 lm.model <- lm(data = DelayMinusFix_TrialLevel_zScored[], absPositionError ~ Gamma.Event_Number)
 print(anova(lm.model))
 
  
 ggplot(data = DelayMinusFix_TrialLevel_zScored[DelayMinusFix_TrialLevel_zScored$age >= 25,], aes(x=Gamma.Event_Number, y=absPositionError)) + stat_smooth(method='lm')  + ggtitle("Gamma Event Numbers (Delay-Fix) vs Position Error: Adults Only")
 lm.model <- lm(data = DelayMinusFix_TrialLevel_zScored[DelayMinusFix_TrialLevel_zScored$age >= 25,], absPositionError ~ Gamma.Event_Number)
 print(anova(lm.model))
 
 gammaTPvsPE <- lmer(absPositionError ~ Gamma.Gamma_Trial_Power * age + (1|Subject), data = DelayMinusFix_TrialLevel_zScored)
 summary(gammaTPvsPE)
 car::Anova(gammaTPvsPE)
 
 
  
}

DelayvsBehavior_TrialLevel_IndividualChannels <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  
  DelayOnly_IndividualChannels_TrialLevel()
  
  #DLPFC
  Gamma_DLPFC <- GammaDelay_Age_Channel %>% filter(Label == "'F3'")
  Gamma_DLPFC$Subject <- Gamma_DLPFC$idvalues
  
  Gamma_DLPFC_behvaior <- merge(Gamma_DLPFC[], Behavior, by = c("Subject", "Trial"), all.x = TRUE, all.y = TRUE)
  Gamma_DLPFC_behvaior_age <- merge(Gamma_DLPFC_behvaior, agefile, by = "Subject")
  Gamma_DLPFC_behvaior_age$Group <- factor(Gamma_DLPFC_behvaior_age$Group)
  
  
  Gamma_DLPFC_behvaior_zScored <- Gamma_DLPFC_behvaior %>% group_by(Subject)  %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number)) < (sd(Gamma_Event_Number) * 2.5)) %>% filter(abs(absPositionError - mean(absPositionError, na.rm= T)) < (sd(absPositionError, na.rm = T) * 2)) %>% mutate_if(is.numeric, scale, center = T)
  
  Gamma_DLPFC_behvaior_zScored <- merge(Gamma_DLPFC_behvaior_zScored, agefile, by = "Subject")
  
  Gamma_DLPFC_behvaior_zScored$Group <- factor(Gamma_DLPFC_behvaior_zScored$Group)
  
  gammaNEvsPE <- lmer(absPositionError ~ Gamma_Event_Number * age + (1|Subject), data = Gamma_DLPFC_behvaior_zScored)
  summary(gammaNEvsPE)
  car::Anova(gammaNEvsPE)
  
  
  ggplot(data = Gamma_DLPFC_behvaior_zScored[], aes(x=Gamma_Event_Number, y=absPositionError)) + stat_smooth(method='lm') + ggtitle("Gamma Event Number in the DLPFC (Delay Only) vs Position Error")
  lm.model <- lm(data = Gamma_DLPFC_behvaior_zScored[], absPositionError ~ Gamma_Event_Number)
  print(anova(lm.model))
  
  
  ggplot(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$age >= 25,], aes(x=Gamma.Gamma_Event_Number, y=absPositionError)) + stat_smooth(method='lm') + ggtitle("Gamma Event Number (Delay Only) vs Position Error: Adults Only")
  lm.model <- lm(data = alldata_TrialLevel_zScored[alldata_TrialLevel_zScored$age >= 25,], absPositionError ~ Gamma.Gamma_Event_Number)
  print(anova(lm.model))
  
  gammaTPvsPE <- lmer(absPositionError ~ Gamma.Gamma_Trial_Power * age + (1|Subject), data = alldata_TrialLevel_zScored)
  summary(gammaTPvsPE)
  car::Anova(gammaTPvsPE)
}

DelayvsBehavior_TrialLevel_AgeGroups_DLPFC <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  
  DelayOnly_IndividualChannels_TrialLevel()
  
  #DLPFC
  Gamma_DLPFC <- GammaDelay_Age_Channel %>% filter(Label == "'F4'")

  Gamma_DLPFC_behvaior <- merge(Gamma_DLPFC[], Behavior, by = c("Subject", "Trial"))
  Gamma_DLPFC_behvaior$Group <- factor(Gamma_DLPFC_behvaior$Group)
  
  
  Gamma_DLPFC_behvaior_zScored <- Gamma_DLPFC_behvaior %>% group_by(Subject)  %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number)) < (sd(Gamma_Event_Number) * 2.5)) %>% filter(abs(absPositionError - mean(absPositionError, na.rm= T)) < (sd(absPositionError, na.rm = T) * 2)) %>% mutate_if(is.numeric, scale, center = T)
  
  Gamma_DLPFC_behvaior_zScored <- merge(Gamma_DLPFC_behvaior_zScored, agefile, by = "Subject")
  
  Gamma_DLPFC_behvaior_zScored$Group <- factor(Gamma_DLPFC_behvaior_zScored$Group)
  
  gammaNEvsPE <- lmer(absPositionError ~ Gamma_Event_Number * age + (1|Subject), data = Gamma_DLPFC_behvaior_zScored)
  summary(gammaNEvsPE)
  car::Anova(gammaNEvsPE)
  
  
  
  # Look at individual age groups dealy vs behavior relationships in the dlpfc 
  #number of events vs position error
  lunaize(ggplot(data = Gamma_DLPFC_behvaior_zScored[], aes(x=Gamma_Event_Number, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Gamma Event Number (Delay Only) vs Accuracy") + xlab("Number of Gamma Events") + ylab("Accuracy") + scale_fill_brewer(palette = "BuGn"))
  
  lm.model <- lm(data = Gamma_DLPFC_behvaior_zScored[Gamma_DLPFC_behvaior_zScored$Group == 4,], absPositionError ~ Gamma_Event_Number)
  print(anova(lm.model))
  
  #same plot, NOT ZSCORED WITH OUTLIER DETECTION 
  
  clean_data  <- Gamma_DLPFC_behvaior %>% group_by(Subject) %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number)) < (sd(Gamma_Event_Number) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = clean_data[], aes(x=Gamma_Event_Number, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Gamma Event Number (Delay Only) vs Accuracy") + xlab("Number of Gamma Events") + ylab("Accuracy (degs)"))
  
  lm.model <- lmer(data = clean_data[clean_data$Group == 4,], absPositionError ~ Gamma_Event_Number + (1|Subject))
  print(car::Anova(lm.model))
  
  # events vs age
  
  lunaize(ggplot(data = clean_data[], aes(x=age, y=Gamma_Event_Number)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("DLPFC Gamma Event Number (Delay Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events"))
  
  lm.model <- lmer(data = clean_data[], Gamma_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  
  #trial power vs position error 
  lunaize(ggplot(data = Gamma_DLPFC_behvaior_zScored[], aes(x=log_Gamma_Power, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Gamma Event Power (Delay Only) vs Accuracy") + xlab("Power") + ylab("Accuracy") + scale_fill_brewer(palette = "BuGn"))
  
  
  lm.model <- lm(data = Gamma_DLPFC_behvaior_zScored[Gamma_DLPFC_behvaior_zScored$Group == 4,], absPositionError ~ log_Gamma_Power)
  print(anova(lm.model))
  
  #trial power vs PE, NOT ZSCORED WITH OUTLIER DETECTION 
  
  clean_data_Power_PE  <- Gamma_DLPFC_behvaior %>% group_by(Subject) %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power)) < (sd(log_Gamma_Power) * 2.5)) %>% filter(abs(absPositionError - mean(absPositionError, na.rm= T)) < (sd(absPositionError, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = clean_data_Power_PE[], aes(x=log_Gamma_Power, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Gamma Power (Delay Only) vs Accuracy") + xlab("log (Gamma Power)") + ylab("Accuracy (degs)"))
  
  lm.model <- lmer(data = clean_data_Power_PE[clean_data_Power_PE$Group == 4,], absPositionError ~ log_Gamma_Power + (1|Subject))
  print(car::Anova(lm.model))
  
  # power vs age
  
  lunaize(ggplot(data = clean_data_Power_PE[], aes(x=age, y=log_Gamma_Power)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("DLPFC Gamma Power (Delay Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)"))
  
  lm.model <- lmer(data = clean_data_Power_PE[], log_Gamma_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  
  #mgs latency vs number of events
  lunaize(ggplot(data = Gamma_DLPFC_behvaior_zScored[], aes(x=Gamma_Event_Number, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Gamma Event Number (Delay Only) vs mgs Latency") + xlab("Number of Gamma Events") + ylab("Latency"))
  
  lm.model <- lm(data = Gamma_DLPFC_behvaior_zScored[Gamma_DLPFC_behvaior_zScored$Group == 4,], mgsLatency ~ Gamma_Event_Number)
  print(anova(lm.model))
  
  # mgs Latency vs number of events, NOT ZSCORED WITH OUTLIER DETECTION 
  clean_data_number_latency  <- Gamma_DLPFC_behvaior %>% group_by(Subject) %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number)) < (sd(Gamma_Event_Number) * 2.5)) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))
  
  lunaize(ggplot(data = clean_data_number_latency[], aes(x=Gamma_Event_Number, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Number of Gamma Events (Delay Only) vs MGS Latency") + xlab("Number of Gamma Events") + ylab("Latency (secs)"))
  
  lm.model <- lmer(data = clean_data_number_latency[clean_data_number_latency$Group == 1,], mgsLatency ~ Gamma_Event_Number + (1|Subject))
  print(car::Anova(lm.model))
  
  
  lunaize(ggplot(data = clean_data_number_latency[], aes(x=age, y=Gamma_Event_Number)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("DLPFC Gamma Event Number (Delay Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events"))
  
  lm.model <- lmer(data = clean_data_number_latency[], Gamma_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  
  
  #mgs latency vs trial power 
  lunaize(ggplot(data = Gamma_DLPFC_behvaior_zScored[], aes(x=log_Gamma_Power, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Gamma Power (Delay Only) vs mgs Latency") + xlab("Power") + ylab("Latency"))
  
  lm.model <- lm(data = Gamma_DLPFC_behvaior_zScored[Gamma_DLPFC_behvaior_zScored$Group == 1,], mgsLatency ~ log_Gamma_Power)
  print(anova(lm.model))
  
  
  #mgs Latency vs trial power, ZSCORED WITH OUTLIER DETECTION 
  clean_data_Power_latency  <- Gamma_DLPFC_behvaior %>% group_by(Subject) %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power)) < (sd(log_Gamma_Power) * 2.5)) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))
  
  lunaize(ggplot(data = clean_data_Power_latency[], aes(x=log_Gamma_Power, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("DLPFC Gamma Power (Delay Only) vs MGS Latency") + xlab("log(Gamma Power)") + ylab("Latency (secs)"))
  
  lm.model <- lmer(data = clean_data_Power_latency[clean_data_Power_latency$Group == 4,], mgsLatency ~ log_Gamma_Power + (1|Subject))
  print(car::Anova(lm.model))
  
  
  lunaize(ggplot(data = clean_data_Power_latency[], aes(x=age, y=log_Gamma_Power)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("DLPFC Gamma Power (Delay Only) vs Age") + xlab("Age") + ylab("log(Gamma Power)"))
  
  lm.model <- lmer(data = clean_data_Power_latency[], log_Gamma_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  
  
  
  
  # Brendans method to see where the age effect actually shows up
  install.packages("tvem")
  library(tvem)
  
  alldata_TrialLevel <- subset(alldata_TrialLevel, select = -c(LunaID,ScanDate,calR2))
  alldata_TrialLevel$Gamma.Gamma_Event_Number <- as.vector(alldata_TrialLevel$Gamma.Gamma_Event_Number)
  alldata_TrialLevel <- alldata_TrialLevel[complete.cases(alldata_TrialLevel),]
  
  alldata_TrialLevel_zScored <- subset(alldata_TrialLevel_zScored, select = -c(LunaID,ScanDate,calR2))
  alldata_TrialLevel_zScored$Gamma.Gamma_Event_Number <- as.vector(alldata_TrialLevel_zScored$Gamma.Gamma_Event_Number)
  alldata_TrialLevel_zScored <- alldata_TrialLevel_zScored[complete.cases(alldata_TrialLevel_zScored),]
  
  test <- subset(alldata_TrialLevel, select = c(Subject, Gamma.Gamma_Event_Number, age, absPositionError))
  
  
  tvemmodelobj <- tvem(data = test,
                       formula = absPositionError ~ Gamma.Gamma_Event_Number,
                       id = Subject,
                       time = age)
  
  bammodel <- tvemmodelobj$back_end_model
  bammodelsum <- summary(bammodel)
  pvi <- grep("primaryvar",row.names(bammodelsum$s.table))
  bammodelpval <- bammodelsum$s.table[pvi,4]
  
  tvem_grid_fitted_coefficients<-tvemmodelobj$grid_fitted_coefficients$primaryvar
  tvem_grid_fitted_coefficients$time_grid<-tvemmodelobj$time_grid
  
}

DelayvsBehavior_TrialLevel_AgeGroups_PFC <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  
  DelayOnly_IndividualChannels_TrialLevel()
  
  #PFC
  Gamma_PFC <- GammaDelay_Age_Channel[grep("Fp|AF", GammaDelay_Age_Channel$Label), ]
  
  Gamma_PFC_behvaior <- merge(Gamma_PFC[], Behavior, by = c("Subject", "Trial"))
  Gamma_PFC_behvaior_age <- merge(Gamma_PFC_behvaior, agefile, by = "Subject")
  Gamma_PFC_behvaior$Group <- factor(Gamma_PFC_behvaior$Group)
  
  
  Gamma_PFC_behvaior_zScored <- Gamma_PFC_behvaior %>% group_by(Subject)  %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number)) < (sd(Gamma_Event_Number) * 2.5)) %>% filter(abs(absPositionError - mean(absPositionError, na.rm= T)) < (sd(absPositionError, na.rm = T) * 2)) %>% mutate_if(is.numeric, scale, center = T)
  
  Gamma_PFC_behvaior_zScored <- merge(Gamma_PFC_behvaior_zScored, agefile, by = "Subject")
  
  Gamma_PFC_behvaior_zScored$Group <- factor(Gamma_PFC_behvaior_zScored$Group)
  
  
  # Look at individual age groups dealy vs behavior relationships in the PFC 
  #number of events vs position error
  lunaize(ggplot(data = Gamma_PFC_behvaior_zScored[], aes(x=Gamma_Event_Number, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Event Number (Delay Only) vs Accuracy") + xlab("Number of Gamma Events") + ylab("Accuracy") + scale_fill_brewer(palette = "BuGn"))
  
  lm.model <- lmer(data = Gamma_PFC_behvaior_zScored[Gamma_PFC_behvaior_zScored$Group == 4,], absPositionError ~ Gamma_Event_Number + (1|Subject))
  print(anova(lm.model))
  
  #number of events vs position error, NOT Z SCORED WITH OUTLIER DETECTION
  
  clean_data  <- Gamma_PFC_behvaior %>% group_by(Subject) %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number)) < (sd(Gamma_Event_Number) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = clean_data[], aes(x=Gamma_Event_Number, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Event Number (Delay Only) vs Accuracy") + xlab("Number of Gamma Events") + ylab("Accuracy"))
  
  lm.model <- lmer(absPositionError ~ Gamma_Event_Number + (1|Subject), data = clean_data[clean_data$Group == 4,])
  (car::Anova(lm.model))
  
  lunaize(ggplot(data = clean_data[], aes(x=age, y=Gamma_Event_Number)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("PFC Gamma Event Number (Delay Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events"))
  
  lm.model <- lmer(Gamma_Event_Number ~ inverseAge + (1|Subject), data = clean_data[])
  (car::Anova(lm.model))
  
  
  #trial power vs position error, Z SCORED WITHOUT OUTLIER DETECTION 
  lunaize(ggplot(data = Gamma_PFC_behvaior_zScored[], aes(x=log_Gamma_Power, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Event Power (Delay Only) vs Accuracy") + xlab("Power") + ylab("Accuracy") + scale_fill_brewer(palette = "BuGn"))
  
  
  lm.model <- lm(data = Gamma_PFC_behvaior_zScored[Gamma_PFC_behvaior_zScored$Group == 4,], absPositionError ~ log_Gamma_Power)
  print(anova(lm.model))
  
  
  #trial power vs position error, NOT Z SCORED WITH OUTLIER DETECTION 
  
  clean_data_power_PE  <- Gamma_PFC_behvaior %>% group_by(Subject) %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power)) < (sd(log_Gamma_Power) * 2)) %>% filter(abs(PositionError - mean(PositionError, na.rm= T)) < (sd(PositionError, na.rm = T) * 2))
 
   lunaize(ggplot(data = clean_data_power_PE[], aes(x=log_Gamma_Power, y=absPositionError, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Power (Delay Only) vs Accuracy") + xlab("Power") + ylab("Accuracy"))
  
  lm.model <- lmer(data = clean_data_power_PE[clean_data_power_PE$Group == 4,], absPositionError ~ log_Gamma_Power + (1|Subject))
  print(car::Anova(lm.model))
  
  lunaize(ggplot(data = clean_data_power_PE[], aes(x=age, y=log_Gamma_Power)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("PFC Gamma Power (Delay Only) vs Age") + xlab("Age") + ylab("Power"))
  
  lm.model <- lmer(data = clean_data_power_PE[], log_Gamma_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  
  #mgs latency vs number of events
  lunaize(ggplot(data = Gamma_PFC_behvaior_zScored[], aes(x=Gamma_Event_Number, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Event Number (Delay Only) vs mgs Latency") + xlab("Number of Gamma Events") + ylab("Latency"))
  
  lm.model <- lm(data = Gamma_PFC_behvaior_zScored[Gamma_PFC_behvaior_zScored$Group == 4,], mgsLatency ~ Gamma_Event_Number)
  print(anova(lm.model))
  
  #mgs latency vs number of events, NOT Z SCORED WITH OUTLIER DETECTION
  clean_data_Latency  <- Gamma_PFC_behvaior %>% group_by(Subject) %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number)) < (sd(Gamma_Event_Number) * 2)) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = clean_data_Latency[], aes(x=Gamma_Event_Number, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Event Number (Delay Only) vs Latency") + xlab("Number of Gamma Events") + ylab("Latency (s)"))
  
  lm.model <- lmer(data = clean_data_Latency[clean_data_Latency$Group == 3,], mgsLatency ~ Gamma_Event_Number + (1|Subject))
  print(car::Anova(lm.model))

  lunaize(ggplot(data = clean_data_Latency[], aes(x=age, y=Gamma_Event_Number)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("PFC Gamma Event Number (Delay Only) vs Age") + xlab("Age") + ylab("Number of Gamma Events"))
  
  lm.model <- lm(data = clean_data_Latency[], Gamma_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  #mgs latency vs trial power, Z SCORED WITHOUT OUTLIER DETECTION 
  lunaize(ggplot(data = Gamma_PFC_behvaior_zScored[], aes(x=log_Gamma_Power, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Power (Delay Only) vs mgs Latency") + xlab("Power") + ylab("Latency"))
  
  lm.model <- lm(data = Gamma_PFC_behvaior_zScored[Gamma_PFC_behvaior_zScored$Group == 1,], mgsLatency ~ log_Gamma_Power)
  print(anova(lm.model))
  
  #mgs Latency vs Trial power, NOT ZSCORED WITH OUTLIER DETECTION 
  
  clean_data_Latency_Power  <- Gamma_PFC_behvaior %>% group_by(Subject) %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power)) < (sd(log_Gamma_Power) * 2)) %>% filter(abs(mgsLatency - mean(mgsLatency, na.rm= T)) < (sd(mgsLatency, na.rm = T) * 2))
  
  
  lunaize(ggplot(data = clean_data_Latency_Power[], aes(x=log_Gamma_Power, y=mgsLatency, group = Group, color = Group)) + stat_smooth(method='lm', se = FALSE) + ggtitle("PFC Gamma Power (Delay Only) vs Latency") + xlab("Gamma Power") + ylab("Latency (s)"))
  
  lm.model <- lmer(data = clean_data_Latency_Power[clean_data_Latency_Power$Group == 4,], mgsLatency ~ log_Gamma_Power + (1|Subject))
  print(car::Anova(lm.model))
  
  lunaize(ggplot(data = clean_data_Latency_Power[], aes(x=age, y=log_Gamma_Power)) + stat_smooth(method='lm', formula='y~I(1/x)') + ggtitle("PFC Gamma Power (Delay Only) vs Age") + xlab("Age") + ylab("Gamma Power"))
  
  lm.model <- lmer(data = clean_data_Latency_Power[], log_Gamma_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  
  
}

