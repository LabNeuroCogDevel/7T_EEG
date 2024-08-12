
#script to look if the length of the delay period matters

DelayLengths <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  
  DelayGammaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_3_4.csv')
  DelayBetaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_data3_4.csv')
  DelayAlphaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_data_3_4.csv')
  DelayThetaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_data3_4.csv')
  
  keepingEverything <- merge(Behavior, DelayGammaTrials[], by = c("Subject", "Trial"), all.x = TRUE, all.y = TRUE)
  keepingEverything <- merge(keepingEverything, DelayBetaTrials[], by = c("Subject", "Trial"), all.x = TRUE, all.y = TRUE)
  keepingEverything <- merge(keepingEverything, DelayAlphaTrials[], by = c("Subject", "Trial"), all.x = TRUE, all.y = TRUE)
  keepingEverything <- merge(keepingEverything, DelayThetaTrials[], by = c("Subject", "Trial"), all.x = TRUE, all.y = TRUE)
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  
  keepingEverything <- merge(keepingEverything,agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  completeSubjects <-  merge(Behavior, DelayGammaTrials[], by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, DelayBetaTrials[], by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, DelayAlphaTrials[], by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, DelayThetaTrials[], by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, agefile, by = "Subject")
  completeSubjects$absPositionError <- abs(completeSubjects$PositionError)
  completeSubjects$Gamma_log1p_Trial_Power <- log1p(completeSubjects$Gamma_Trial_Power)
  completeSubjects$Beta_log1p_Trial_Power <- log1p(completeSubjects$Beta_Trial_Power)
  completeSubjects$Theta_log1p_Trial_Power <- log1p(completeSubjects$Theta_Trial_Power)
  completeSubjects$Alpha_log1p_Trial_Power <- log1p(completeSubjects$Alpha_Trial_Power)
  
  completeSubjects$Subject <- as.factor(completeSubjects$Subject)
  
  test_df <- na.omit(completeSubjects) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  
  avg_ageGroup <- aggregate(test_df, by = list(test_df$Delay, test_df$Group), mean)
  var_ageGroup <- test_df %>% group_by(Delay,Group) %>% summarise_if(is.numeric,sd)
  
 
   #position error
   ## average 
  
lunaize(ggplot(data = avg_ageGroup[], aes(x = Delay, y = absPositionError, color = Group, group = Group))  + geom_point() + stat_smooth(method = 'lm') + theme(plot.title = element_text(hjust = 0.5)) + xlab("Delay Period Length (s)") + ylab("Avg. Position Error"))
  
  avg_ageGroup %>% split(paste(.$Delay,.$Group)) %>% lapply(lm,formula=absPositionError~Delay) %>% lapply(summary)
  
  pos_err_by_dly <- avg_ageGroup %>% 
    split(.$Group) %>% 
    lapply(lm,formula=absPositionError~Delay)
  slope_i <- 2
  pr_name <- "Pr(>|t|)"
  slope_pvals <-
    pos_err_by_dly %>%
    sapply(function(.) summary(.) %>% coef %>% `[`(slope_i, pr_name))
  
  ped_pvals_avg <- sapply(pos_err_by_dly, function(m) coef(summary(m))["Delay", "Pr(>|t|)"])
  
  ## Variability 
  
  lunaize(ggplot(data = var_ageGroup[], aes(x = Delay, y = absPositionError, color = Group, group = Group))  + geom_point() + stat_smooth(method = 'lm') + theme(plot.title = element_text(hjust = 0.5)) + xlab("Delay Period Length (s)") + ylab("Position Error Variability"))
  
  var_ageGroup %>% split(paste(.$Delay,.$Group)) %>% lapply(lm,formula=absPositionError~Delay) %>% lapply(summary)
  
  pos_err_var_by_dly <- var_ageGroup %>% 
    split(.$Group) %>% 
    lapply(lm,formula=absPositionError~Delay)
  slope_i <- 2
  pr_name <- "Pr(>|t|)"
  slope_pvals <-
    pos_err_var_by_dly %>%
    sapply(function(.) summary(.) %>% coef %>% `[`(slope_i, pr_name))
  
  ped_pvals_var <- sapply(pos_err_var_by_dly, function(m) coef(summary(m))["Delay", "Pr(>|t|)"])
  
  # Latency 
  ## average 
  
  lunaize(ggplot(data = avg_ageGroup[], aes(x = Delay, y = mgsLatency, color = Group, group = Group))  + geom_point() + stat_smooth(method = 'lm') + theme(plot.title = element_text(hjust = 0.5)) + xlab("Delay Period Length (s)") + ylab("Avg. MGS Latency"))
  
  avg_ageGroup %>% split(paste(.$Delay,.$Group)) %>% lapply(lm,formula=mgsLatency~Delay) %>% lapply(summary)
  
  Lat_by_dly <- avg_ageGroup %>% 
    split(.$Group) %>% 
    lapply(lm,formula=mgsLatency~Delay)
  slope_i <- 2
  pr_name <- "Pr(>|t|)"
  slope_pvals <-
    Lat_by_dly %>%
    sapply(function(.) summary(.) %>% coef %>% `[`(slope_i, pr_name))
  
  latd_pvals <- sapply(Lat_by_dly, function(m) coef(summary(m))["Delay", "Pr(>|t|)"])
  
  ## Variability
  
  lunaize(ggplot(data = var_ageGroup[], aes(x = Delay, y = mgsLatency, color = Group, group = Group))  + geom_point() + stat_smooth(method = 'lm') + theme(plot.title = element_text(hjust = 0.5)) + xlab("Delay Period Length (s)") + ylab("Avg. MGS Latency"))
  
  var_ageGroup %>% split(paste(.$Delay,.$Group)) %>% lapply(lm,formula=mgsLatency~Delay) %>% lapply(summary)
  
  Latvar_by_dly <- var_ageGroup %>% 
    split(.$Group) %>% 
    lapply(lm,formula=mgsLatency~Delay)
  slope_i <- 2
  pr_name <- "Pr(>|t|)"
  slope_pvals <-
    Latvar_by_dly %>%
    sapply(function(.) summary(.) %>% coef %>% `[`(slope_i, pr_name))
  
  latd_pvals <- sapply(Latvar_by_dly, function(m) coef(summary(m))["Delay", "Pr(>|t|)"])
  
  
  
  
  
  
  
}