

Behavior_SubLevel_shane <- function() {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200929.csv')
  
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate) #behavioral data contains each trial 
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  # BehaviorTable <- data.table(Behavior)
  # trialsPerSub <- BehaviorTable[, .(rowCount = .N), by = Subject]
  # 
  # 
  # Behavior %>% group_by(Subject) %>% rowSums(Behavior$Trial)

  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv')

  Behavior <- merge(Behavior, agefile, by = "Subject")
  
  # remove trials where the abs value of PE is > 23, remove express saccades
  behaviorClean <- Behavior %>%
    mutate(
      mgsLatency = ifelse(mgsLatency < .1, NA, mgsLatency),
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<2, mgsLatency, NA),
      absPositionError = ifelse(abs(scale(mgsLatency)[,1])<2, absPositionError, NA)
    )
  
  behaviorClean <- subset(behaviorClean, visitno == 1)
  
  
  # Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 
  
  outliers <- function(x) {
    
   (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
  }
  
  behaviorClean_new <- behaviorClean
  
  cols = names(behaviorClean[c(11,13)])
  for (col in cols) {
    
    indx <- outliers(behaviorClean[[col]])
    
    behaviorClean_new[[col]] <- Map(replace, behaviorClean_new[[col]], indx, NA)
    behaviorClean_new[[col]] <- as.numeric(behaviorClean_new[[col]])
    
  }
  
  
  # See how many NA values each trial has and then sum by subject 
  # LatencyNAs <- rowSums(is.na(behaviorClean_new[c(10)])) # only look at the latency and abs position error columns 
  # PE_NAs <- rowSums(is.na(behaviorClean_new[c(13)])) # only look at the latency and abs position error columns 
  # 
  # Subject <- behaviorClean_new$Subject # define a column of all the subject values
  # subNAs <- data.frame(cbind(Subject, LatencyNAs, PE_NAs)) # bind the columns together
  # subNAs$LatencyNAs <- as.numeric(as.character(subNAs$LatencyNAs)) # change values to numeric 
  # subNAs$PE_NAs <- as.numeric(as.character(subNAs$PE_NAs))
  # 
  # SubjectNAs <- aggregate(. ~ Subject, data = subNAs, sum) # sum the NAs of each trial together for all the subjects 
  # SubjectNAs <- merge(SubjectNAs, trialsPerSub, by = "Subject")
  # 
  # SubjectNAs$Lat_percentNA <- (SubjectNAs$LatencyNAs / SubjectNAs$rowCount) * 100
  # SubjectNAs$PE_percentNA <- (SubjectNAs$PE_NAs / SubjectNAs$rowCount) * 100
  # 
  # 
  # SubjectNAs <- merge(SubjectNAs, agefile, by = "Subject")


  avgBehavior <- aggregate(.~Subject, data = behaviorClean_new, mean, na.action = na.omit)
  sdBehavior <- aggregate(.~Subject, data = behaviorClean_new, sd, na.action = na.omit)
  subjectBehavior <- merge(avgBehavior, sdBehavior, by = "Subject", suffix = c("", "_sd"))
  subjectBehavior$inverseAge <- 1/subjectBehavior$age
  
  # Clean outliers on subject level 
  
  subjectBehavior_new <- subjectBehavior
  cols = names(subjectBehavior[c(11,13,27,29)]) # only going to do latency and PE and the variability of them
  for (col in cols) {
    
    indx <- outliers(subjectBehavior[[col]])
    
    subjectBehavior_new[[col]] <- Map(replace, subjectBehavior_new[[col]], indx, NA)
    subjectBehavior_new[[col]] <- as.numeric(subjectBehavior_new[[col]])
    
  }


  return(subjectBehavior_new)
  
  
}

Behavior_Sublevel_Maria <- function() {
  
  library("dplyr")
  # load behavior
  z_thresh = 2
  data <- merge(
    read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200929.csv') %>% mutate(Subject = paste0(LunaID, '_', ScanDate)), 
    read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv'), 
    by='Subject')
  #filter out entire trial if it is not to be trusted
  #MP - added best error and mgs latency filtering
  Behavior <- data %>%
    filter(visitno == 1) %>%
    mutate(absPositionError = abs(PositionError)) %>%
    mutate(absBestError = abs(BestError)) %>%
    filter(absPositionError < 23) %>% 
    filter(absBestError < 23) %>%
    filter(mgsLatency > .1) 
  # Compute per-trial z-scores, and remove abs(z)>2 for absPosErr & mgsLat
  #MP - added best error, fixed position error typo
  behaviorClean <- Behavior %>%
    mutate(
      mgsLatency = ifelse(mgsLatency < .1, NA, mgsLatency), 
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA)
    )
  # Merge with age file
  # MP - added best error
  m_mean <- behaviorClean %>%
    group_by(Subject, age) %>%
    summarize_at(vars(absPositionError, absBestError, mgsLatency, vgsLatency), mean, na.rm=T) %>%  mutate(inverseAge = 1/age)
  
  m_sd <- behaviorClean %>%
    group_by(Subject, age) %>%
    summarize_at(vars(absPositionError, absBestError, mgsLatency, vgsLatency), sd, na.rm=T) %>%  mutate(inverseAge = 1/age)
  
  m <- merge(m_mean, m_sd, by = "Subject", suffixes = c("","_sd"))
  
  # Check for subject-level outliers
  # MP - added best error
  m_clean <- m %>% ungroup() %>%
    dplyr::mutate(
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA),
      mgsLatency_sd = ifelse(abs(scale(mgsLatency_sd)[,1])<z_thresh, mgsLatency_sd, NA),
      absPositionError_sd = ifelse(abs(scale(absPositionError_sd)[,1])<z_thresh, absPositionError_sd, NA), 
      absBestError_sd = ifelse(abs(scale(absBestError_sd)[,1])<z_thresh, absBestError_sd, NA)
    )
  summary(lm(absPositionError ~ inverseAge, data=m_clean))
  ggplot(data=m_clean, aes(x=age, y=absPositionError)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
  
  return(m_clean)
  
  
  
}

Behavior_TrialLevel_shane <- function () {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200929.csv')
  
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate) #behavioral data contains each trial 
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  BehaviorTable <- data.table(Behavior)
  trialsPerSub <- BehaviorTable[, .(rowCount = .N), by = Subject]
  
  
  Behavior %>% group_by(Subject) %>% rowSums(Behavior$Trial)
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  
  # remove trials where the abs value of PE is > 23 
  behaviorClean <- Behavior[Behavior$absPositionError < 23,]
  behaviorClean <- subset(behaviorClean, visitno == 1)
  
  
  # Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
  }
  
  behaviorClean_new <- behaviorClean
  
  cols = names(behaviorClean[c(7:11,13)])
  for (col in cols) {
    
    indx <- outliers(behaviorClean[[col]])
    
    behaviorClean_new[[col]] <- Map(replace, behaviorClean_new[[col]], indx, NA)
    behaviorClean_new[[col]] <- as.numeric(behaviorClean_new[[col]])
    
  }
  
  
return(behaviorClean_new)
  
}

Behavior_TrialLevel_Maria <- function() {
  
  library("dplyr")
  # load behavior
  z_thresh = 2
  data <- merge(
    read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200929.csv') %>% mutate(Subject = paste0(LunaID, '_', ScanDate)), 
    read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv'), 
    by='Subject')
  #filter out entire trial if it is not to be trusted
  #MP - added best error and mgs latency filtering
  Behavior <- data %>%
    filter(visitno == 1) %>%
    mutate(absPositionError = abs(PositionError)) %>%
    mutate(absBestError = abs(BestError)) %>%
    filter(absPositionError < 23) %>% 
    filter(absBestError < 23) %>%
    filter(mgsLatency > .1) 
  # Compute per-trial z-scores, and remove abs(z)>2 for absPosErr & mgsLat
  #MP - added best error, fixed position error typo
  behaviorClean <- Behavior %>%
    mutate(
      mgsLatency = ifelse(mgsLatency < .1, NA, mgsLatency), 
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA)
    )
  
}

Combine_New_Subs <- function () {
  
  # Gamma Frequency Band 
  newsubsGammaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_newSubs_1_2.csv')
  newsubsGammaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_newSubs_2_3.csv')
  newsubsGammaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_newSubs_3_4.csv')
  newsubsGammaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_newSubs_4_5.csv')
  newsubsGammaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_newSubs_5_6.csv')
  
  GammaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_data_1_2.csv')
  GammaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_all_data_2_3.csv')
  GammaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_all_data_3_4.csv')
  GammaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_all_data_4_5.csv')
  GammaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_all_data_5_6.csv')
  
  GammaDelay1_2 <- rbind(GammaDelay1_2, newsubsGammaDelay1_2)
  GammaDelay2_3 <- rbind(GammaDelay2_3, newsubsGammaDelay2_3)
  GammaDelay3_4 <- rbind(GammaDelay3_4, newsubsGammaDelay3_4)
  GammaDelay4_5 <- rbind(GammaDelay4_5, newsubsGammaDelay4_5)
  GammaDelay5_6 <- rbind(GammaDelay5_6, newsubsGammaDelay5_6)
  
  write.csv(GammaDelay1_2, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_data_1_2_20200810.csv')
  write.csv(GammaDelay2_3, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_all_data_2_3_20200810.csv')
  write.csv(GammaDelay3_4, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_all_data_3_4_20200810.csv')
  write.csv(GammaDelay4_5, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_all_data_4_5_20200810.csv')
  write.csv(GammaDelay5_6, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_all_data_5_6_20200810.csv')
  
  
  #alpha frequency band
  newsubsAlphaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_newSubs_data1_2.csv')
  newsubsAlphaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_newSubs_data2_3.csv')
  newsubsAlphaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_newSubs_data3_4.csv')
  newsubsAlphaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_newSubs_data4_5.csv')
  newsubsAlphaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_newSubs_data5_6.csv')
  
  AlphaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_all_data1_2.csv')
  AlphaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_all_data2_3.csv')
  AlphaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_all_data3_4.csv')
  AlphaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_all_data4_5.csv')
  AlphaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_all_data5_6.csv')
  
  AlphaDelay1_2 <- rbind(AlphaDelay1_2, newsubsAlphaDelay1_2)
  AlphaDelay2_3 <- rbind(AlphaDelay2_3, newsubsAlphaDelay2_3)
  AlphaDelay3_4 <- rbind(AlphaDelay3_4, newsubsAlphaDelay3_4)
  AlphaDelay4_5 <- rbind(AlphaDelay4_5, newsubsAlphaDelay4_5)
  AlphaDelay5_6 <- rbind(AlphaDelay5_6, newsubsAlphaDelay5_6)
  
  write.csv(AlphaDelay1_2, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_all_data_1_2_20200810.csv')
  write.csv(AlphaDelay2_3, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_all_data_2_3_20200810.csv')
  write.csv(AlphaDelay3_4, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_all_data_3_4_20200810.csv')
  write.csv(AlphaDelay4_5, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_all_data_4_5_20200810.csv')
  write.csv(AlphaDelay5_6, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_all_data_5_6_20200810.csv')
  
  
  
  #theta frequency band 
  newsubsThetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_newSubs_data1_2.csv')
  newsubsThetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_newSubs_data2_3.csv')
  newsubsThetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_newSubs_data3_4.csv')
  newsubsThetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_newSubs_data4_5.csv')
  newsubsThetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_newSubs_data5_6.csv')
  
  ThetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_all_data1_2.csv')
  ThetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_all_data2_3.csv')
  ThetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_all_data3_4.csv')
  ThetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_all_data4_5.csv')
  ThetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_all_data5_6.csv')
  
  ThetaDelay1_2 <- rbind(ThetaDelay1_2, newsubsThetaDelay1_2)
  ThetaDelay2_3 <- rbind(ThetaDelay2_3, newsubsThetaDelay2_3)
  ThetaDelay3_4 <- rbind(ThetaDelay3_4, newsubsThetaDelay3_4)
  ThetaDelay4_5 <- rbind(ThetaDelay4_5, newsubsThetaDelay4_5)
  ThetaDelay5_6 <- rbind(ThetaDelay5_6, newsubsThetaDelay5_6)
  
  
  write.csv(ThetaDelay1_2, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_all_data_1_2_20200810.csv')
  write.csv(ThetaDelay2_3, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_all_data_2_3_20200810.csv')
  write.csv(ThetaDelay3_4, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_all_data_3_4_20200810.csv')
  write.csv(ThetaDelay4_5, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_all_data_4_5_20200810.csv')
  write.csv(ThetaDelay5_6, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_all_data_5_6_20200810.csv')
  
  
  #beta frequency band
  
  newsubsBetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_newSubs_data1_2.csv')
  newsubsBetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_newSubs_data2_3.csv')
  newsubsBetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_newSubs_data3_4.csv')
  newsubsBetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_newSubs_data4_5.csv')
  newsubsBetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_newSubs_data5_6.csv')
  
  BetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_data1_2.csv')
  BetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_all_data2_3.csv')
  BetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_all_data3_4.csv')
  BetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_all_data4_5.csv')
  BetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_all_data5_6.csv')
  
  BetaDelay1_2 <- rbind(BetaDelay1_2, newsubsBetaDelay1_2)
  BetaDelay2_3 <- rbind(BetaDelay2_3, newsubsBetaDelay2_3)
  BetaDelay3_4 <- rbind(BetaDelay3_4, newsubsBetaDelay3_4)
  BetaDelay4_5 <- rbind(BetaDelay4_5, newsubsBetaDelay4_5)
  BetaDelay5_6 <- rbind(BetaDelay5_6, newsubsBetaDelay5_6)
  
  
  write.csv(BetaDelay1_2, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_data_1_2_20200810.csv')
  write.csv(BetaDelay2_3, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_all_data_2_3_20200810.csv')
  write.csv(BetaDelay3_4, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_all_data_3_4_20200810.csv')
  write.csv(BetaDelay4_5, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_all_data_4_5_20200810.csv')
  write.csv(BetaDelay5_6, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_all_data_5_6_20200810.csv')
  
  
  
  # FIXATION PERIOD
  newGammaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_newSubs_FIX.csv')
  newAlphaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_newSubs_datafix.csv')
  newBetaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_newSUbs_datafix.csv')
  newThetaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_newSubs_datafix.csv')
  
  GammaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_all_data_fix.csv')
  AlphaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_all_datafix.csv')
  BetaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_all_datafix.csv')
  ThetaT  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_all_datafix.csv')
  
  GammaT <- rbind(GammaT, newGammaT)
  AlphaT <- rbind(AlphaT, newAlphaT)
  BetaT <- rbind(BetaT, newBetaT)
  ThetaT <- rbind(ThetaT, newThetaT)
  
  
  write.csv(GammaT, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_all_data_fix_20200810.csv')
  write.csv(AlphaT, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_all_datafix_20200810.csv')
  write.csv(BetaT, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_all_datafix_20200810.csv')
  write.csv(ThetaT, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_all_datafix_20200810.csv')

}

Combine_Delay_Windows <- function(){
  

# Comparing the Delay Periods 
## Gamma Frequency Band 
agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
agefile$Subject <- agefile$idvalues 


GammaDelay1_2$log1p_Gamma_Trial_Power <- log1p(GammaDelay1_2$Gamma_Trial_Power)
GammaDelay2_3$log1p_Gamma_Trial_Power <- log1p(GammaDelay2_3$Gamma_Trial_Power)
GammaDelay3_4$log1p_Gamma_Trial_Power <- log1p(GammaDelay3_4$Gamma_Trial_Power)
GammaDelay4_5$log1p_Gamma_Trial_Power <- log1p(GammaDelay4_5$Gamma_Trial_Power)
GammaDelay5_6$log1p_Gamma_Trial_Power <- log1p(GammaDelay5_6$Gamma_Trial_Power)

Avg_power1_2 <- aggregate(log1p_Gamma_Trial_Power~ Subject, GammaDelay1_2[], mean)
Avg_power2_3 <- aggregate(log1p_Gamma_Trial_Power~ Subject, GammaDelay2_3[], mean)
Avg_power3_4 <- aggregate(log1p_Gamma_Trial_Power~ Subject, GammaDelay3_4[], mean)
Avg_power4_5 <- aggregate(log1p_Gamma_Trial_Power~ Subject, GammaDelay4_5[], mean)
Avg_power5_6 <- aggregate(log1p_Gamma_Trial_Power~ Subject, GammaDelay5_6[], mean)

sd_power1_2 <- aggregate(log1p_Gamma_Trial_Power ~ Subject, GammaDelay1_2[], sd)
sd_power2_3 <- aggregate(log1p_Gamma_Trial_Power ~ Subject, GammaDelay2_3[], sd)
sd_power3_4 <- aggregate(log1p_Gamma_Trial_Power ~ Subject, GammaDelay3_4[], sd)
sd_power4_5 <- aggregate(log1p_Gamma_Trial_Power ~ Subject, GammaDelay4_5[], sd)
sd_power5_6 <- aggregate(log1p_Gamma_Trial_Power ~ Subject, GammaDelay5_6[], sd)

#delay 1-2
Avg_GammaDelay1_2 <- setNames(aggregate(list(GammaDelay1_2$Gamma_Event_Number, GammaDelay1_2$Gamma_Event_Duration), by = list(GammaDelay1_2$Subject),mean), c("Subject", "Gamma_Event_Number","Gamma_Event_Duration"))

SD_GammaDelay1_2 <- setNames(aggregate(list(GammaDelay1_2$Gamma_Event_Number, GammaDelay1_2$Gamma_Event_Duration), by = list(GammaDelay1_2$Subject),sd), c("Subject", "Gamma_Event_Number_SD","Gamma_Event_Duration_SD"))

Avg_GammaDelay1_2 <- merge(Avg_GammaDelay1_2, Avg_power1_2, by = "Subject", suffix = c("", "_mean"))
Avg_GammaDelay1_2 <- merge(Avg_GammaDelay1_2, sd_power1_2, by = "Subject", suffix = c("", "_SD"))
Avg_GammaDelay1_2 <- merge(Avg_GammaDelay1_2, SD_GammaDelay1_2, by = "Subject")


#delay 2-3
Avg_GammaDelay2_3 <- setNames(aggregate(list(GammaDelay2_3$Gamma_Event_Number, GammaDelay2_3$Gamma_Event_Duration), by = list(GammaDelay2_3$Subject),mean), c("Subject", "Gamma_Event_Number","Gamma_Event_Duration"))
Avg_GammaDelay2_3 <- merge(Avg_GammaDelay2_3, Avg_power2_3, by = "Subject")
Avg_GammaDelay2_3 <- merge(Avg_GammaDelay2_3, sd_power2_3, by = "Subject",suffix = c("", "_SD"))

SD_GammaDelay2_3 <- setNames(aggregate(list(GammaDelay2_3$Gamma_Event_Number, GammaDelay2_3$Gamma_Event_Duration), by = list(GammaDelay2_3$Subject),sd), c("Subject", "Gamma_Event_Number_SD","Gamma_Event_Duration_SD"))
Avg_GammaDelay2_3 <- merge(Avg_GammaDelay2_3, SD_GammaDelay2_3, by = "Subject")


#delay 3-4
Avg_GammaDelay3_4 <- setNames(aggregate(list(GammaDelay3_4$Gamma_Event_Number, GammaDelay3_4$Gamma_Event_Duration), by = list(GammaDelay3_4$Subject),mean), c("Subject", "Gamma_Event_Number","Gamma_Event_Duration"))
Avg_GammaDelay3_4 <- merge(Avg_GammaDelay3_4, Avg_power3_4, by = "Subject")
Avg_GammaDelay3_4 <- merge(Avg_GammaDelay3_4, sd_power3_4, by = "Subject", suffix = c("", "_SD"))

SD_GammaDelay3_4 <- setNames(aggregate(list(GammaDelay3_4$Gamma_Event_Number, GammaDelay3_4$Gamma_Event_Duration), by = list(GammaDelay3_4$Subject),sd), c("Subject", "Gamma_Event_Number_SD","Gamma_Event_Duration_SD"))
Avg_GammaDelay3_4 <- merge(Avg_GammaDelay3_4, SD_GammaDelay3_4, by = "Subject")


#delay 4-5
Avg_GammaDelay4_5 <- setNames(aggregate(list(GammaDelay4_5$Gamma_Event_Number, GammaDelay4_5$Gamma_Event_Duration), by = list(GammaDelay4_5$Subject),mean), c("Subject",  "Gamma_Event_Number","Gamma_Event_Duration"))
Avg_GammaDelay4_5 <- merge(Avg_GammaDelay4_5, Avg_power4_5, by = "Subject")
Avg_GammaDelay4_5 <- merge(Avg_GammaDelay4_5, sd_power4_5, by = "Subject", suffix = c("", "_SD"))


SD_GammaDelay4_5 <- setNames(aggregate(list(GammaDelay4_5$Gamma_Event_Number, GammaDelay4_5$Gamma_Event_Duration), by = list(GammaDelay4_5$Subject),sd), c("Subject", "Gamma_Event_Number_SD","Gamma_Event_Duration_SD"))
Avg_GammaDelay4_5 <- merge(Avg_GammaDelay4_5, SD_GammaDelay4_5, by = "Subject")


#delay 5-6
Avg_GammaDelay5_6 <- setNames(aggregate(list( GammaDelay5_6$Gamma_Event_Number, GammaDelay5_6$Gamma_Event_Duration), by = list(GammaDelay5_6$Subject),mean), c("Subject", "Gamma_Event_Number","Gamma_Event_Duration"))
Avg_GammaDelay5_6 <- merge(Avg_GammaDelay5_6, Avg_power5_6, by = "Subject")
Avg_GammaDelay5_6 <- merge(Avg_GammaDelay5_6, sd_power5_6, by = "Subject", suffix = c("", "_SD"))

SD_GammaDelay5_6 <- setNames(aggregate(list(GammaDelay5_6$Gamma_Event_Number, GammaDelay5_6$Gamma_Event_Duration), by = list(GammaDelay5_6$Subject),sd), c("Subject", "Gamma_Event_Number_SD","Gamma_Event_Duration_SD"))
Avg_GammaDelay5_6 <- merge(Avg_GammaDelay5_6, SD_GammaDelay5_6, by = "Subject")



Avg_GammaDelay1_2$epoch <- "1-2"
Avg_GammaDelay2_3$epoch <- "2-3"
Avg_GammaDelay3_4$epoch <- "3-4"
Avg_GammaDelay4_5$epoch <- "4-5"
Avg_GammaDelay5_6$epoch <- "5-6"

allDelayGamma <- rbind(Avg_GammaDelay1_2, Avg_GammaDelay2_3)
allDelayGamma <- rbind(allDelayGamma, Avg_GammaDelay3_4)
allDelayGamma <- rbind(allDelayGamma, Avg_GammaDelay4_5)
allDelayGamma <- rbind(allDelayGamma, Avg_GammaDelay5_6)
allDelayGamma <- merge(allDelayGamma, agefile, by = "Subject")

delayGamma <- setNames(aggregate(list(allDelayGamma$log1p_Gamma_Trial_Power,allDelayGamma$log1p_Gamma_Trial_Power_SD, allDelayGamma$Gamma_Event_Number, allDelayGamma$Gamma_Event_Number_SD, allDelayGamma$Gamma_Event_Duration, allDelayGamma$Gamma_Event_Duration_SD), by = list(allDelayGamma$Subject),mean), c("Subject", "Trial_Power","Trial_Power_Variability", "Event_Number","Event_Number_Variability", "Event_Duration", "Event_Duration_Variability"))

write.csv(delayGamma, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/delayGamma.csv')



## Alpha Frequency Band 

AlphaDelay1_2$log1p_Alpha_Trial_Power <- log1p(AlphaDelay1_2$Alpha_Trial_Power)
AlphaDelay2_3$log1p_Alpha_Trial_Power <- log1p(AlphaDelay2_3$Alpha_Trial_Power)
AlphaDelay3_4$log1p_Alpha_Trial_Power <- log1p(AlphaDelay3_4$Alpha_Trial_Power)
AlphaDelay4_5$log1p_Alpha_Trial_Power <- log1p(AlphaDelay4_5$Alpha_Trial_Power)
AlphaDelay5_6$log1p_Alpha_Trial_Power <- log1p(AlphaDelay5_6$Alpha_Trial_Power)

Avg_power1_2 <- aggregate(log1p_Alpha_Trial_Power~ Subject, AlphaDelay1_2[], mean)
Avg_power2_3 <- aggregate(log1p_Alpha_Trial_Power~ Subject, AlphaDelay2_3[], mean)
Avg_power3_4 <- aggregate(log1p_Alpha_Trial_Power~ Subject, AlphaDelay3_4[], mean)
Avg_power4_5 <- aggregate(log1p_Alpha_Trial_Power~ Subject, AlphaDelay4_5[], mean)
Avg_power5_6 <- aggregate(log1p_Alpha_Trial_Power~ Subject, AlphaDelay5_6[], mean)

sd_power1_2 <- aggregate(log1p_Alpha_Trial_Power ~ Subject, AlphaDelay1_2[], sd)
sd_power2_3 <- aggregate(log1p_Alpha_Trial_Power ~ Subject, AlphaDelay2_3[], sd)
sd_power3_4 <- aggregate(log1p_Alpha_Trial_Power ~ Subject, AlphaDelay3_4[], sd)
sd_power4_5 <- aggregate(log1p_Alpha_Trial_Power ~ Subject, AlphaDelay4_5[], sd)
sd_power5_6 <- aggregate(log1p_Alpha_Trial_Power ~ Subject, AlphaDelay5_6[], sd)

#delay 1-2
Avg_AlphaDelay1_2 <- setNames(aggregate(list(AlphaDelay1_2$Alpha_Event_Number, AlphaDelay1_2$Alpha_Event_Duration), by = list(AlphaDelay1_2$Subject),mean), c("Subject", "Alpha_Event_Number","Alpha_Event_Duration"))

SD_AlphaDelay1_2 <- setNames(aggregate(list(AlphaDelay1_2$Alpha_Event_Number, AlphaDelay1_2$Alpha_Event_Duration), by = list(AlphaDelay1_2$Subject),sd), c("Subject", "Alpha_Event_Number_SD","Alpha_Event_Duration_SD"))

Avg_AlphaDelay1_2 <- merge(Avg_AlphaDelay1_2, Avg_power1_2, by = "Subject", suffix = c("", "_mean"))
Avg_AlphaDelay1_2 <- merge(Avg_AlphaDelay1_2, sd_power1_2, by = "Subject", suffix = c("", "_SD"))
Avg_AlphaDelay1_2 <- merge(Avg_AlphaDelay1_2, SD_AlphaDelay1_2, by = "Subject")


#delay 2-3
Avg_AlphaDelay2_3 <- setNames(aggregate(list(AlphaDelay2_3$Alpha_Event_Number, AlphaDelay2_3$Alpha_Event_Duration), by = list(AlphaDelay2_3$Subject),mean), c("Subject", "Alpha_Event_Number","Alpha_Event_Duration"))
Avg_AlphaDelay2_3 <- merge(Avg_AlphaDelay2_3, Avg_power2_3, by = "Subject")
Avg_AlphaDelay2_3 <- merge(Avg_AlphaDelay2_3, sd_power2_3, by = "Subject",suffix = c("", "_SD"))

SD_AlphaDelay2_3 <- setNames(aggregate(list(AlphaDelay2_3$Alpha_Event_Number, AlphaDelay2_3$Alpha_Event_Duration), by = list(AlphaDelay2_3$Subject),sd), c("Subject", "Alpha_Event_Number_SD","Alpha_Event_Duration_SD"))
Avg_AlphaDelay2_3 <- merge(Avg_AlphaDelay2_3, SD_AlphaDelay2_3, by = "Subject")


#delay 3-4
Avg_AlphaDelay3_4 <- setNames(aggregate(list(AlphaDelay3_4$Alpha_Event_Number, AlphaDelay3_4$Alpha_Event_Duration), by = list(AlphaDelay3_4$Subject),mean), c("Subject", "Alpha_Event_Number","Alpha_Event_Duration"))
Avg_AlphaDelay3_4 <- merge(Avg_AlphaDelay3_4, Avg_power3_4, by = "Subject")
Avg_AlphaDelay3_4 <- merge(Avg_AlphaDelay3_4, sd_power3_4, by = "Subject", suffix = c("", "_SD"))

SD_AlphaDelay3_4 <- setNames(aggregate(list(AlphaDelay3_4$Alpha_Event_Number, AlphaDelay3_4$Alpha_Event_Duration), by = list(AlphaDelay3_4$Subject),sd), c("Subject", "Alpha_Event_Number_SD","Alpha_Event_Duration_SD"))
Avg_AlphaDelay3_4 <- merge(Avg_AlphaDelay3_4, SD_AlphaDelay3_4, by = "Subject")


#delay 4-5
Avg_AlphaDelay4_5 <- setNames(aggregate(list(AlphaDelay4_5$Alpha_Event_Number, AlphaDelay4_5$Alpha_Event_Duration), by = list(AlphaDelay4_5$Subject),mean), c("Subject",  "Alpha_Event_Number","Alpha_Event_Duration"))
Avg_AlphaDelay4_5 <- merge(Avg_AlphaDelay4_5, Avg_power4_5, by = "Subject")
Avg_AlphaDelay4_5 <- merge(Avg_AlphaDelay4_5, sd_power4_5, by = "Subject", suffix = c("", "_SD"))


SD_AlphaDelay4_5 <- setNames(aggregate(list(AlphaDelay4_5$Alpha_Event_Number, AlphaDelay4_5$Alpha_Event_Duration), by = list(AlphaDelay4_5$Subject),sd), c("Subject", "Alpha_Event_Number_SD","Alpha_Event_Duration_SD"))
Avg_AlphaDelay4_5 <- merge(Avg_AlphaDelay4_5, SD_AlphaDelay4_5, by = "Subject")


#delay 5-6
Avg_AlphaDelay5_6 <- setNames(aggregate(list( AlphaDelay5_6$Alpha_Event_Number, AlphaDelay5_6$Alpha_Event_Duration), by = list(AlphaDelay5_6$Subject),mean), c("Subject", "Alpha_Event_Number","Alpha_Event_Duration"))
Avg_AlphaDelay5_6 <- merge(Avg_AlphaDelay5_6, Avg_power5_6, by = "Subject")
Avg_AlphaDelay5_6 <- merge(Avg_AlphaDelay5_6, sd_power5_6, by = "Subject", suffix = c("", "_SD"))

SD_AlphaDelay5_6 <- setNames(aggregate(list(AlphaDelay5_6$Alpha_Event_Number, AlphaDelay5_6$Alpha_Event_Duration), by = list(AlphaDelay5_6$Subject),sd), c("Subject", "Alpha_Event_Number_SD","Alpha_Event_Duration_SD"))
Avg_AlphaDelay5_6 <- merge(Avg_AlphaDelay5_6, SD_AlphaDelay5_6, by = "Subject")

Avg_AlphaDelay1_2$epoch <- "1-2"
Avg_AlphaDelay2_3$epoch <- "2-3"
Avg_AlphaDelay3_4$epoch <- "3-4"
Avg_AlphaDelay4_5$epoch <- "4-5"
Avg_AlphaDelay5_6$epoch <- "5-6"

allDelayAlpha <- rbind(Avg_AlphaDelay1_2, Avg_AlphaDelay2_3)
allDelayAlpha <- rbind(allDelayAlpha, Avg_AlphaDelay3_4)
allDelayAlpha <- rbind(allDelayAlpha, Avg_AlphaDelay4_5)
allDelayAlpha <- rbind(allDelayAlpha, Avg_AlphaDelay5_6)
allDelayAlpha <- merge(allDelayAlpha, agefile, by = "Subject")

delayAlpha <- setNames(aggregate(list(allDelayAlpha$log1p_Alpha_Trial_Power,allDelayAlpha$log1p_Alpha_Trial_Power_SD, allDelayAlpha$Alpha_Event_Number, allDelayAlpha$Alpha_Event_Number_SD, allDelayAlpha$Alpha_Event_Duration, allDelayAlpha$Alpha_Event_Duration_SD), by = list(allDelayAlpha$Subject),mean), c("Subject", "Trial_Power","Trial_Power_Variability", "Event_Number","Event_Number_Variability", "Event_Duration", "Event_Duration_Variability"))
write.csv(delayAlpha, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/delayAlpha.csv')


## Theta Frequency Band 

ThetaDelay1_2$log1p_Theta_Trial_Power <- log1p(ThetaDelay1_2$Theta_Trial_Power)
ThetaDelay2_3$log1p_Theta_Trial_Power <- log1p(ThetaDelay2_3$Theta_Trial_Power)
ThetaDelay3_4$log1p_Theta_Trial_Power <- log1p(ThetaDelay3_4$Theta_Trial_Power)
ThetaDelay4_5$log1p_Theta_Trial_Power <- log1p(ThetaDelay4_5$Theta_Trial_Power)
ThetaDelay5_6$log1p_Theta_Trial_Power <- log1p(ThetaDelay5_6$Theta_Trial_Power)

Avg_power1_2 <- aggregate(log1p_Theta_Trial_Power~ Subject, ThetaDelay1_2[], mean)
Avg_power2_3 <- aggregate(log1p_Theta_Trial_Power~ Subject, ThetaDelay2_3[], mean)
Avg_power3_4 <- aggregate(log1p_Theta_Trial_Power~ Subject, ThetaDelay3_4[], mean)
Avg_power4_5 <- aggregate(log1p_Theta_Trial_Power~ Subject, ThetaDelay4_5[], mean)
Avg_power5_6 <- aggregate(log1p_Theta_Trial_Power~ Subject, ThetaDelay5_6[], mean)


sd_power1_2 <- aggregate(log1p_Theta_Trial_Power ~ Subject, ThetaDelay1_2[], sd)
sd_power2_3 <- aggregate(log1p_Theta_Trial_Power ~ Subject, ThetaDelay2_3[], sd)
sd_power3_4 <- aggregate(log1p_Theta_Trial_Power ~ Subject, ThetaDelay3_4[], sd)
sd_power4_5 <- aggregate(log1p_Theta_Trial_Power ~ Subject, ThetaDelay4_5[], sd)
sd_power5_6 <- aggregate(log1p_Theta_Trial_Power ~ Subject, ThetaDelay5_6[,], sd)

#delay 1-2
Avg_ThetaDelay1_2 <- setNames(aggregate(list(ThetaDelay1_2$Theta_Event_Number, ThetaDelay1_2$Theta_Event_Duration), by = list(ThetaDelay1_2$Subject),mean), c("Subject", "Theta_Event_Number","Theta_Event_Duration"))

SD_ThetaDelay1_2 <- setNames(aggregate(list(ThetaDelay1_2$Theta_Event_Number, ThetaDelay1_2$Theta_Event_Duration), by = list(ThetaDelay1_2$Subject),sd), c("Subject", "Theta_Event_Number_SD","Theta_Event_Duration_SD"))

Avg_ThetaDelay1_2 <- merge(Avg_ThetaDelay1_2, Avg_power1_2, by = "Subject", suffix = c("", "_mean"))
Avg_ThetaDelay1_2 <- merge(Avg_ThetaDelay1_2, sd_power1_2, by = "Subject", suffix = c("", "_SD"))
Avg_ThetaDelay1_2 <- merge(Avg_ThetaDelay1_2, SD_ThetaDelay1_2, by = "Subject")


#delay 2-3
Avg_ThetaDelay2_3 <- setNames(aggregate(list(ThetaDelay2_3$Theta_Event_Number, ThetaDelay2_3$Theta_Event_Duration), by = list(ThetaDelay2_3$Subject),mean), c("Subject", "Theta_Event_Number","Theta_Event_Duration"))
Avg_ThetaDelay2_3 <- merge(Avg_ThetaDelay2_3, Avg_power2_3, by = "Subject")
Avg_ThetaDelay2_3 <- merge(Avg_ThetaDelay2_3, sd_power2_3, by = "Subject",suffix = c("", "_SD"))

SD_ThetaDelay2_3 <- setNames(aggregate(list(ThetaDelay2_3$Theta_Event_Number, ThetaDelay2_3$Theta_Event_Duration), by = list(ThetaDelay2_3$Subject),sd), c("Subject", "Theta_Event_Number_SD","Theta_Event_Duration_SD"))
Avg_ThetaDelay2_3 <- merge(Avg_ThetaDelay2_3, SD_ThetaDelay2_3, by = "Subject")


#delay 3-4
Avg_ThetaDelay3_4 <- setNames(aggregate(list(ThetaDelay3_4$Theta_Event_Number, ThetaDelay3_4$Theta_Event_Duration), by = list(ThetaDelay3_4$Subject),mean), c("Subject", "Theta_Event_Number","Theta_Event_Duration"))
Avg_ThetaDelay3_4 <- merge(Avg_ThetaDelay3_4, Avg_power3_4, by = "Subject")
Avg_ThetaDelay3_4 <- merge(Avg_ThetaDelay3_4, sd_power3_4, by = "Subject", suffix = c("", "_SD"))

SD_ThetaDelay3_4 <- setNames(aggregate(list(ThetaDelay3_4$Theta_Event_Number, ThetaDelay3_4$Theta_Event_Duration), by = list(ThetaDelay3_4$Subject),sd), c("Subject", "Theta_Event_Number_SD","Theta_Event_Duration_SD"))
Avg_ThetaDelay3_4 <- merge(Avg_ThetaDelay3_4, SD_ThetaDelay3_4, by = "Subject")


#delay 4-5
Avg_ThetaDelay4_5 <- setNames(aggregate(list(ThetaDelay4_5$Theta_Event_Number, ThetaDelay4_5$Theta_Event_Duration), by = list(ThetaDelay4_5$Subject),mean), c("Subject",  "Theta_Event_Number","Theta_Event_Duration"))
Avg_ThetaDelay4_5 <- merge(Avg_ThetaDelay4_5, Avg_power4_5, by = "Subject")
Avg_ThetaDelay4_5 <- merge(Avg_ThetaDelay4_5, sd_power4_5, by = "Subject", suffix = c("", "_SD"))


SD_ThetaDelay4_5 <- setNames(aggregate(list(ThetaDelay4_5$Theta_Event_Number, ThetaDelay4_5$Theta_Event_Duration), by = list(ThetaDelay4_5$Subject),sd), c("Subject", "Theta_Event_Number_SD","Theta_Event_Duration_SD"))
Avg_ThetaDelay4_5 <- merge(Avg_ThetaDelay4_5, SD_ThetaDelay4_5, by = "Subject")


#delay 5-6
Avg_ThetaDelay5_6 <- setNames(aggregate(list( ThetaDelay5_6$Theta_Event_Number, ThetaDelay5_6$Theta_Event_Duration), by = list(ThetaDelay5_6$Subject),mean), c("Subject", "Theta_Event_Number","Theta_Event_Duration"))
Avg_ThetaDelay5_6 <- merge(Avg_ThetaDelay5_6, Avg_power5_6, by = "Subject")
Avg_ThetaDelay5_6 <- merge(Avg_ThetaDelay5_6, sd_power5_6, by = "Subject", suffix = c("", "_SD"))

SD_ThetaDelay5_6 <- setNames(aggregate(list(ThetaDelay5_6$Theta_Event_Number, ThetaDelay5_6$Theta_Event_Duration), by = list(ThetaDelay5_6$Subject),sd), c("Subject", "Theta_Event_Number_SD","Theta_Event_Duration_SD"))
Avg_ThetaDelay5_6 <- merge(Avg_ThetaDelay5_6, SD_ThetaDelay5_6, by = "Subject")


Avg_ThetaDelay1_2$epoch <- "1-2"
Avg_ThetaDelay2_3$epoch <- "2-3"
Avg_ThetaDelay3_4$epoch <- "3-4"
Avg_ThetaDelay4_5$epoch <- "4-5"
Avg_ThetaDelay5_6$epoch <- "5-6"

allDelayTheta <- rbind(Avg_ThetaDelay1_2, Avg_ThetaDelay2_3)
allDelayTheta <- rbind(allDelayTheta, Avg_ThetaDelay3_4)
allDelayTheta <- rbind(allDelayTheta, Avg_ThetaDelay4_5)
allDelayTheta <- rbind(allDelayTheta, Avg_ThetaDelay5_6)
allDelayTheta <- merge(allDelayTheta, agefile, by = "Subject")

delayTheta <- setNames(aggregate(list(allDelayTheta$log1p_Theta_Trial_Power,allDelayTheta$log1p_Theta_Trial_Power_SD, allDelayTheta$Theta_Event_Number, allDelayTheta$Theta_Event_Number_SD, allDelayTheta$Theta_Event_Duration, allDelayTheta$Theta_Event_Duration_SD), by = list(allDelayTheta$Subject),mean), c("Subject", "Trial_Power","Trial_Power_Variability", "Event_Number","Event_Number_Variability", "Event_Duration", "Event_Duration_Variability"))
write.csv(delayTheta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/delayTheta.csv')


## Beta Frequency Band 
BetaDelay1_2$log1p_Beta_Trial_Power <- log1p(BetaDelay1_2$Beta_Trial_Power)
BetaDelay2_3$log1p_Beta_Trial_Power <- log1p(BetaDelay2_3$Beta_Trial_Power)
BetaDelay3_4$log1p_Beta_Trial_Power <- log1p(BetaDelay3_4$Beta_Trial_Power)
BetaDelay4_5$log1p_Beta_Trial_Power <- log1p(BetaDelay4_5$Beta_Trial_Power)
BetaDelay5_6$log1p_Beta_Trial_Power <- log1p(BetaDelay5_6$Beta_Trial_Power)

Avg_power1_2 <- aggregate(log1p_Beta_Trial_Power~ Subject, BetaDelay1_2[], mean)
Avg_power2_3 <- aggregate(log1p_Beta_Trial_Power~ Subject, BetaDelay2_3[], mean)
Avg_power3_4 <- aggregate(log1p_Beta_Trial_Power~ Subject, BetaDelay3_4[], mean)
Avg_power4_5 <- aggregate(log1p_Beta_Trial_Power~ Subject, BetaDelay4_5[], mean)
Avg_power5_6 <- aggregate(log1p_Beta_Trial_Power~ Subject, BetaDelay5_6[], mean)


sd_power1_2 <- aggregate(log1p_Beta_Trial_Power ~ Subject, BetaDelay1_2[], sd)
sd_power2_3 <- aggregate(log1p_Beta_Trial_Power ~ Subject, BetaDelay2_3[], sd)
sd_power3_4 <- aggregate(log1p_Beta_Trial_Power ~ Subject, BetaDelay3_4[], sd)
sd_power4_5 <- aggregate(log1p_Beta_Trial_Power ~ Subject, BetaDelay4_5[], sd)
sd_power5_6 <- aggregate(log1p_Beta_Trial_Power ~ Subject, BetaDelay5_6[], sd)

#delay 1-2
Avg_BetaDelay1_2 <- setNames(aggregate(list(BetaDelay1_2$Beta_Event_Number, BetaDelay1_2$Beta_Event_Duration), by = list(BetaDelay1_2$Subject),mean), c("Subject", "Beta_Event_Number","Beta_Event_Duration"))

SD_BetaDelay1_2 <- setNames(aggregate(list(BetaDelay1_2$Beta_Event_Number, BetaDelay1_2$Beta_Event_Duration), by = list(BetaDelay1_2$Subject),sd), c("Subject", "Beta_Event_Number_SD","Beta_Event_Duration_SD"))

Avg_BetaDelay1_2 <- merge(Avg_BetaDelay1_2, Avg_power1_2, by = "Subject", suffix = c("", "_mean"))
Avg_BetaDelay1_2 <- merge(Avg_BetaDelay1_2, sd_power1_2, by = "Subject", suffix = c("", "_SD"))
Avg_BetaDelay1_2 <- merge(Avg_BetaDelay1_2, SD_BetaDelay1_2, by = "Subject")


#delay 2-3
Avg_BetaDelay2_3 <- setNames(aggregate(list(BetaDelay2_3$Beta_Event_Number, BetaDelay2_3$Beta_Event_Duration), by = list(BetaDelay2_3$Subject),mean), c("Subject", "Beta_Event_Number","Beta_Event_Duration"))
Avg_BetaDelay2_3 <- merge(Avg_BetaDelay2_3, Avg_power2_3, by = "Subject")
Avg_BetaDelay2_3 <- merge(Avg_BetaDelay2_3, sd_power2_3, by = "Subject",suffix = c("", "_SD"))

SD_BetaDelay2_3 <- setNames(aggregate(list(BetaDelay2_3$Beta_Event_Number, BetaDelay2_3$Beta_Event_Duration), by = list(BetaDelay2_3$Subject),sd), c("Subject", "Beta_Event_Number_SD","Beta_Event_Duration_SD"))
Avg_BetaDelay2_3 <- merge(Avg_BetaDelay2_3, SD_BetaDelay2_3, by = "Subject")


#delay 3-4
Avg_BetaDelay3_4 <- setNames(aggregate(list(BetaDelay3_4$Beta_Event_Number, BetaDelay3_4$Beta_Event_Duration), by = list(BetaDelay3_4$Subject),mean), c("Subject", "Beta_Event_Number","Beta_Event_Duration"))
Avg_BetaDelay3_4 <- merge(Avg_BetaDelay3_4, Avg_power3_4, by = "Subject")
Avg_BetaDelay3_4 <- merge(Avg_BetaDelay3_4, sd_power3_4, by = "Subject", suffix = c("", "_SD"))

SD_BetaDelay3_4 <- setNames(aggregate(list(BetaDelay3_4$Beta_Event_Number, BetaDelay3_4$Beta_Event_Duration), by = list(BetaDelay3_4$Subject),sd), c("Subject", "Beta_Event_Number_SD","Beta_Event_Duration_SD"))
Avg_BetaDelay3_4 <- merge(Avg_BetaDelay3_4, SD_BetaDelay3_4, by = "Subject")


#delay 4-5
Avg_BetaDelay4_5 <- setNames(aggregate(list(BetaDelay4_5$Beta_Event_Number, BetaDelay4_5$Beta_Event_Duration), by = list(BetaDelay4_5$Subject),mean), c("Subject",  "Beta_Event_Number","Beta_Event_Duration"))
Avg_BetaDelay4_5 <- merge(Avg_BetaDelay4_5, Avg_power4_5, by = "Subject")
Avg_BetaDelay4_5 <- merge(Avg_BetaDelay4_5, sd_power4_5, by = "Subject", suffix = c("", "_SD"))


SD_BetaDelay4_5 <- setNames(aggregate(list(BetaDelay4_5$Beta_Event_Number, BetaDelay4_5$Beta_Event_Duration), by = list(BetaDelay4_5$Subject),sd), c("Subject", "Beta_Event_Number_SD","Beta_Event_Duration_SD"))
Avg_BetaDelay4_5 <- merge(Avg_BetaDelay4_5, SD_BetaDelay4_5, by = "Subject")


#delay 5-6
Avg_BetaDelay5_6 <- setNames(aggregate(list( BetaDelay5_6$Beta_Event_Number, BetaDelay5_6$Beta_Event_Duration), by = list(BetaDelay5_6$Subject),mean), c("Subject", "Beta_Event_Number","Beta_Event_Duration"))
Avg_BetaDelay5_6 <- merge(Avg_BetaDelay5_6, Avg_power5_6, by = "Subject")
Avg_BetaDelay5_6 <- merge(Avg_BetaDelay5_6, sd_power5_6, by = "Subject", suffix = c("", "_SD"))

SD_BetaDelay5_6 <- setNames(aggregate(list(BetaDelay5_6$Beta_Event_Number, BetaDelay5_6$Beta_Event_Duration), by = list(BetaDelay5_6$Subject),sd), c("Subject", "Beta_Event_Number_SD","Beta_Event_Duration_SD"))
Avg_BetaDelay5_6 <- merge(Avg_BetaDelay5_6, SD_BetaDelay5_6, by = "Subject")

Avg_BetaDelay1_2$epoch <- "1-2"
Avg_BetaDelay2_3$epoch <- "2-3"
Avg_BetaDelay3_4$epoch <- "3-4"
Avg_BetaDelay4_5$epoch <- "4-5"
Avg_BetaDelay5_6$epoch <- "5-6"

allDelayBeta <- rbind(Avg_BetaDelay1_2, Avg_BetaDelay2_3)
allDelayBeta <- rbind(allDelayBeta, Avg_BetaDelay3_4)
allDelayBeta <- rbind(allDelayBeta, Avg_BetaDelay4_5)
allDelayBeta <- rbind(allDelayBeta, Avg_BetaDelay5_6)
allDelayBeta <- merge(allDelayBeta, agefile, by = "Subject")

delayBeta <- setNames(aggregate(list(allDelayBeta$log1p_Beta_Trial_Power,allDelayBeta$log1p_Beta_Trial_Power_SD, allDelayBeta$Beta_Event_Number, allDelayBeta$Beta_Event_Number_SD, allDelayBeta$Beta_Event_Duration, allDelayBeta$Beta_Event_Duration_SD), by = list(allDelayBeta$Subject),mean), c("Subject", "Trial_Power","Trial_Power_Variability", "Event_Number","Event_Number_Variability", "Event_Duration", "Event_Duration_Variability"))
write.csv(delayBeta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/delayBeta.csv')



# All Delay Period Epochs Averaged Together 

## All Frequency Bands Across Age 

delayGamma$Band <- "Gamma"
delayBeta$Band <- "Beta"
delayTheta$Band <- "Theta"
delayAlpha$Band <- "Alpha"


allBands <- rbind(delayGamma, delayBeta)
allBands <- rbind(allBands, delayTheta)
allBands <- rbind(allBands, delayAlpha)

allBands <- merge(allBands, agefile, by = "Subject")
write.csv(allBands, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/allBands_DelayPeriod.csv')


# Peak Power and Frequency 
## Delay Period
### Gamma Frequency Band 

GammaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_PeakFreq_Power1_2.csv')
GammaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_all_PeakFreq_Power2_3.csv')
GammaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_all_PeakFreq_Power3_4.csv')
GammaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_all_PeakFreq_Power4_5.csv')
GammaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_all_PeakFreq_Power5_6.csv')

newSubsGammaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_newSubs_PeakFreq_Power1_2.csv')
newSubsGammaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_newSubs_PeakFreq_Power2_3.csv')
newSubsGammaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_newSubs_PeakFreq_Power3_4.csv')
newSubsGammaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_newSubs_PeakFreq_Power4_5.csv')
newSubsGammaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_newSubs_PeakFreq_Power5_6.csv')

GammaDelay1_2 <- rbind(GammaDelay1_2, newSubsGammaDelay1_2)
GammaDelay2_3 <- rbind(GammaDelay2_3, newSubsGammaDelay2_3)
GammaDelay3_4 <- rbind(GammaDelay3_4, newSubsGammaDelay3_4)
GammaDelay4_5 <- rbind(GammaDelay4_5, newSubsGammaDelay4_5)
GammaDelay5_6 <- rbind(GammaDelay5_6, newSubsGammaDelay5_6)


GammaDelay1_2$log_Peak_Power <- log1p(GammaDelay1_2$Gamma_Peak_Power)
GammaDelay2_3$log_Peak_Power <- log1p(GammaDelay2_3$Gamma_Peak_Power)
GammaDelay3_4$log_Peak_Power <- log1p(GammaDelay3_4$Gamma_Peak_Power)
GammaDelay4_5$log_Peak_Power <- log1p(GammaDelay4_5$Gamma_Peak_Power)
GammaDelay5_6$log_Peak_Power <- log1p(GammaDelay5_6$Gamma_Peak_Power)

write.csv(GammaDelay1_2, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_data_peaks_1_2_20200810.csv')
write.csv(GammaDelay2_3, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_all_data_peaks_2_3_20200810.csv')
write.csv(GammaDelay3_4, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_all_data_peaks_3_4_20200810.csv')
write.csv(GammaDelay4_5, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_all_data_peaks_4_5_20200810.csv')
write.csv(GammaDelay5_6, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_all_data_peaks_5_6_20200810.csv')




Avg_Gamma_Peak_Frequency1_2 <- aggregate(Gamma_Peak_Frequency ~ Subject, GammaDelay1_2, mean)
Avg_Gamma_Peak_Power1_2 <- aggregate(log_Peak_Power ~ Subject, GammaDelay1_2, mean)
Avg_GammaDelay1_2 <- merge (Avg_Gamma_Peak_Frequency1_2, Avg_Gamma_Peak_Power1_2, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Gamma_Peak_Frequency2_3 <- aggregate(Gamma_Peak_Frequency ~ Subject, GammaDelay2_3, mean)
Avg_Gamma_Peak_Power2_3 <- aggregate(log_Peak_Power ~ Subject, GammaDelay2_3, mean)
Avg_GammaDelay2_3 <- merge (Avg_Gamma_Peak_Frequency2_3, Avg_Gamma_Peak_Power2_3, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Gamma_Peak_Frequency3_4 <- aggregate(Gamma_Peak_Frequency ~ Subject, GammaDelay3_4, mean)
Avg_Gamma_Peak_Power3_4 <- aggregate(log_Peak_Power ~ Subject, GammaDelay3_4, mean)
Avg_GammaDelay3_4 <- merge (Avg_Gamma_Peak_Frequency3_4, Avg_Gamma_Peak_Power3_4, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Gamma_Peak_Frequency4_5 <- aggregate(Gamma_Peak_Frequency ~ Subject, GammaDelay4_5, mean)
Avg_Gamma_Peak_Power4_5 <- aggregate(log_Peak_Power ~ Subject, GammaDelay4_5, mean)
Avg_GammaDelay4_5 <- merge (Avg_Gamma_Peak_Frequency4_5, Avg_Gamma_Peak_Power4_5, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Gamma_Peak_Frequency5_6 <- aggregate(Gamma_Peak_Frequency ~ Subject, GammaDelay5_6, mean)
Avg_Gamma_Peak_Power5_6 <- aggregate(log_Peak_Power ~ Subject, GammaDelay5_6, mean)
Avg_GammaDelay5_6 <- merge (Avg_Gamma_Peak_Frequency5_6, Avg_Gamma_Peak_Power5_6, by = "Subject",  all.x = TRUE, all.y = TRUE)

allGammaPeakInfo <- rbind(Avg_GammaDelay1_2, Avg_GammaDelay2_3)
allGammaPeakInfo <- rbind(allGammaPeakInfo, Avg_GammaDelay3_4)
allGammaPeakInfo <- rbind(allGammaPeakInfo,Avg_GammaDelay4_5)
allGammaPeakInfo <- rbind(allGammaPeakInfo, Avg_GammaDelay5_6)
allGammaPeakInfo <- merge(allGammaPeakInfo, agefile, by= "Subject")

Gamma_Avg_log_Peak_Power <- aggregate(log_Peak_Power ~ Subject, allGammaPeakInfo, mean)
Gamma_Avg_Peak_Frequency <- aggregate(Gamma_Peak_Frequency ~ Subject, allGammaPeakInfo, mean)

gammaPeaks <-  merge(Gamma_Avg_log_Peak_Power, Gamma_Avg_Peak_Frequency, by = "Subject")
write.csv(gammaPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/delayGamma_Peaks_20200810.csv')



### Alpha Frequency Band 

AlphaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_all_PeakFreq_1_2.csv')
AlphaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_all_PeakFreq_2_3.csv')
AlphaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_all_PeakFreq_3_4.csv')
AlphaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_all_PeakFreq_4_5.csv')
AlphaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_all_PeakFreq_5_6.csv')

AlphaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_newSubs_PeakFreq_1_2.csv')
AlphaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_newSubs_PeakFreq_2_3.csv')
AlphaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_newSubs_PeakFreq_3_4.csv')
AlphaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_newSubs_PeakFreq_4_5.csv')
AlphaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_newSubs_PeakFreq_5_6.csv')

AlphaDelay1_2$log_Peak_Power <- log1p(AlphaDelay1_2$Alpha_Peak_Power)
AlphaDelay2_3$log_Peak_Power <- log1p(AlphaDelay2_3$Alpha_Peak_Power)
AlphaDelay3_4$log_Peak_Power <- log1p(AlphaDelay3_4$Alpha_Peak_Power)
AlphaDelay4_5$log_Peak_Power <- log1p(AlphaDelay4_5$Alpha_Peak_Power)
AlphaDelay5_6$log_Peak_Power <- log1p(AlphaDelay5_6$Alpha_Peak_Power)


Avg_Alpha_Peak_Frequency1_2 <- aggregate(Alpha_Peak_Frequency ~ Subject, AlphaDelay1_2, mean)
Avg_Alpha_Peak_Power1_2 <- aggregate(log_Peak_Power ~ Subject, AlphaDelay1_2, mean)
Avg_AlphaDelay1_2 <- merge (Avg_Alpha_Peak_Frequency1_2, Avg_Alpha_Peak_Power1_2, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Alpha_Peak_Frequency2_3 <- aggregate(Alpha_Peak_Frequency ~ Subject, AlphaDelay2_3, mean)
Avg_Alpha_Peak_Power2_3 <- aggregate(log_Peak_Power ~ Subject, AlphaDelay2_3, mean)
Avg_AlphaDelay2_3 <- merge (Avg_Alpha_Peak_Frequency2_3, Avg_Alpha_Peak_Power2_3, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Alpha_Peak_Frequency3_4 <- aggregate(Alpha_Peak_Frequency ~ Subject, AlphaDelay3_4, mean)
Avg_Alpha_Peak_Power3_4 <- aggregate(log_Peak_Power ~ Subject, AlphaDelay3_4, mean)
Avg_AlphaDelay3_4 <- merge (Avg_Alpha_Peak_Frequency3_4, Avg_Alpha_Peak_Power3_4, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Alpha_Peak_Frequency4_5 <- aggregate(Alpha_Peak_Frequency ~ Subject, AlphaDelay4_5, mean)
Avg_Alpha_Peak_Power4_5 <- aggregate(log_Peak_Power ~ Subject, AlphaDelay4_5, mean)
Avg_AlphaDelay4_5 <- merge (Avg_Alpha_Peak_Frequency4_5, Avg_Alpha_Peak_Power4_5, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Alpha_Peak_Frequency5_6 <- aggregate(Alpha_Peak_Frequency ~ Subject, AlphaDelay5_6, mean)
Avg_Alpha_Peak_Power5_6 <- aggregate(log_Peak_Power ~ Subject, AlphaDelay5_6, mean)
Avg_AlphaDelay5_6 <- merge (Avg_Alpha_Peak_Frequency5_6, Avg_Alpha_Peak_Power5_6, by = "Subject",  all.x = TRUE, all.y = TRUE)

allAlphaPeakInfo <- rbind(Avg_AlphaDelay1_2, Avg_AlphaDelay2_3)
allAlphaPeakInfo <- rbind(allAlphaPeakInfo, Avg_AlphaDelay3_4)
allAlphaPeakInfo <- rbind(allAlphaPeakInfo,Avg_AlphaDelay4_5)
allAlphaPeakInfo <- rbind(allAlphaPeakInfo, Avg_AlphaDelay5_6)
allAlphaPeakInfo <- merge(allAlphaPeakInfo, agefile, by= "Subject")

Alpha_Avg_log_Peak_Power <- aggregate(log_Peak_Power ~ Subject, allAlphaPeakInfo, mean)
Alpha_Avg_Peak_Frequency <- aggregate(Alpha_Peak_Frequency ~ Subject, allAlphaPeakInfo, mean)

AlphaPeaks <-  merge(Alpha_Avg_log_Peak_Power, Alpha_Avg_Peak_Frequency, by = "Subject")
write.csv(AlphaPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/delayAlpha_Peaks_newSubs.csv')


### Theta Frequency Band 

ThetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_all_PeakFreq_Power1_2.csv')
ThetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_all_PeakFreq_Power2_3.csv')
ThetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_all_PeakFreq_Power3_4.csv')
ThetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_all_PeakFreq_Power4_5.csv')
ThetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_all_PeakFreq_Power5_6.csv')

ThetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_all_PeakFreq_Power1_2.csv')
ThetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_all_PeakFreq_Power2_3.csv')
ThetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_all_PeakFreq_Power3_4.csv')
ThetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_all_PeakFreq_Power4_5.csv')
ThetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_all_PeakFreq_Power5_6.csv')


ThetaDelay1_2$log_Peak_Power <- log1p(ThetaDelay1_2$Theta_Peak_Power)
ThetaDelay2_3$log_Peak_Power <- log1p(ThetaDelay2_3$Theta_Peak_Power)
ThetaDelay3_4$log_Peak_Power <- log1p(ThetaDelay3_4$Theta_Peak_Power)
ThetaDelay4_5$log_Peak_Power <- log1p(ThetaDelay4_5$Theta_Peak_Power)
ThetaDelay5_6$log_Peak_Power <- log1p(ThetaDelay5_6$Theta_Peak_Power)


Avg_Theta_Peak_Frequency1_2 <- aggregate(Theta_Peak_Frequency ~ Subject, ThetaDelay1_2, mean)
Avg_Theta_Peak_Power1_2 <- aggregate(log_Peak_Power ~ Subject, ThetaDelay1_2, mean)
Avg_ThetaDelay1_2 <- merge (Avg_Theta_Peak_Frequency1_2, Avg_Theta_Peak_Power1_2, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Theta_Peak_Frequency2_3 <- aggregate(Theta_Peak_Frequency ~ Subject, ThetaDelay2_3, mean)
Avg_Theta_Peak_Power2_3 <- aggregate(log_Peak_Power ~ Subject, ThetaDelay2_3, mean)
Avg_ThetaDelay2_3 <- merge (Avg_Theta_Peak_Frequency2_3, Avg_Theta_Peak_Power2_3, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Theta_Peak_Frequency3_4 <- aggregate(Theta_Peak_Frequency ~ Subject, ThetaDelay3_4, mean)
Avg_Theta_Peak_Power3_4 <- aggregate(log_Peak_Power ~ Subject, ThetaDelay3_4, mean)
Avg_ThetaDelay3_4 <- merge (Avg_Theta_Peak_Frequency3_4, Avg_Theta_Peak_Power3_4, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Theta_Peak_Frequency4_5 <- aggregate(Theta_Peak_Frequency ~ Subject, ThetaDelay4_5, mean)
Avg_Theta_Peak_Power4_5 <- aggregate(log_Peak_Power ~ Subject, ThetaDelay4_5, mean)
Avg_ThetaDelay4_5 <- merge (Avg_Theta_Peak_Frequency4_5, Avg_Theta_Peak_Power4_5, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Theta_Peak_Frequency5_6 <- aggregate(Theta_Peak_Frequency ~ Subject, ThetaDelay5_6, mean)
Avg_Theta_Peak_Power5_6 <- aggregate(log_Peak_Power ~ Subject, ThetaDelay5_6, mean)
Avg_ThetaDelay5_6 <- merge (Avg_Theta_Peak_Frequency5_6, Avg_Theta_Peak_Power5_6, by = "Subject",  all.x = TRUE, all.y = TRUE)

allThetaPeakInfo <- rbind(Avg_ThetaDelay1_2, Avg_ThetaDelay2_3)
allThetaPeakInfo <- rbind(allThetaPeakInfo, Avg_ThetaDelay3_4)
allThetaPeakInfo <- rbind(allThetaPeakInfo,Avg_ThetaDelay4_5)
allThetaPeakInfo <- rbind(allThetaPeakInfo, Avg_ThetaDelay5_6)
allThetaPeakInfo <- merge(allThetaPeakInfo, agefile, by= "Subject")

Theta_Avg_log_Peak_Power <- aggregate(log_Peak_Power ~ Subject, allThetaPeakInfo, mean)
Theta_Avg_Peak_Frequency <- aggregate(Theta_Peak_Frequency ~ Subject, allThetaPeakInfo, mean)

ThetaPeaks <-  merge(Theta_Avg_log_Peak_Power, Theta_Avg_Peak_Frequency, by = "Subject")
write.csv(ThetaPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/delayTheta_Peaks_newSUbs.csv')


### Beta Frequency Band 

BetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_PeakFreq_Power1_2.csv')
BetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_all_PeakFreq_Power2_3.csv')
BetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_all_PeakFreq_Power3_4.csv')
BetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_all_PeakFreq_Power4_5.csv')
BetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_all_PeakFreq_Power5_6.csv')

newSubsBetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_newSubs_PeakFreq_Power1_2.csv')
newSubsBetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_newSubs_PeakFreq_Power2_3.csv')
newSubsBetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_newSubs_PeakFreq_Power3_4.csv')
newSubsBetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_newSubs_PeakFreq_Power4_5.csv')
newSubsBetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_newSubs_PeakFreq_Power5_6.csv')

BetaDelay1_2 <- rbind(BetaDelay1_2, newSubsBetaDelay1_2)
BetaDelay2_3 <- rbind(BetaDelay2_3, newSubsBetaDelay2_3)
BetaDelay3_4 <- rbind(BetaDelay3_4, newSubsBetaDelay3_4)
BetaDelay4_5 <- rbind(BetaDelay4_5, newSubsBetaDelay4_5)
BetaDelay5_6 <- rbind(BetaDelay5_6, newSubsBetaDelay5_6)

BetaDelay1_2$log_Peak_Power <- log1p(BetaDelay1_2$Beta_Peak_Power)
BetaDelay2_3$log_Peak_Power <- log1p(BetaDelay2_3$Beta_Peak_Power)
BetaDelay3_4$log_Peak_Power <- log1p(BetaDelay3_4$Beta_Peak_Power)
BetaDelay4_5$log_Peak_Power <- log1p(BetaDelay4_5$Beta_Peak_Power)
BetaDelay5_6$log_Peak_Power <- log1p(BetaDelay5_6$Beta_Peak_Power)


write.csv(BetaDelay1_2, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_data_peaks_1_2_20200810.csv')
write.csv(BetaDelay2_3, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_all_data_peaks_2_3_20200810.csv')
write.csv(BetaDelay3_4, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_all_data_peaks_3_4_20200810.csv')
write.csv(BetaDelay4_5, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_all_data_peaks_4_5_20200810.csv')
write.csv(BetaDelay5_6, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_all_data_peaks_5_6_20200810.csv')




Avg_Beta_Peak_Frequency1_2 <- aggregate(Beta_Peak_Frequency ~ Subject, BetaDelay1_2, mean)
Avg_Beta_Peak_Power1_2 <- aggregate(log_Peak_Power ~ Subject, BetaDelay1_2, mean)
Avg_BetaDelay1_2 <- merge (Avg_Beta_Peak_Frequency1_2, Avg_Beta_Peak_Power1_2, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Beta_Peak_Frequency2_3 <- aggregate(Beta_Peak_Frequency ~ Subject, BetaDelay2_3, mean)
Avg_Beta_Peak_Power2_3 <- aggregate(log_Peak_Power ~ Subject, BetaDelay2_3, mean)
Avg_BetaDelay2_3 <- merge (Avg_Beta_Peak_Frequency2_3, Avg_Beta_Peak_Power2_3, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Beta_Peak_Frequency3_4 <- aggregate(Beta_Peak_Frequency ~ Subject, BetaDelay3_4, mean)
Avg_Beta_Peak_Power3_4 <- aggregate(log_Peak_Power ~ Subject, BetaDelay3_4, mean)
Avg_BetaDelay3_4 <- merge (Avg_Beta_Peak_Frequency3_4, Avg_Beta_Peak_Power3_4, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Beta_Peak_Frequency4_5 <- aggregate(Beta_Peak_Frequency ~ Subject, BetaDelay4_5, mean)
Avg_Beta_Peak_Power4_5 <- aggregate(log_Peak_Power ~ Subject, BetaDelay4_5, mean)
Avg_BetaDelay4_5 <- merge (Avg_Beta_Peak_Frequency4_5, Avg_Beta_Peak_Power4_5, by = "Subject",  all.x = TRUE, all.y = TRUE)


Avg_Beta_Peak_Frequency5_6 <- aggregate(Beta_Peak_Frequency ~ Subject, BetaDelay5_6, mean)
Avg_Beta_Peak_Power5_6 <- aggregate(log_Peak_Power ~ Subject, BetaDelay5_6, mean)
Avg_BetaDelay5_6 <- merge (Avg_Beta_Peak_Frequency5_6, Avg_Beta_Peak_Power5_6, by = "Subject",  all.x = TRUE, all.y = TRUE)

allBetaPeakInfo <- rbind(Avg_BetaDelay1_2, Avg_BetaDelay2_3)
allBetaPeakInfo <- rbind(allBetaPeakInfo, Avg_BetaDelay3_4)
allBetaPeakInfo <- rbind(allBetaPeakInfo,Avg_BetaDelay4_5)
allBetaPeakInfo <- rbind(allBetaPeakInfo, Avg_BetaDelay5_6)
allBetaPeakInfo <- merge(allBetaPeakInfo, agefile, by= "Subject")

Beta_Avg_log_Peak_Power <- aggregate(log_Peak_Power ~ Subject, allBetaPeakInfo, mean)
Beta_Avg_Peak_Frequency <- aggregate(Beta_Peak_Frequency ~ Subject, allBetaPeakInfo, mean)

BetaPeaks <-  merge(Beta_Avg_log_Peak_Power, Beta_Avg_Peak_Frequency, by = "Subject")
write.csv(BetaPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/delayBeta_Peaks_20200810.csv')

}

Create_CSV_files_delay_and_fix <- function () {
  #gamma 
  GammaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_3_4.csv')
  GammaDelay3_4$log1p_Gamma_Trial_Power <- log1p(GammaDelay3_4$Gamma_Trial_Power)
  
  GammaDelay3_4_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_PeakFreq_Power3_4.csv')
  GammaDelay3_4_peaks$log1p_Gamma_Peak_Power <- log1p(GammaDelay3_4_peaks$Gamma_Peak_Power)
  
  GammaDelay3_4_mean <- setnames(aggregate(. ~ Subject, data = GammaDelay3_4, mean), c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  GammaDelay3_4_peaks_mean <- setnames(aggregate(. ~ Subject, data = GammaDelay3_4_peaks, mean), c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power"))
  
  allGammaDelay_3_4 <- merge(GammaDelay3_4_mean, GammaDelay3_4_peaks_mean, by = "Subject")
  
  GammaDelay3_4_var <- setnames(aggregate(. ~ Subject, data = GammaDelay3_4, sd), c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  GammaDelay3_4_peaks_var <- setnames(aggregate(. ~ Subject, data = GammaDelay3_4_peaks, sd), c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))
  
  allGammaDelay_3_4_var <- merge(GammaDelay3_4_var, GammaDelay3_4_peaks_var, by = "Subject")
  
  allGammaDelay_3_4 <- merge(allGammaDelay_3_4, allGammaDelay_3_4_var, by = "Subject")
  
  GammaFix <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Subs_Fix.csv')
  GammaFix$log1p_Gamma_Trial_Power <- log1p(GammaFix$Gamma_Trial_Power)
  GammaFix_mean <- setnames(aggregate(. ~ Subject, data = GammaFix, mean), c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  GammaFix_var <- setnames(aggregate(. ~ Subject, data = GammaFix, sd), c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  
  
  GammaFix_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Subs_PeakFreq_PowerFix.csv')
  GammaFix_peaks$log1p_Gamma_Peak_Power <- log1p(GammaFix_peaks$Gamma_Peak_Power)
  GammaFix_peaks_mean <- setnames(aggregate(. ~ Subject, data = GammaFix_peaks, mean), c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power"))
  GammaFix_peaks_var <- setnames(aggregate(. ~ Subject, data = GammaFix_peaks, sd), c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))
  
  
  allGammaFix_mean <- merge(GammaFix_mean, GammaFix_peaks_mean, by = "Subject")
  allGammaFix_var <- merge(GammaFix_var, GammaFix_peaks_var, by = "Subject")
  
  allGammaFix <- merge(allGammaFix_mean, allGammaFix_var, by = "Subject")
  
  allGamma_3_4 <- merge(allGammaDelay_3_4, allGammaFix, by = "Subject", suffix = c("_Delay", "_Fix"))

  write.csv(allGamma_3_4,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/allGamma_3_4_20200821.csv')
  
  
  #alpha 
  AlphaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_data_3_4.csv')
  AlphaDelay3_4$log1p_Alpha_Trial_Power <- log1p(AlphaDelay3_4$Alpha_Trial_Power)
  
  AlphaDelay3_4_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_PeakFreq_3_4.csv')
  AlphaDelay3_4_peaks$log1p_Alpha_Peak_Power <- log1p(AlphaDelay3_4_peaks$Alpha_Peak_Power)
  
  AlphaDelay3_4_mean <- setnames(aggregate(. ~ Subject, data = AlphaDelay3_4, mean), c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  AlphaDelay3_4_var <- setnames(aggregate(. ~ Subject, data = AlphaDelay3_4, sd), c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  
  AlphaDelay3_4_peaks_mean <- setnames(aggregate(. ~ Subject, data = AlphaDelay3_4_peaks, mean), c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power")) 
  AlphaDelay3_4_peaks_var <- setnames(aggregate(. ~ Subject, data = AlphaDelay3_4_peaks, sd), c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))

  
  allAlphaDelay_3_4 <- merge(AlphaDelay3_4_mean, AlphaDelay3_4_peaks_mean, by = "Subject")
  allAlphaDelay_3_4_var <- merge(AlphaDelay3_4_var, AlphaDelay3_4_peaks_var, by = "Subject")
  
  allAlphaDelay_3_4 <- merge(allAlphaDelay_3_4, allAlphaDelay_3_4_var, by = "Subject")
  
  
  
  AlphaFix <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Subs_datafix.csv')
  AlphaFix$log1p_Alpha_Trial_Power <- log1p(AlphaFix$Alpha_Trial_Power)
  AlphaFix_mean <- setnames(aggregate(. ~ Subject, data = AlphaFix, mean), c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  AlphaFix_var <- setnames(aggregate(. ~ Subject, data = AlphaFix, sd), c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  
  
  AlphaFix_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Subs_PeakFreq_Fix.csv')
  AlphaFix_peaks$log1p_Alpha_Peak_Power <- log1p(AlphaFix_peaks$Alpha_Peak_Power)
  AlphaFix_peaks_mean <- setnames(aggregate(. ~ Subject, data = AlphaFix_peaks, mean),  c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power"))
  AlphaFix_peaks_var <- setnames(aggregate(. ~ Subject, data = AlphaFix_peaks, sd),  c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))
  
  
  allAlphaFix <- merge(AlphaFix_mean, AlphaFix_peaks_mean, by = "Subject")
  allAlphaFix_Var <- merge(AlphaFix_var, AlphaFix_peaks_var, by = "Subject")
  
  allAlphaFix <- merge(allAlphaFix, allAlphaFix_Var, by = "Subject")
  
  allAlpha_3_4 <- merge(allAlphaDelay_3_4, allAlphaFix, by = "Subject", suffix = c("_Delay", "_Fix"))
  write.csv(allAlpha_3_4,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/allAlpha_3_4_20200821.csv')
  
  
  #Theta 
  ThetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_data3_4.csv')
  ThetaDelay3_4$log1p_Theta_Trial_Power <- log1p(ThetaDelay3_4$Theta_Trial_Power)
  
  ThetaDelay3_4_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_PeakFreq_Power3_4.csv')
  ThetaDelay3_4_peaks$log1p_Theta_Peak_Power <- log1p(ThetaDelay3_4_peaks$Theta_Peak_Power)
  
  ThetaDelay3_4_mean <- setnames(aggregate(. ~ Subject, data = ThetaDelay3_4, mean),  c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  ThetaDelay3_4_var <- setnames(aggregate(. ~ Subject, data = ThetaDelay3_4, sd),  c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  
  ThetaDelay3_4_peaks_mean <- setnames(aggregate(. ~ Subject, data = ThetaDelay3_4_peaks, mean),  c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power"))
  ThetaDelay3_4_peaks_var <- setnames(aggregate(. ~ Subject, data = ThetaDelay3_4_peaks, sd),  c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))
  
  
  allThetaDelay_3_4 <- merge(ThetaDelay3_4_mean, ThetaDelay3_4_peaks_mean, by = "Subject")
  allThetaDelay_3_4_var <- merge(ThetaDelay3_4_var, ThetaDelay3_4_peaks_var, by = "Subject")
  allThetaDelay_3_4 <- merge(allThetaDelay_3_4, allThetaDelay_3_4_var, by = "Subject")
  
  
  
  ThetaFix <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Subs_dataFix.csv')
  ThetaFix$log1p_Theta_Trial_Power <- log1p(ThetaFix$Theta_Trial_Power)
  ThetaFix_mean <- setnames(aggregate(. ~ Subject, data = ThetaFix, mean),  c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  ThetaFix_var <- setnames(aggregate(. ~ Subject, data = ThetaFix, sd),  c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  
  
  ThetaFix_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Subs_PeakFreq_PowerFix.csv')
  ThetaFix_peaks$log1p_Theta_Peak_Power <- log1p(ThetaFix_peaks$Theta_Peak_Power)
  ThetaFix_peaks_mean <- setnames(aggregate(. ~ Subject, data = ThetaFix_peaks, mean),  c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power"))
  ThetaFix_peaks_var <- setnames(aggregate(. ~ Subject, data = ThetaFix_peaks, sd),  c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))
  
  allThetaFix <- merge(ThetaFix_mean, ThetaFix_peaks_mean, by = "Subject")
  allThetaFix_var <- merge(ThetaFix_var, ThetaFix_peaks_var, by = "Subject")
  allThetaFix <- merge(allThetaFix, allThetaFix_var, by = "Subject")
  
  
  allTheta_3_4 <- merge(allThetaDelay_3_4, allThetaFix, by = "Subject", suffix = c("_Delay", "_Fix"))
  write.csv(allTheta_3_4,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/allTheta_3_4_20200821.csv')
  
  #Beta 
  BetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_data3_4.csv')
  BetaDelay3_4$log1p_Beta_Trial_Power <- log1p(BetaDelay3_4$Beta_Trial_Power)
  
  BetaDelay3_4_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_PeakFreq_Power_3_4.csv')
  BetaDelay3_4_peaks$log1p_Beta_Peak_Power <- log1p(BetaDelay3_4_peaks$Beta_Peak_Power)
  
  BetaDelay3_4_mean <- setnames(aggregate(. ~ Subject, data = BetaDelay3_4, mean), c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  BetaDelay3_4_var <- setnames(aggregate(. ~ Subject, data = BetaDelay3_4, sd), c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  
  BetaDelay3_4_peaks_mean <- setnames(aggregate(. ~ Subject, data = BetaDelay3_4_peaks, mean), c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power"))
  BetaDelay3_4_peaks_var <- setnames(aggregate(. ~ Subject, data = BetaDelay3_4_peaks, sd), c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))
  
  
  allBetaDelay_3_4 <- merge(BetaDelay3_4_mean, BetaDelay3_4_peaks_mean, by = "Subject")
  allBetaDelay_3_4_var <- merge(BetaDelay3_4_var, BetaDelay3_4_peaks_var, by = "Subject")
  allBetaDelay_3_4 <- merge(allBetaDelay_3_4, allBetaDelay_3_4_var, by = "Subject")
  
  
  BetaFix <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Subs_dataFix.csv')
  BetaFix$log1p_Beta_Trial_Power <- log1p(BetaFix$Beta_Trial_Power)
  BetaFix_mean <- setnames(aggregate(. ~ Subject, data = BetaFix, mean), c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power"))
  BetaFix_var <- setnames(aggregate(. ~ Subject, data = BetaFix, sd), c("Subject", "Trial", "Trial_Power_Variability", "Event_Number_Variability", "Event_Duration_Variability", "log1p_Trial_Power_Variability"))
  
  BetaFix_peaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Subs_PeakFreq_PowerFix.csv')
  BetaFix_peaks$log1p_Beta_Peak_Power <- log1p(BetaFix_peaks$Beta_Peak_Power)
  BetaFix_peaks_mean <- setnames(aggregate(. ~ Subject, data = BetaFix_peaks, mean), c("Subject", "Trial", "Peak_Frequency", "Peak_Power", "log1p_Peak_Power"))
  BetaFix_peaks_var <- setnames(aggregate(. ~ Subject, data = BetaFix_peaks, sd), c("Subject", "Trial", "Peak_Frequency_Variability", "Peak_Power_Variability", "log1p_Peak_Power_Variability"))
  
  allBetaFix <- merge(BetaFix_mean, BetaFix_peaks_mean, by = "Subject")
  allBetaFix_var <- merge(BetaFix_var, BetaFix_peaks_var, by = "Subject")
  allBetaFix <- merge(allBetaFix, allBetaFix_var, by = "Subject")
  
  
  allBeta_3_4 <- merge(allBetaDelay_3_4, allBetaFix, by = "Subject", suffix = c("_Delay", "_Fix"))
  write.csv(allBeta_3_4,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/allBeta_3_4_20200821.csv')
  
}

Only_Take_One_Delay_Bin_TrialLevel <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_data3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_data_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_data3_4.csv' )
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
  ## ---- alldata ----
  alldata_TrialLevel <- merge(alpha[,c("Subject","Trial",alphavars)],beta[,c("Subject","Trial",betavars)],by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,theta[,c("Subject","Trial",thetavars)],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,gamma[,c("Subject","Trial",gammavars)],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  alldata_TrialLevel$Gamma.log1p_Trial_Power <- log1p(alldata_TrialLevel$Gamma.Trial_Power)
  alldata_TrialLevel$Beta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Beta.Trial_Power)
  alldata_TrialLevel$Alpha.log1p_Trial_Power <- log1p(alldata_TrialLevel$Alpha.Trial_Power)
  alldata_TrialLevel$Theta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Theta.Trial_Power)
  
  
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv')
  agefile$inverseAge <- 1/agefile$age
  
  alldata_TrialLevel <- merge(alldata_TrialLevel, agefile, by = c("Subject"))
  alldata_TrialLevel <- alldata_TrialLevel[alldata_TrialLevel$visitno < 2,]
  alldata_TrialLevel <- na.omit(alldata_TrialLevel)
  
  
  
  write.csv(alldata_TrialLevel,"H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/alldata_TrialLevel_20210528.csv")
  
  
  
  # Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_TrialLevel_new <- alldata_TrialLevel
  
  cols = names(alldata_TrialLevel[3:18])
  for (col in cols) {
    
    alldata_TrialLevel_grouped <- alldata_TrialLevel %>% group_by(Subject)
    
    indx <- outliers(alldata_TrialLevel_grouped[[col]])
    
    alldata_TrialLevel_new[[col]] <- Map(replace, alldata_TrialLevel_new[[col]], indx, NA)
    alldata_TrialLevel_new[[col]] <- as.numeric(alldata_TrialLevel_new[[col]])
    
  }  
  
  alldf$alldata_TrialLevel_new <- alldata_TrialLevel_new
  
  # Find the percent of NAs per subject per variable 

  percentNAs <- alldata_TrialLevel_new %>% group_by(Subject) %>% dplyr::select(Subject, Alpha.Trial_Power:Theta.log1p_Trial_Power) %>% summarise_all(funs(pctNAs = sum(is.na(.))/length(.)* 100) )
  
  return(alldf)
  
}


DelayOnly_Sublevel <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_data3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_data_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_data3_4.csv' )
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
  ## ---- alldata ----
  alldata_TrialLevel <- merge(alpha[,c("Subject","Trial",alphavars)],beta[,c("Subject","Trial",betavars)],by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,theta[,c("Subject","Trial",thetavars)],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,gamma[,c("Subject","Trial",gammavars)],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  alldata_TrialLevel$Gamma.log1p_Trial_Power <- log1p(alldata_TrialLevel$Gamma.Trial_Power)
  alldata_TrialLevel$Beta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Beta.Trial_Power)
  alldata_TrialLevel$Alpha.log1p_Trial_Power <- log1p(alldata_TrialLevel$Alpha.Trial_Power)
  alldata_TrialLevel$Theta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Theta.Trial_Power)
  
  # Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_TrialLevel_new <- alldata_TrialLevel
  
  cols = names(alldata_TrialLevel[3:18])
  for (col in cols) {
    
    alldata_TrialLevel_grouped <- alldata_TrialLevel %>% group_by(Subject)
    
    indx <- outliers(alldata_TrialLevel_grouped[[col]])
    
    alldata_TrialLevel_new[[col]] <- Map(replace, alldata_TrialLevel_new[[col]], indx, NA)
    alldata_TrialLevel_new[[col]] <- as.numeric(alldata_TrialLevel_new[[col]])
    
  }  
  
  alldata_subLevel_avg <- aggregate(.~ Subject, alldata_TrialLevel_new , mean)
  alldata_subLevel_sd <- aggregate(.~ Subject, alldata_TrialLevel_new , sd)
  alldata_subLevel <- merge(alldata_subLevel_avg, alldata_subLevel_sd, by = "Subject", suffixes = c("", "_Variability"))
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  alldata_delayOnly_age <- merge(alldata_subLevel, agefile, by = "Subject")
  
  alldata_delayOnly_age <- subset(alldata_delayOnly_age, visitno == 1)
  
  alldata_delayOnly_age$inverseAge <- 1/alldata_delayOnly_age$age
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_delayOnly_age_new <- alldata_delayOnly_age
  
  cols = names(alldata_delayOnly_age[3:35])
  for (col in cols) {
    
    indx <- outliers(alldata_delayOnly_age[[col]])
    
    alldata_delayOnly_age_new[[col]] <- Map(replace, alldata_delayOnly_age_new[[col]], indx, NA)
    alldata_delayOnly_age_new[[col]] <- as.numeric(alldata_delayOnly_age_new[[col]])
  }
  
  
  
  alldf$alldata_delayOnly <- alldata_delayOnly_age_new
  return(alldf)
  
}


Only_Take_One_Delay_Bin_and_Fix <- function () {
  
  library(knitr)
  
  alldf <- list()
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200821.csv')
  
  Behavior$log_absPositionError <- log1p(abs(Behavior$PositionError))
  Behavior$log_mgsLatency <- log1p(Behavior$mgsLatency)
  behavior <- setnames(aggregate(list(Behavior$log_mgsLatency, Behavior$log_absPositionError), list(Subject =  paste0(Behavior$LunaID, '_', Behavior$ScanDate)), na.rm = TRUE, mean), c("Subject", "log_mgsLatency", "log_absPositionError"))
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/allTheta_3_4_20200821.csv')
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta)[names(theta) %in% analysisvars] <- paste("Theta",names(theta[,analysisvars]),sep=".")
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/allAlpha_3_4_20200821.csv')
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha)[names(alpha) %in% analysisvars] <- paste("Alpha",names(alpha[,analysisvars]),sep=".")
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/allBeta_3_4_20200821.csv' )
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta)[names(beta) %in% analysisvars] <- paste("Beta",names(beta[,analysisvars]),sep=".")
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/allGamma_3_4_20200821.csv')
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma)[names(gamma) %in% analysisvars] <- paste("Gamma",names(gamma[,analysisvars]),sep=".") ##add Gamma to analysis var names
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  alldf$allGamma <- allGamma
  
  ## ---- alldata ----
  alldata <- merge(alpha[,c("Subject",alphavars)],beta[,c("Subject",betavars)],by="Subject", all.x = TRUE, all.y = TRUE) %>% 
    merge(.,theta[,c("Subject",thetavars)],by="Subject", all.x = TRUE, all.y = TRUE) %>% merge(.,gamma[,c("Subject",gammavars)],by="Subject", all.x = TRUE, all.y = TRUE)
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  agefile$inverseAge <- 1/agefile$age
  
  alldata <- merge(alldata, agefile, by = "Subject",  all.x = TRUE, all.y = TRUE)
  alldata <- alldata[alldata$visitno <2,]
  
  alldata_complete <- alldata[complete.cases(alldata),] #this will delete incomplete cases 

  alldf$alldata <- alldata_complete
  
  return(alldf)
  
  
}

IndividualTrialData <- function() {
  
  #gamma
  GammaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_data_1_2.csv')
  GammaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_all_data_2_3.csv')
  GammaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_all_data_3_4.csv')
  GammaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_all_data_4_5.csv')
  GammaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_all_data_5_6.csv')
  
  GammaDelay1_2$log1p_Gamma_Trial_Power <- log1p(GammaDelay1_2$Gamma_Trial_Power)
  GammaDelay2_3$log1p_Gamma_Trial_Power <- log1p(GammaDelay2_3$Gamma_Trial_Power)
  GammaDelay3_4$log1p_Gamma_Trial_Power <- log1p(GammaDelay3_4$Gamma_Trial_Power)
  GammaDelay4_5$log1p_Gamma_Trial_Power <- log1p(GammaDelay4_5$Gamma_Trial_Power)
  GammaDelay5_6$log1p_Gamma_Trial_Power <- log1p(GammaDelay5_6$Gamma_Trial_Power)
  
  GammaDelay1_2$epoch <- "1-2"
  GammaDelay2_3$epoch <- "2-3"
  GammaDelay3_4$epoch <- "3-4"
  GammaDelay4_5$epoch <- "4-5"
  GammaDelay5_6$epoch <- "5-6"
  
  allGammaTrialData <- rbind(GammaDelay1_2, GammaDelay2_3)
  allGammaTrialData <- rbind(allGammaTrialData, GammaDelay3_4)
  allGammaTrialData <- rbind(allGammaTrialData, GammaDelay4_5)
  allGammaTrialData <- rbind(allGammaTrialData, GammaDelay5_6)
  
  
  Avg_allGammaTrialData <- allGammaTrialData[] %>% dplyr::group_by(Subject,Trial) %>% summarize_all(mean)
  
  setnames(Avg_allGammaTrialData, c("Subject","Trial", "Trial_Power", "Event_Number","Event_Duration", "log1p_Trial_Power", "epoch"))
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(Avg_allGammaTrialData),value=TRUE) ###matches names

  names(Avg_allGammaTrialData)[names(Avg_allGammaTrialData) %in% analysisvars] <- paste("Gamma",names(Avg_allGammaTrialData[,analysisvars]),sep=".") ##add Gamma to analysis var names
  gammavars <- grep("Gamma",names(Avg_allGammaTrialData),value=TRUE)
  
  allGammaTrials <- Avg_allGammaTrialData[,c("Subject","Trial",gammavars)]
  allGammaTrials <- allGammaTrials[complete.cases(allGammaTrials),]
  allGammaTrials <- allGammaTrials[!duplicated(allGammaTrials),]
  
  write.csv(allGammaTrials, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/delayGamma_Trial.csv')
  
  #Beta
  BetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_data1_2.csv')
  BetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_all_data2_3.csv')
  BetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_all_data3_4.csv')
  BetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_all_data4_5.csv')
  BetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_all_data5_6.csv')
  
  BetaDelay1_2$log1p_Beta_Trial_Power <- log1p(BetaDelay1_2$Beta_Trial_Power)
  BetaDelay2_3$log1p_Beta_Trial_Power <- log1p(BetaDelay2_3$Beta_Trial_Power)
  BetaDelay3_4$log1p_Beta_Trial_Power <- log1p(BetaDelay3_4$Beta_Trial_Power)
  BetaDelay4_5$log1p_Beta_Trial_Power <- log1p(BetaDelay4_5$Beta_Trial_Power)
  BetaDelay5_6$log1p_Beta_Trial_Power <- log1p(BetaDelay5_6$Beta_Trial_Power)
  
  BetaDelay1_2$epoch <- "1-2"
  BetaDelay2_3$epoch <- "2-3"
  BetaDelay3_4$epoch <- "3-4"
  BetaDelay4_5$epoch <- "4-5"
  BetaDelay5_6$epoch <- "5-6"
  
  allBetaTrialData <- rbind(BetaDelay1_2, BetaDelay2_3)
  allBetaTrialData <- rbind(allBetaTrialData, BetaDelay3_4)
  allBetaTrialData <- rbind(allBetaTrialData, BetaDelay4_5)
  allBetaTrialData <- rbind(allBetaTrialData, BetaDelay5_6)
  
  
  Avg_allBetaTrialData <- allBetaTrialData[] %>% dplyr::group_by(Subject,Trial) %>% summarize_all(mean)
  
  setnames(Avg_allBetaTrialData, c("Subject","Trial", "Trial_Power", "Event_Number","Event_Duration", "log1p_Trial_Power", "epoch"))
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(Avg_allBetaTrialData),value=TRUE) ###matches names
  
  names(Avg_allBetaTrialData)[names(Avg_allBetaTrialData) %in% analysisvars] <- paste("Beta",names(Avg_allBetaTrialData[,analysisvars]),sep=".") ##add Beta to analysis var names
  Betavars <- grep("Beta",names(Avg_allBetaTrialData),value=TRUE)
  
  allBetaTrials <- Avg_allBetaTrialData[,c("Subject","Trial",Betavars)]
  allBetaTrials <- allBetaTrials[complete.cases(allBetaTrials),]
  allBetaTrials <- allBetaTrials[!duplicated(allBetaTrials),]
  
  write.csv(allBetaTrials, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/delayBeta_Trial.csv')  
  
  #Alpha
  AlphaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_all_data1_2.csv')
  AlphaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_all_data2_3.csv')
  AlphaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_all_data3_4.csv')
  AlphaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_all_data4_5.csv')
  AlphaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_all_data5_6.csv')
  
  AlphaDelay1_2$log1p_Alpha_Trial_Power <- log1p(AlphaDelay1_2$Alpha_Trial_Power)
  AlphaDelay2_3$log1p_Alpha_Trial_Power <- log1p(AlphaDelay2_3$Alpha_Trial_Power)
  AlphaDelay3_4$log1p_Alpha_Trial_Power <- log1p(AlphaDelay3_4$Alpha_Trial_Power)
  AlphaDelay4_5$log1p_Alpha_Trial_Power <- log1p(AlphaDelay4_5$Alpha_Trial_Power)
  AlphaDelay5_6$log1p_Alpha_Trial_Power <- log1p(AlphaDelay5_6$Alpha_Trial_Power)
  
  AlphaDelay1_2$epoch <- "1-2"
  AlphaDelay2_3$epoch <- "2-3"
  AlphaDelay3_4$epoch <- "3-4"
  AlphaDelay4_5$epoch <- "4-5"
  AlphaDelay5_6$epoch <- "5-6"
  
  allAlphaTrialData <- rbind(AlphaDelay1_2, AlphaDelay2_3)
  allAlphaTrialData <- rbind(allAlphaTrialData, AlphaDelay3_4)
  allAlphaTrialData <- rbind(allAlphaTrialData, AlphaDelay4_5)
  allAlphaTrialData <- rbind(allAlphaTrialData, AlphaDelay5_6)
  
  
  Avg_allAlphaTrialData <- allAlphaTrialData[] %>% dplyr::group_by(Subject,Trial) %>% summarize_all(mean)

  setnames(Avg_allAlphaTrialData, c("Subject","Trial", "Trial_Power", "Event_Number","Event_Duration", "log1p_Trial_Power", "epoch"))
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(Avg_allAlphaTrialData),value=TRUE) ###matches names
  
  names(Avg_allAlphaTrialData)[names(Avg_allAlphaTrialData) %in% analysisvars] <- paste("Alpha",names(Avg_allAlphaTrialData[,analysisvars]),sep=".") ##add Alpha to analysis var names
  Alphavars <- grep("Alpha",names(Avg_allAlphaTrialData),value=TRUE)
  
  allAlphaTrials <- Avg_allAlphaTrialData[,c("Subject","Trial",Alphavars)]
  allAlphaTrials <- allAlphaTrials[complete.cases(allAlphaTrials),]
  allAlphaTrials <- allAlphaTrials[!duplicated(allAlphaTrials),]
  
  write.csv(allAlphaTrials, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/delayAlpha_Trial.csv')
 
   #Theta
  ThetaDelay1_2 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_all_data1_2.csv')
  ThetaDelay2_3 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_all_data2_3.csv')
  ThetaDelay3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_all_data3_4.csv')
  ThetaDelay4_5 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_all_data4_5.csv')
  ThetaDelay5_6 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_all_data5_6.csv')
  
  ThetaDelay1_2$log1p_Theta_Trial_Power <- log1p(ThetaDelay1_2$Theta_Trial_Power)
  ThetaDelay2_3$log1p_Theta_Trial_Power <- log1p(ThetaDelay2_3$Theta_Trial_Power)
  ThetaDelay3_4$log1p_Theta_Trial_Power <- log1p(ThetaDelay3_4$Theta_Trial_Power)
  ThetaDelay4_5$log1p_Theta_Trial_Power <- log1p(ThetaDelay4_5$Theta_Trial_Power)
  ThetaDelay5_6$log1p_Theta_Trial_Power <- log1p(ThetaDelay5_6$Theta_Trial_Power)
  
  ThetaDelay1_2$epoch <- "1-2"
  ThetaDelay2_3$epoch <- "2-3"
  ThetaDelay3_4$epoch <- "3-4"
  ThetaDelay4_5$epoch <- "4-5"
  ThetaDelay5_6$epoch <- "5-6"
  
  allThetaTrialData <- rbind(ThetaDelay1_2, ThetaDelay2_3)
  allThetaTrialData <- rbind(allThetaTrialData, ThetaDelay3_4)
  allThetaTrialData <- rbind(allThetaTrialData, ThetaDelay4_5)
  allThetaTrialData <- rbind(allThetaTrialData, ThetaDelay5_6)
  
  
  Avg_allThetaTrialData <- allThetaTrialData[] %>% dplyr::group_by(Subject,Trial) %>% summarize_all(mean)
  
  setnames(Avg_allThetaTrialData, c("Subject","Trial", "Trial_Power", "Event_Number","Event_Duration", "log1p_Trial_Power", "epoch"))
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(Avg_allThetaTrialData),value=TRUE) ###matches names
  
  names(Avg_allThetaTrialData)[names(Avg_allThetaTrialData) %in% analysisvars] <- paste("Theta",names(Avg_allThetaTrialData[,analysisvars]),sep=".") ##add Theta to analysis var names
  Thetavars <- grep("Theta",names(Avg_allThetaTrialData),value=TRUE)
  
  allThetaTrials <- Avg_allThetaTrialData[,c("Subject","Trial",Thetavars)]
  allThetaTrials <- allThetaTrials[complete.cases(allThetaTrials),]
  allThetaTrials <- allThetaTrials[!duplicated(allThetaTrials),]
  
  write.csv(allThetaTrials, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/delayTheta_Trial.csv')
  
  }

AverageFixationTrials <- function() {
  
 
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  GrandTable <- merge(GammaT, AlphaT, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
  GrandTable <- merge(GrandTable, BetaT, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
  GrandTable <- merge(GrandTable, ThetaT, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
  GrandTable <- merge(GrandTable, agefile, by = 'Subject', all.x = TRUE, all.y = TRUE)
  
  GrandTable$log1p_Gamma_Trial_Power <- log1p(GrandTable$Gamma_Trial_Power)
  GrandTable$log1p_Theta_Trial_Power <- log1p(GrandTable$Theta_Trial_Power)
  GrandTable$log1p_Beta_Trial_Power <- log1p(GrandTable$Beta_Trial_Power)
  GrandTable$log1p_Alpha_Trial_Power <- log1p(GrandTable$Alpha_Trial_Power)
  
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior.csv')
  
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior <- merge(Behavior, agefile, by = 'Subject', all.x = TRUE, all.y = TRUE)
  
  Everything <- merge(GrandTable, Behavior, by = c('Subject', 'Trial'), all.x = TRUE, all.y = TRUE)
  Everything$absPositionError <- abs(Everything$PositionError)
  
  ## average across trials 
  #gamma
  Avg_log1p_gamma_power <- aggregate(log1p_Gamma_Trial_Power ~ Subject, Everything[Everything$log1p_Gamma_Trial_Power != '-Inf',], mean)
  TrialVar_log1p_gamma_power <- aggregate(log1p_Gamma_Trial_Power ~ Subject, Everything[Everything$log1p_Gamma_Trial_Power != '-Inf',], sd)
  PerSubjectGamma <- merge(Avg_log1p_gamma_power, TrialVar_log1p_gamma_power, by = 'Subject', suffix = c("_mean", "_SD") )
  PerSubjectGamma <- merge(PerSubjectGamma, agefile, by= 'Subject')
  
  Avg_gamma_duration <- aggregate(Gamma_Event_Duration ~ Subject, Everything, mean)
  TrialVar_gamma_duration <- aggregate(Gamma_Event_Duration ~ Subject, Everything, sd)
  PerSubjectGamma <- merge(PerSubjectGamma, Avg_gamma_duration, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectGamma <- merge(PerSubjectGamma, TrialVar_gamma_duration, by= 'Subject', suffix = c("", "_SD"))
  
  Avg_gamma_event_number <- aggregate(Gamma_Event_Number ~ Subject, Everything, mean)
  TrialVar_gamma_event_number <- aggregate(Gamma_Event_Number ~ Subject, Everything, sd)
  PerSubjectGamma <- merge(PerSubjectGamma, Avg_gamma_event_number, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectGamma <- merge(PerSubjectGamma, TrialVar_gamma_event_number, by= 'Subject', suffix = c("", "_SD"))
  PerSubjectGamma$ageInverse = 1/PerSubjectGamma$age
  
  
  
  PerSubjectGamma <- setnames(PerSubjectGamma, c("Subject", "Trial_Power","Trial_Power_Variability", "idvalues", "age","visitno","Event_Duration", "Event_Duration_Variability", "Event_Number","Event_Number_Variability", "ageInverse")) 
  write.csv(PerSubjectGamma, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/fixGamma_20200810.csv')
  
  #alpha
  Avg_log1p_alpha_power <- aggregate(log1p_Alpha_Trial_Power ~ Subject, Everything[Everything$log1p_Alpha_Trial_Power != '-Inf',], mean)
  TrialVar_log1p_alpha_power <- aggregate(log1p_Alpha_Trial_Power ~ Subject, Everything[Everything$log1p_Alpha_Trial_Power != '-Inf',], sd)
  PerSubjectAlpha <- merge(Avg_log1p_alpha_power, TrialVar_log1p_alpha_power, by = 'Subject', suffix = c("_mean", "_SD") )
  PerSubjectAlpha <- merge(PerSubjectAlpha, agefile, by= 'Subject')
  
  Avg_alpha_duration <- aggregate(Alpha_Event_Duration ~ Subject, Everything, mean)
  TrialVar_alpha_duration <- aggregate(Alpha_Event_Duration ~ Subject, Everything, sd)
  PerSubjectAlpha <- merge(PerSubjectAlpha, Avg_alpha_duration, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectAlpha <- merge(PerSubjectAlpha, TrialVar_alpha_duration, by= 'Subject', suffix = c("", "_SD"))
  
  Avg_alpha_event_number <- aggregate(Alpha_Event_Number ~ Subject, Everything, mean)
  TrialVar_alpha_event_number <- aggregate(Alpha_Event_Number ~ Subject, Everything, sd)
  PerSubjectAlpha <- merge(PerSubjectAlpha, Avg_alpha_event_number, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectAlpha <- merge(PerSubjectAlpha, TrialVar_alpha_event_number, by= 'Subject', suffix = c("", "_SD"))
  PerSubjectAlpha$ageInverse = 1/PerSubjectAlpha$age
  
  PerSubjectAlpha <- setnames(PerSubjectAlpha, c("Subject", "Trial_Power","Trial_Power_Variability", "idvalues", "age","visitno", "Event_Duration", "Event_Duration_Variability", "Event_Number","Event_Number_Variability", "ageInverse")) 
  write.csv(PerSubjectAlpha, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/fixAlpha_20200810.csv')
  
  
  #theta
  Avg_log1p_theta_power <- aggregate(log1p_Theta_Trial_Power ~ Subject, Everything[Everything$log1p_Theta_Trial_Power != '-Inf',], mean)
  TrialVar_log1p_theta_power <- aggregate(log1p_Theta_Trial_Power ~ Subject, Everything[Everything$log1p_Theta_Trial_Power != '-Inf',], sd)
  PerSubjectTheta <- merge(Avg_log1p_theta_power, TrialVar_log1p_theta_power, by = 'Subject', suffix = c("_mean", "_SD") )
  PerSubjectTheta <- merge(PerSubjectTheta, agefile, by= 'Subject')
  
  Avg_theta_duration <- aggregate(Theta_Event_Duration ~ Subject, Everything, mean)
  TrialVar_theta_duration <- aggregate(Theta_Event_Duration ~ Subject, Everything, sd)
  PerSubjectTheta <- merge(PerSubjectTheta, Avg_theta_duration, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectTheta <- merge(PerSubjectTheta, TrialVar_theta_duration, by= 'Subject', suffix = c("", "_SD"))
  
  Avg_theta_event_number <- aggregate(Theta_Event_Number ~ Subject, Everything, mean)
  TrialVar_theta_event_number <- aggregate(Theta_Event_Number ~ Subject, Everything, sd)
  PerSubjectTheta <- merge(PerSubjectTheta, Avg_theta_event_number, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectTheta <- merge(PerSubjectTheta, TrialVar_theta_event_number, by= 'Subject', suffix = c("", "_SD"))
  PerSubjectTheta$ageInverse = 1/PerSubjectTheta$age
  
  PerSubjectTheta <- setnames(PerSubjectTheta, c("Subject", "Trial_Power","Trial_Power_Variability", "idvalues", "age","visitno","Event_Duration", "Event_Duration_Variability", "Event_Number","Event_Number_Variability", "ageInverse")) 
  write.csv(PerSubjectTheta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/fixTheta_20200810.csv')
  
  #beta
  Avg_log1p_beta_power <- aggregate(log1p_Beta_Trial_Power ~ Subject, Everything[Everything$log1p_Beta_Trial_Power != '-Inf',], mean)
  TrialVar_log1p_beta_power <- aggregate(log1p_Beta_Trial_Power ~ Subject, Everything[Everything$log1p_Beta_Trial_Power != '-Inf',], sd)
  PerSubjectBeta <- merge(Avg_log1p_beta_power, TrialVar_log1p_beta_power, by = 'Subject', suffix = c("_mean", "_SD") )
  PerSubjectBeta <- merge(PerSubjectBeta, agefile, by= 'Subject')
  
  Avg_beta_duration <- aggregate(Beta_Event_Duration ~ Subject, Everything, mean)
  TrialVar_beta_duration <- aggregate(Beta_Event_Duration ~ Subject, Everything, sd)
  PerSubjectBeta <- merge(PerSubjectBeta, Avg_beta_duration, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectBeta <- merge(PerSubjectBeta, TrialVar_beta_duration, by= 'Subject', suffix = c("", "_SD"))
  
  Avg_beta_event_number <- aggregate(Beta_Event_Number ~ Subject, Everything, mean)
  TrialVar_beta_event_number <- aggregate(Beta_Event_Number ~ Subject, Everything, sd)
  PerSubjectBeta <- merge(PerSubjectBeta, Avg_beta_event_number, by= 'Subject', suffix = c("", "_mean"))
  PerSubjectBeta <- merge(PerSubjectBeta, TrialVar_beta_event_number, by= 'Subject', suffix = c("", "_SD"))
  PerSubjectBeta$ageInverse = 1/PerSubjectBeta$age
  
  PerSubjectBeta <- setnames(PerSubjectBeta, c("Subject", "Trial_Power","Trial_Power_Variability", "idvalues", "age","visitno", "Event_Duration", "Event_Duration_Variability", "Event_Number","Event_Number_Variability", "ageInverse")) 
  write.csv(PerSubjectBeta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/fixBeta_20200810.csv')

  
  # peak frequency and power
  newSubsfixGammaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_newSubs_PeakFreq_PowerFix.csv')
  newSubsfixBetaPeaks  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_newSubs_PeakFreq_PowerFix.csv')
  newSubsfixAlphaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_newSubs_PeakFreq_Fix.csv')
  newSubsfixThetaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_newSubs_PeakFreq_PowerFix.csv')
  
  fixGammaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_all_PeakFreq_PowerFix.csv')
  fixBetaPeaks  <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_all_PeakFreq_PowerFix.csv')
  fixAlphaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_all_PeakFreq_Fix.csv')
  fixThetaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_all_PeakFreq_PowerFix.csv')
  
  
  fixGammaPeaks$log1p_Gamma_Trial_Power <- log1p(fixGammaPeaks$Gamma_Peak_Power)
  fixBetaPeaks$log1p_Beta_Trial_Power <- log1p(fixBetaPeaks$Beta_Peak_Power)
  fixAlphaPeaks$log1p_Alpha_Trial_Power <- log1p(fixAlphaPeaks$Alpha_Peak_Power)
  fixThetaPeaks$log1p_Theta_Trial_Power <- log1p(fixThetaPeaks$Theta_Peak_Power)
  
  PerSubjectGammaFixPeaks <- aggregate(.~Subject, fixGammaPeaks, mean)
  PerSubjectBetaFixPeaks <- aggregate(.~Subject, fixBetaPeaks, mean)
  PerSubjectAlphaFixPeaks <- aggregate(.~Subject, fixAlphaPeaks, mean)
  PerSubjectThetaFixPeaks <- aggregate(.~Subject, fixThetaPeaks, mean)
  
  write.csv(PerSubjectGammaFixPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/fixGamma_Peak_20200810.csv')
  write.csv(PerSubjectBetaFixPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/fixBeta_Peak_20200810.csv')
  write.csv(PerSubjectAlphaFixPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/fixAlpha_Peak_20200810.csv')
  write.csv(PerSubjectThetaFixPeaks, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/fixTheta_Peak_20200810.csv')
  
}

Read_Organize_Delay_Fix_Data <- function(){
  
  library(knitr)
  
  alldf <- list()
  
  delayGamma <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/delayGamma.csv')
  delayAlpha <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/delayAlpha.csv')
  delayTheta <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/delayTheta.csv')
  delayBeta <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/delayBeta.csv')
  
  delayGammaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/delayGamma_Peaks_20200810.csv')
  delayAlphaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/delayAlpha_Peaks_20200810.csv')
  delayThetaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/delayTheta_Peaks_20200810.csv')
  delayBetaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/delayBeta_Peaks_20200810.csv')
  
  
  all_delayGamma <- merge(delayGamma, delayGammaPeaks, by= "Subject", suffix = c("",""))
  all_delayAlpha <- merge(delayAlpha, delayAlphaPeaks,  by= "Subject", suffix = c("",""))
  all_delayTheta <- merge(delayTheta, delayThetaPeaks,  by= "Subject", suffix = c("",""))
  all_delayBeta <- merge(delayBeta, delayBetaPeaks,   by= "Subject", suffix = c("",""))
  
  write.csv(all_delayGamma, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/delayGamma_allVariables_20200810.csv')
  write.csv(all_delayAlpha, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/delayAlpha_allVariables_20200810.csv')
  write.csv(all_delayTheta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/delayTheta_allVariables_20200810.csv')
  write.csv(all_delayBeta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/delayBeta_allVariables_20200810.csv')
  
  
  fixGamma <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/fixGamma_20200810.csv')
  fixAlpha <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/fixAlpha_20200810.csv')
  fixTheta <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/fixTheta_20200810.csv')
  fixBeta <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/fixBeta_20200810.csv')
  
  fixGammaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/fixGamma_Peak_20200810.csv')
  fixAlphaPeaks <-read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/fixAlpha_Peak_20200810.csv')
  fixThetaPeaks <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/fixTheta_Peak_20200810.csv')
  fixBetaPeaks <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/fixBeta_Peak_20200810.csv')
  
  fixGamma <- merge(fixGamma, fixGammaPeaks, by = "Subject", suffix = c("",""))
  fixAlpha <- merge(fixAlpha, fixAlphaPeaks,by = "Subject", suffix = c("",""))
  fixTheta <- merge(fixTheta, fixThetaPeaks,by = "Subject", suffix = c("",""))
  fixBeta <- merge(fixBeta, fixBetaPeaks,by = "Subject", suffix = c("",""))
  
  write.csv(fixGamma, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/fixGamma_allVariables_20200810.csv')
  write.csv(fixAlpha, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/fixAlpha_allVariables_20200810.csv')
  write.csv(fixTheta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/fixTheta_allVariables_20200810.csv')
  write.csv(fixBeta, 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/fixBeta_allVariables_20200810.csv')
  
  
  #combine fixation and delay 
  allGamma <- merge(all_delayGamma[,c(1,3,4,5,6,7,8, 10, 11)], fixGamma[, c(1,3,4,5,6,7,8, 9,10, 11, 12, 14, 15, 16, 17)], by = "Subject", suffix = c("_Delay", "_Fix"))
  allAlpha <- merge(all_delayAlpha[,c(1,3,4,5,6,7,8, 10, 11)], fixAlpha[, c(1,3,4,5,6,7,8, 9,10, 11, 12, 14, 15, 16, 17)], by = "Subject", suffix = c("_Delay", "_Fix"))
  allTheta <- merge(all_delayTheta[,c(1,3,4,5,6,7,8, 10, 11)], fixTheta[, c(1,3,4,5,6,7,8, 9,10, 11, 12, 14, 15, 16, 17)], by = "Subject", suffix = c("_Delay", "_Fix"))
  allBeta <- merge(all_delayBeta[,c(1,3,4,5,6,7,8, 10, 11)], fixBeta[, c(1,3,4,5,6,7,8, 9,10, 11, 12, 14, 15, 16, 17)], by = "Subject", suffix = c("_Delay", "_Fix"))
  
  write.csv(allGamma,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allVariables_20200810.csv' )
  write.csv(allAlpha,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allVariables_20200810.csv' )
  write.csv(allTheta,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allVariables_20200810.csv' )
  write.csv(allBeta,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allVariables_20200810.csv' )
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior.csv')
  
  Behavior$log_absPositionError <- log1p(abs(Behavior$PositionError))
  Behavior$log_mgsLatency <- log1p(Behavior$mgsLatency)
  behavior <- setnames(aggregate(list(Behavior$log_mgsLatency, Behavior$log_absPositionError), list(Subject =  paste0(Behavior$LunaID, '_', Behavior$ScanDate)), na.rm = TRUE, mean), c("Subject", "log_mgsLatency", "log_absPositionError"))
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allVariables_20200810.csv')
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta)[names(theta) %in% analysisvars] <- paste("Theta",names(theta[,analysisvars]),sep=".")
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  alldf$allTheta <- allTheta

  #alpha frequency band 
  alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allVariables_20200810.csv')
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha)[names(alpha) %in% analysisvars] <- paste("Alpha",names(alpha[,analysisvars]),sep=".")
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  alldf$allAlpha <- allAlpha

  #beta frequency band 
  beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allVariables_20200810.csv' )
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta)[names(beta) %in% analysisvars] <- paste("Beta",names(beta[,analysisvars]),sep=".")
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  alldf$allBeta <- allBeta

  #gamma frequency band 
  gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allVariables_20200810.csv')
  analysisvars <- grep("Trial|Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma)[names(gamma) %in% analysisvars] <- paste("Gamma",names(gamma[,analysisvars]),sep=".") ##add Gamma to analysis var names
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  alldf$allGamma <- allGamma

  ## ---- alldata ----
  alldata <- merge(alpha[,c("Subject",alphavars)],beta[,c("Subject",betavars)],by="Subject") %>% 
    merge(.,theta[,c("Subject",thetavars)],by="Subject") %>% merge(.,gamma[,c("Subject",gammavars)],by="Subject")
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  
  alldata <- merge(alldata, agefile, by = "Subject")
  
  
  alldata <- alldata[complete.cases(alldata),]
  alldata <- alldata[!duplicated(alldata),]
  alldata$inverseAge <- 1/alldata$age
  alldf$alldata <- alldata
  
  return(alldf)
  
  
}


Only_Take_Fix_TrialLevel <- function () {
  
  alldf_Fix <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Subs_dataFix.csv')
 
  alldf_Fix$allTheta <- theta
  
  #alpha frequency band 
  alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Subs_datafix.csv')

    alldf_Fix$allAlpha <- alpha
  
  #beta frequency band 
  beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Subs_dataFix.csv' )

  alldf_Fix$allBeta <- beta
  
  #gamma frequency band 
  gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Subs_FIX.csv')
 
  alldf_Fix$allGamma <- gamma
  
  ## ---- alldata ----
  alldata_Fix_TrialLevel <- merge(alpha[],beta[],by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% 
    merge(.,theta[],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,gamma[],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  agefile$inverseAge <- 1/agefile$age
  
  alldata_Fix_TrialLevel <- merge(alldata_Fix_TrialLevel, agefile, by = "Subject",  all.x = TRUE, all.y = TRUE)
  
  alldf_Fix$alldata_Fix_TrialLevel <- alldata_Fix_TrialLevel
  
  return(alldf_Fix)
  
  
}

Resting_State_Data_TrialLevel <- function () {
  
  alldf_RS <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_trialLevel_data_RS.csv')
  theta$log1p_Trial_Power <- log1p(theta$Theta_Trial_Power)
  theta <- setnames(theta, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(theta),value=TRUE) 
  names(theta)[names(theta) %in% analysisvars] <- paste("Theta",names(theta[,analysisvars]),sep=".")
  
  theta <- theta[complete.cases(theta),]
  theta <- theta[!duplicated(theta),]
  alldf_RS$theta <- theta
  
  #alpha frequency band 
  alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_trialLevel_data_RS.csv')
  alpha$log1p_Trial_Power <- log1p(alpha$Alpha_Trial_Power)
  alpha <- setnames(alpha, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(alpha),value=TRUE) 
  
  names(alpha)[names(alpha) %in% analysisvars] <- paste("Alpha",names(alpha[,analysisvars]),sep=".")
  
  alpha <- alpha[complete.cases(alpha),]
  alpha <- alpha[!duplicated(alpha),]
  alldf_RS$alpha <- alpha
  
  #beta frequency band 
  beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_trialLevel_data_RS.csv' )
  beta$log1p_Trial_Power <- log1p(beta$Beta_Trial_Power)
  beta <- setnames(beta, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(beta),value=TRUE) 
  
  names(beta)[names(beta) %in% analysisvars] <- paste("Beta",names(beta[,analysisvars]),sep=".")
  beta <- beta[complete.cases(beta),]
  beta <- beta[!duplicated(beta),]
  alldf_RS$beta <- beta
  
  #gamma frequency band 
  gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_trialLevel_data_RS.csv')
  gamma$log1p_Trial_Power <- log1p(gamma$Gamma_Trial_Power)
  gamma <- setnames(gamma, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(gamma),value=TRUE) 
  
  names(gamma)[names(gamma) %in% analysisvars] <- paste("Gamma",names(gamma[,analysisvars]),sep=".")
  gamma <- gamma[complete.cases(gamma),]
  gamma <- gamma[!duplicated(gamma),]
  alldf_RS$gamma <- gamma
  
  ## ---- alldata ----
  alldata_TrialLevel <- merge(alpha,beta,by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% 
    merge(.,theta, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,gamma, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  agefile$inverseAge <- 1/agefile$age
  
  alldata_TrialLevel <- merge(alldata_TrialLevel, agefile, by = "Subject",  all.x = TRUE, all.y = TRUE)
  
  alldf_RS$alldata_RS_TrialLevel <- alldata_TrialLevel
  
  return(alldf_RS)
  
  
  
  
  
  
}

Resting_State_Data_SubjectLevel <- function () {
  
  alldf_SubjectLevel_RS <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_trialLevel_data_RS.csv')
  theta$log1p_Trial_Power <- log1p(theta$Theta_Trial_Power)
  theta <- setnames(theta, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(theta),value=TRUE) 
  names(theta)[names(theta) %in% analysisvars] <- paste("Theta",names(theta[,analysisvars]),sep=".")
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  theta <- theta[complete.cases(theta),]
  theta <- theta[!duplicated(theta),]
  
  firstDF <- theta[seq(1, nrow(theta),by = 3),]
  firstDF_sd <- aggregate(.~Subject, firstDF, sd)
  
  secondDF <- theta[seq(2, nrow(theta),by = 3),]
  secondDF_sd <- aggregate(.~Subject, secondDF, sd)
  
  thirdDF <- theta[seq(3, nrow(theta),by = 3),]
  thirdDF_sd <- aggregate(.~Subject, thirdDF, sd)
  
  
  Theta.log1p_Trial_Power <- (firstDF_sd$Theta.log1p_Trial_Power + secondDF_sd$Theta.log1p_Trial_Power + thirdDF_sd$Theta.log1p_Trial_Power)/3
  Theta.Event_Number <- (firstDF_sd$Theta.Event_Number + secondDF_sd$Theta.Event_Number + thirdDF_sd$Theta.Event_Number)/3
  Theta.Event_Duration <- (firstDF_sd$Theta.Event_Duration + secondDF_sd$Theta.Event_Duration + thirdDF_sd$Theta.Event_Duration)/3
  Subject <- firstDF_sd$Subject
  
  avg_SDs <- data.frame(Subject, Theta.log1p_Trial_Power, Theta.Event_Number, Theta.Event_Duration)
 
  theta_mean <- aggregate(.~Subject, theta, mean)

  allTheta <- merge(theta_mean, avg_SDs, by = "Subject", suffixes = c("","_Variability"))
  
  alldf_SubjectLevel_RS$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_trialLevel_data_RS.csv')
  alpha$log1p_Trial_Power <- log1p(alpha$Alpha_Trial_Power)
  alpha <- setnames(alpha, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(alpha),value=TRUE) 
  
  names(alpha)[names(alpha) %in% analysisvars] <- paste("Alpha",names(alpha[,analysisvars]),sep=".")
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  alpha <- alpha[complete.cases(alpha),]
  alpha <- alpha[!duplicated(alpha),]
  
  firstDF <- alpha[seq(1, nrow(alpha),by = 3),]
  firstDF_sd <- aggregate(.~Subject, firstDF, sd)
  
  secondDF <- alpha[seq(2, nrow(alpha),by = 3),]
  secondDF_sd <- aggregate(.~Subject, secondDF, sd)
  
  thirdDF <- alpha[seq(3, nrow(alpha),by = 3),]
  thirdDF_sd <- aggregate(.~Subject, thirdDF, sd)
  
  
  Alpha.log1p_Trial_Power <- (firstDF_sd$Alpha.log1p_Trial_Power + secondDF_sd$Alpha.log1p_Trial_Power + thirdDF_sd$Alpha.log1p_Trial_Power)/3
  Alpha.Event_Number <- (firstDF_sd$Alpha.Event_Number + secondDF_sd$Alpha.Event_Number + thirdDF_sd$Alpha.Event_Number)/3
  Alpha.Event_Duration <- (firstDF_sd$Alpha.Event_Duration + secondDF_sd$Alpha.Event_Duration + thirdDF_sd$Alpha.Event_Duration)/3
  Subject <- firstDF_sd$Subject
  
  avg_SDs <- data.frame(Subject, Alpha.log1p_Trial_Power, Alpha.Event_Number, Alpha.Event_Duration)
  
  alpha_mean <- aggregate(.~Subject, alpha, mean)
  
  allAlpha <- merge(alpha_mean, avg_SDs, by = "Subject", suffixes = c("","_Variability"))
  
  
  alldf_SubjectLevel_RS$allAlpha <- allAlpha

  #beta frequency band 
  beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_trialLevel_data_RS.csv' )
  beta$log1p_Trial_Power <- log1p(beta$Beta_Trial_Power)
  beta <- setnames(beta, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(beta),value=TRUE) 
  
  names(beta)[names(beta) %in% analysisvars] <- paste("Beta",names(beta[,analysisvars]),sep=".")
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  beta <- beta[complete.cases(beta),]
  beta <- beta[!duplicated(beta),]
  
  firstDF <- beta[seq(1, nrow(beta),by = 3),]
  firstDF_sd <- aggregate(.~Subject, firstDF, sd)
  
  secondDF <- beta[seq(2, nrow(beta),by = 3),]
  secondDF_sd <- aggregate(.~Subject, secondDF, sd)
  
  thirdDF <- beta[seq(3, nrow(beta),by = 3),]
  thirdDF_sd <- aggregate(.~Subject, thirdDF, sd)
  
  
  Beta.log1p_Trial_Power <- (firstDF_sd$Beta.log1p_Trial_Power + secondDF_sd$Beta.log1p_Trial_Power + thirdDF_sd$Beta.log1p_Trial_Power)/3
  Beta.Event_Number <- (firstDF_sd$Beta.Event_Number + secondDF_sd$Beta.Event_Number + thirdDF_sd$Beta.Event_Number)/3
  Beta.Event_Duration <- (firstDF_sd$Beta.Event_Duration + secondDF_sd$Beta.Event_Duration + thirdDF_sd$Beta.Event_Duration)/3
  Subject <- firstDF_sd$Subject
  
  avg_SDs <- data.frame(Subject, Beta.log1p_Trial_Power, Beta.Event_Number, Beta.Event_Duration)
  
  Beta_mean <- aggregate(.~Subject, beta, mean)
  
  allBeta <- merge(Beta_mean, avg_SDs, by = "Subject", suffixes = c("","_Variability"))
  
  alldf_SubjectLevel_RS$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_trialLevel_data_RS.csv')
  gamma$log1p_Trial_Power <- log1p(gamma$Gamma_Trial_Power)
  gamma <- setnames(gamma, c("Subject", "Trial", "Trial_Power", "Event_Number", "Event_Duration", "log1p_Trial_Power") )
  analysisvars <- grep("Power|Duration|Number",names(gamma),value=TRUE) 
  
  names(gamma)[names(gamma) %in% analysisvars] <- paste("Gamma",names(gamma[,analysisvars]),sep=".")
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  gamma <- gamma[complete.cases(gamma),]
  gamma <- gamma[!duplicated(gamma),]
  
  firstDF <- gamma[seq(1, nrow(gamma),by = 3),]
  firstDF_sd <- aggregate(.~Subject, firstDF, sd)
  
  secondDF <- gamma[seq(2, nrow(gamma),by = 3),]
  secondDF_sd <- aggregate(.~Subject, secondDF, sd)
  
  thirdDF <- gamma[seq(3, nrow(gamma),by = 3),]
  thirdDF_sd <- aggregate(.~Subject, thirdDF, sd)
  
  
  Gamma.log1p_Trial_Power <- (firstDF_sd$Gamma.log1p_Trial_Power + secondDF_sd$Gamma.log1p_Trial_Power + thirdDF_sd$Gamma.log1p_Trial_Power)/3
  Gamma.Event_Number <- (firstDF_sd$Gamma.Event_Number + secondDF_sd$Gamma.Event_Number + thirdDF_sd$Gamma.Event_Number)/3
  Gamma.Event_Duration <- (firstDF_sd$Gamma.Event_Duration + secondDF_sd$Gamma.Event_Duration + thirdDF_sd$Gamma.Event_Duration)/3
  Subject <- firstDF_sd$Subject
  
  avg_SDs <- data.frame(Subject, Gamma.log1p_Trial_Power, Gamma.Event_Number, Gamma.Event_Duration)
  
  Gamma_mean <- aggregate(.~Subject, gamma, mean)
  
  allGamma <- merge(Gamma_mean, avg_SDs, by = "Subject", suffixes = c("","_Variability"))
  
  alldf_SubjectLevel_RS$allGamma <- allGamma
  
  ## ---- alldata ----
  rest_SubLevel <- merge(allAlpha,allBeta,by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% 
    merge(.,allTheta, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,allGamma, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  agefile$inverseAge <- 1/agefile$age
  
  rest_SubLevel <- merge(rest_SubLevel, agefile, by = "Subject")
  rest_SubLevel <- rest_SubLevel[rest_SubLevel$visitno <2,]
  
  alldf_SubjectLevel_RS$rest_SubLevel <- rest_SubLevel
  
  return(alldf_SubjectLevel_RS)
  
  
  
  
  
  
}

RestDelay_IndividualChannels_TrialLevel <- function () {

agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')

channelLocations <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')

# Gamma

GammaResting <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Resting.csv')  

GammaDelay <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 

colnames(GammaResting)[1] <- "idvalues"
colnames(GammaDelay)[1] <- "idvalues"

GammaResting_Age <- merge(GammaResting, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
GammaResting_Age$log_Gamma_Power <- log1p(GammaResting_Age$Gamma_Trial_Power)
GammaResting_Age$inverseAge <- 1/GammaResting_Age$age
GammaResting_Age_Channel <- merge(GammaResting_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
GammaResting_Age_Channel$Task <- 'Rest'

GammaDelay_Age <- merge(GammaDelay, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
GammaDelay_Age$log_Gamma_Power <- log1p(GammaDelay_Age$Gamma_Trial_Power)
GammaDelay_Age$inverseAge <- 1/GammaDelay_Age$age
GammaDelay_Age_Channel <- merge(GammaDelay_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
GammaDelay_Age_Channel$Task <- 'Delay'

Gamma_DelayRest <- rbind(GammaResting_Age_Channel, GammaDelay_Age_Channel)

}

DelayOnly_IndividualChannels_TrialLevel <- function () {
  
  individualChannelDF <- list()
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv')
  
  channelLocations <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')
  
  # Gamma

  GammaDelay <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
  

  GammaDelay_Age <- merge(GammaDelay, agefile, by = "Subject")
  GammaDelay_Age$log_Gamma_Power <- log1p(GammaDelay_Age$Gamma_Trial_Power)
  GammaDelay_Age$inverseAge <- 1/GammaDelay_Age$age
  GammaDelay_Age_Channel <- merge(GammaDelay_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  GammaDelay_Age_Channel$Task <- 'Delay'
  
  GammaDelay_Age_Channel <- subset(GammaDelay_Age_Channel, visitno == 1) 
    
  # outlier detection 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  GammaDelay_Age_Channel_new <- GammaDelay_Age_Channel
  
  cols = names(GammaDelay_Age_Channel[c(4:6,10)])
  for (col in cols) {
    
    GammaDelay_Age_Channel_group <- GammaDelay_Age_Channel %>% group_by(Subject) %>% group_by(Channel)
      indx <- outliers(GammaDelay_Age_Channel_group[[col]])
    
    GammaDelay_Age_Channel_new[[col]] <- Map(replace, GammaDelay_Age_Channel_new[[col]], indx, NA)
    GammaDelay_Age_Channel_new[[col]] <- as.numeric(GammaDelay_Age_Channel_new[[col]])
    
  }  
  
  individualChannelDF$GammaDelay_Age_Channel_new <- GammaDelay_Age_Channel_new
  
  
  # Beta
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  

  Beta_Age <- merge(Beta, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  Beta_Age$log_Beta_Power <- log1p(Beta_Age$Beta_Trial_Power)
  Beta_Age$inverseAge <- 1/Beta_Age$age
  
  Beta_Age_Channel <- merge(Beta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  Beta_Age_Channel <- subset(Beta_Age_Channel, visitno == 1) 
  
  
  Beta_Age_Channel_new <- Beta_Age_Channel
  
  cols = names(Beta_Age_Channel[c(4:6,10)])
  for (col in cols) {
    
    Beta_Age_Channel_group <- Beta_Age_Channel %>% group_by(Subject) %>% group_by(Channel)
    indx <- outliers(Beta_Age_Channel_group[[col]])
    
    Beta_Age_Channel_new[[col]] <- Map(replace, Beta_Age_Channel_new[[col]], indx, NA)
    Beta_Age_Channel_new[[col]] <- as.numeric(Beta_Age_Channel_new[[col]])
    
  }  
  
  individualChannelDF$Beta_Age_Channel_new <- Beta_Age_Channel_new
  
  # Alpha
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  
  Alpha_Age <- merge(Alpha, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  Alpha_Age$log_Alpha_Power <- log1p(Alpha_Age$Alpha_Trial_Power)
  Alpha_Age$inverseAge <- 1/Alpha_Age$age
  
  Alpha_Age_Channel <- merge(Alpha_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  Alpha_Age_Channel <- subset(Alpha_Age_Channel, visitno == 1) 
  
  
  Alpha_Age_Channel_new <- Alpha_Age_Channel
  
  cols = names(Alpha_Age_Channel[c(4:6,10)])
  for (col in cols) {
    
    Alpha_Age_Channel_group <- Alpha_Age_Channel %>% group_by(Subject) %>% group_by(Channel)
    indx <- outliers(Alpha_Age_Channel_group[[col]])
    
    Alpha_Age_Channel_new[[col]] <- Map(replace, Alpha_Age_Channel_new[[col]], indx, NA)
    Alpha_Age_Channel_new[[col]] <- as.numeric(Alpha_Age_Channel_new[[col]])
    
  }  
  
  individualChannelDF$Alpha_Age_Channel_new <- Alpha_Age_Channel_new
  
  # Theta
  Theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  
  Theta_Age <- merge(Theta, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  Theta_Age$log_Theta_Power <- log1p(Theta_Age$Theta_Trial_Power)
  Theta_Age$inverseAge <- 1/Theta_Age$age
  
  Theta_Age_Channel <- merge(Theta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  Theta_Age_Channel <- subset(Theta_Age_Channel, visitno == 1) 
  
  
  Theta_Age_Channel_new <- Theta_Age_Channel
  
  cols = names(Theta_Age_Channel[c(4:6,10)])
  for (col in cols) {
    
    Theta_Age_Channel_group <- Theta_Age_Channel %>% group_by(Subject) %>% group_by(Channel)
    indx <- outliers(Theta_Age_Channel_group[[col]])
    
    Theta_Age_Channel_new[[col]] <- Map(replace, Theta_Age_Channel_new[[col]], indx, NA)
    Theta_Age_Channel_new[[col]] <- as.numeric(Theta_Age_Channel_new[[col]])
    
  }  
  
  individualChannelDF$Theta_Age_Channel_new <- Theta_Age_Channel_new
  
  return(individualChannelDF)
  

}

IndividualChannels_GroupedintoRegions <- function () {
  
  alldf <- list()
  #DLPFC
  #Gamma
  
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  Gamma_DLPFC <- filter(GammaDelay_Age_Channel, GammaDelay_Age_Channel$Label == "'F3'" | GammaDelay_Age_Channel$Label == "'F4'")
  
  #keep the data on the trial level, but aggregate all the DLPFC channels together
  Gamma_DLPFC_avgChannel_Events <- Gamma_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  Gamma_DLPFC_avgChannel_logPower <- Gamma_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  Gamma_DLPFC_avgChannel_Power <- Gamma_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = mean(Gamma_Trial_Power, na.rm = T))
  
  Gamma_DLPFC_avgChannel_Duration <- Gamma_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_DLPFC_avgChannel <- merge(Gamma_DLPFC_avgChannel_Events, Gamma_DLPFC_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_DLPFC_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Gamma_DLPFC_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  #variability
  Gamma_DLPFC_sdChannel_Events <- Gamma_DLPFC_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  Gamma_DLPFC_sdChannel_logPower <- Gamma_DLPFC_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  Gamma_DLPFC_sdChannel_Power <- Gamma_DLPFC_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = sd(Gamma_Trial_Power, na.rm = T))
  Gamma_DLPFC_sdChannel_Duration <- Gamma_DLPFC_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_DLPFC_sdChannel <- merge(Gamma_DLPFC_sdChannel_Events, Gamma_DLPFC_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_DLPFC_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Gamma_DLPFC_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  #Beta
  
  BetaDelay_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
  Beta_DLPFC <- filter(BetaDelay_Age_Channel, BetaDelay_Age_Channel$Label == "'F3'" | BetaDelay_Age_Channel$Label == "'F4'")
  
  #keep the data on the trial level, but aggregate all the DLPFC channels together
  Beta_DLPFC_avgChannel_Events <- Beta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))
  Beta_DLPFC_avgChannel_logPower <- Beta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))
  Beta_DLPFC_avgChannel_Power <- Beta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = mean(Beta_Trial_Power, na.rm = T))
  
  Beta_DLPFC_avgChannel_Duration <- Beta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))
  
  Beta_DLPFC_avgChannel <- merge(Beta_DLPFC_avgChannel_Events, Beta_DLPFC_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_DLPFC_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Beta_DLPFC_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  #variability
  Beta_DLPFC_sdChannel_Events <- Beta_DLPFC_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))
  Beta_DLPFC_sdChannel_logPower <- Beta_DLPFC_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))
  Beta_DLPFC_sdChannel_Power <- Beta_DLPFC_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = sd(Beta_Trial_Power, na.rm = T))
  Beta_DLPFC_sdChannel_Duration <- Beta_DLPFC_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))
  
  Beta_DLPFC_sdChannel <- merge(Beta_DLPFC_sdChannel_Events, Beta_DLPFC_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_DLPFC_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Beta_DLPFC_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  #Alpha
  
  AlphaDelay_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
  Alpha_DLPFC <- filter(AlphaDelay_Age_Channel, AlphaDelay_Age_Channel$Label == "'F3'" | AlphaDelay_Age_Channel$Label == "'F4'")
  
  #keep the data on the trial level, but aggregate all the DLPFC channels together
  Alpha_DLPFC_avgChannel_Events <- Alpha_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
  Alpha_DLPFC_avgChannel_logPower <- Alpha_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))
  Alpha_DLPFC_avgChannel_Power <- Alpha_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = mean(Alpha_Trial_Power, na.rm = T))
  
  Alpha_DLPFC_avgChannel_Duration <- Alpha_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
  
  Alpha_DLPFC_avgChannel <- merge(Alpha_DLPFC_avgChannel_Events, Alpha_DLPFC_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_DLPFC_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Alpha_DLPFC_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  #variability
  Alpha_DLPFC_sdChannel_Events <- Alpha_DLPFC_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))
  Alpha_DLPFC_sdChannel_logPower <- Alpha_DLPFC_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))
  Alpha_DLPFC_sdChannel_Power <- Alpha_DLPFC_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = sd(Alpha_Trial_Power, na.rm = T))
  Alpha_DLPFC_sdChannel_Duration <- Alpha_DLPFC_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))
  
  Alpha_DLPFC_sdChannel <- merge(Alpha_DLPFC_sdChannel_Events, Alpha_DLPFC_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_DLPFC_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Alpha_DLPFC_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  #Theta
  
  ThetaDelay_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
  Theta_DLPFC <- filter(ThetaDelay_Age_Channel, ThetaDelay_Age_Channel$Label == "'F3'" | ThetaDelay_Age_Channel$Label == "'F4'")
  
  #keep the data on the trial level, but aggregate all the DLPFC channels together
  Theta_DLPFC_avgChannel_Events <- Theta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))
  Theta_DLPFC_avgChannel_logPower <- Theta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))
  Theta_DLPFC_avgChannel_Power <- Theta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = mean(Theta_Trial_Power, na.rm = T))
  
  Theta_DLPFC_avgChannel_Duration <- Theta_DLPFC %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))
  
  Theta_DLPFC_avgChannel <- merge(Theta_DLPFC_avgChannel_Events, Theta_DLPFC_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_DLPFC_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Theta_DLPFC_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  #variability
  Theta_DLPFC_sdChannel_Events <- Theta_DLPFC_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))
  Theta_DLPFC_sdChannel_logPower <- Theta_DLPFC_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))
  Theta_DLPFC_sdChannel_Power <- Theta_DLPFC_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = sd(Theta_Trial_Power, na.rm = T))
  Theta_DLPFC_sdChannel_Duration <- Theta_DLPFC_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))
  
  Theta_DLPFC_sdChannel <- merge(Theta_DLPFC_sdChannel_Events, Theta_DLPFC_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_DLPFC_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Theta_DLPFC_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  
  allDLPFC<- merge(Gamma_DLPFC_avgChannel, Beta_DLPFC_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_DLPFC_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_DLPFC_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))
  
  allDLPFC_sd <- merge(Gamma_DLPFC_sdChannel, Beta_DLPFC_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_DLPFC_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_DLPFC_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
  
  alldf$allDLPFC <- allDLPFC
  alldf$allDLPFC_sd <- allDLPFC_sd
  
  
  #Frontal 
  #Gamma
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  Gamma_Frontal <- filter(GammaDelay_Age_Channel, str_detect(GammaDelay_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the frontal channels together
  Gamma_Frontal_avgChannel_Events <- Gamma_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  Gamma_Frontal_avgChannel_logPower <- Gamma_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  Gamma_Frontal_avgChannel_Power <- Gamma_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = mean(Gamma_Trial_Power, na.rm = T))
  
  Gamma_Frontal_avgChannel_Duration <- Gamma_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Frontal_avgChannel <- merge(Gamma_Frontal_avgChannel_Events, Gamma_Frontal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Gamma_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Gamma_Frontal_sdChannel_Events <- Gamma_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  Gamma_Frontal_sdChannel_logPower <- Gamma_Frontal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  Gamma_Frontal_sdChannel_Power <- Gamma_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = sd(Gamma_Trial_Power, na.rm = T))
  Gamma_Frontal_sdChannel_Duration <- Gamma_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Frontal_sdChannel <- merge(Gamma_Frontal_sdChannel_Events, Gamma_Frontal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Gamma_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  #Beta
  Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
  Beta_Frontal <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the frontal channels together
  Beta_Frontal_avgChannel_Events <- Beta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))
  Beta_Frontal_avgChannel_logPower <- Beta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))
  Beta_Frontal_avgChannel_Power <- Beta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = mean(Beta_Trial_Power, na.rm = T))
  Beta_Frontal_avgChannel_Duration <- Beta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))
  
  Beta_Frontal_avgChannel <- merge(Beta_Frontal_avgChannel_Events, Beta_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Beta_Frontal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Beta_Frontal_sdChannel_Events <- Beta_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))
  Beta_Frontal_sdChannel_logPower <- Beta_Frontal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))
  Beta_Frontal_sdChannel_Power <- Beta_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = sd(Beta_Trial_Power, na.rm = T))
  Beta_Frontal_sdChannel_Duration <- Beta_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))
  
  Beta_Frontal_sdChannel <- merge(Beta_Frontal_sdChannel_Events, Beta_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Beta_Frontal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Alpha
  Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
  Alpha_Frontal <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the frontal channels together
  Alpha_Frontal_avgChannel_Events <- Alpha_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
  Alpha_Frontal_avgChannel_logPower <- Alpha_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))
  Alpha_Frontal_avgChannel_Power <- Alpha_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = mean(Alpha_Trial_Power, na.rm = T))
  Alpha_Frontal_avgChannel_Duration <- Alpha_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Frontal_avgChannel <- merge(Alpha_Frontal_avgChannel_Events, Alpha_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Frontal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the frontal channels and THEN take SD of the trials 
  Alpha_Frontal_sdChannel_Events <- Alpha_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))
  Alpha_Frontal_sdChannel_logPower <- Alpha_Frontal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))
  Alpha_Frontal_sdChannel_Power <- Alpha_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = sd(Alpha_Trial_Power, na.rm = T))
  Alpha_Frontal_sdChannel_Duration <- Alpha_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Frontal_sdChannel <- merge(Alpha_Frontal_sdChannel_Events, Alpha_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Frontal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Theta
  Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
  Theta_Frontal <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the frontal channels together
  Theta_Frontal_avgChannel_Events <- Theta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))
  Theta_Frontal_avgChannel_logPower <- Theta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))
  Theta_Frontal_avgChannel_Power <- Theta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = mean(Theta_Trial_Power, na.rm = T))
  Theta_Frontal_avgChannel_Duration <- Theta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))
  
  Theta_Frontal_avgChannel <- merge(Theta_Frontal_avgChannel_Events, Theta_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Frontal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the frontal channels and THEN take SD of the trials 
  Theta_Frontal_sdChannel_Events <- Theta_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))
  Theta_Frontal_sdChannel_logPower <- Theta_Frontal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))
  Theta_Frontal_sdChannel_Power <- Theta_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = sd(Theta_Trial_Power, na.rm = T))
  Theta_Frontal_sdChannel_Duration <- Theta_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))
  
  Theta_Frontal_sdChannel <- merge(Theta_Frontal_sdChannel_Events, Theta_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Theta_Frontal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  allFrontal <- merge(Gamma_Frontal_avgChannel, Beta_Frontal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_Frontal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_Frontal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))
  
  allFrontal_sd <- merge(Gamma_Frontal_sdChannel, Beta_Frontal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_Frontal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_Frontal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
  
  alldf$allFrontal <- allFrontal
  alldf$allFrontal_sd <- allFrontal_sd
  
  
  
  
  #Parietal 
  #Gamma
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  Gamma_Parietal <- filter(GammaDelay_Age_Channel, str_detect(GammaDelay_Age_Channel$Label, "P"))
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Gamma_Parietal_avgChannel_Events <- Gamma_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  Gamma_Parietal_avgChannel_logPower <- Gamma_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  Gamma_Parietal_avgChannel_Power <- Gamma_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = mean(Gamma_Trial_Power, na.rm = T))
  
  Gamma_Parietal_avgChannel_Duration <- Gamma_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Parietal_avgChannel <- merge(Gamma_Parietal_avgChannel_Events, Gamma_Parietal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Gamma_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Gamma_Parietal_sdChannel_Events <- Gamma_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  Gamma_Parietal_sdChannel_logPower <- Gamma_Parietal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  Gamma_Parietal_sdChannel_Power <- Gamma_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = sd(Gamma_Trial_Power, na.rm = T))
  Gamma_Parietal_sdChannel_Duration <- Gamma_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Parietal_sdChannel <- merge(Gamma_Parietal_sdChannel_Events, Gamma_Parietal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Gamma_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  #Beta
  Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
  Beta_Parietal <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "P"))
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Beta_Parietal_avgChannel_Events <- Beta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))
  Beta_Parietal_avgChannel_logPower <- Beta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))
  Beta_Parietal_avgChannel_Power <- Beta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = mean(Beta_Trial_Power, na.rm = T))
  Beta_Parietal_avgChannel_Duration <- Beta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))
  
  Beta_Parietal_avgChannel <- merge(Beta_Parietal_avgChannel_Events, Beta_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Beta_Parietal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Beta_Parietal_sdChannel_Events <- Beta_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))
  Beta_Parietal_sdChannel_logPower <- Beta_Parietal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))
  Beta_Parietal_sdChannel_Power <- Beta_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = sd(Beta_Trial_Power, na.rm = T))
  Beta_Parietal_sdChannel_Duration <- Beta_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))
  
  Beta_Parietal_sdChannel <- merge(Beta_Parietal_sdChannel_Events, Beta_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Beta_Parietal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Alpha
  Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
  Alpha_Parietal <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "P"))
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Alpha_Parietal_avgChannel_Events <- Alpha_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
  Alpha_Parietal_avgChannel_logPower <- Alpha_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))
  Alpha_Parietal_avgChannel_Power <- Alpha_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = mean(Alpha_Trial_Power, na.rm = T))
  Alpha_Parietal_avgChannel_Duration <- Alpha_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Parietal_avgChannel <- merge(Alpha_Parietal_avgChannel_Events, Alpha_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Parietal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Parietal channels and THEN take SD of the trials 
  Alpha_Parietal_sdChannel_Events <- Alpha_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))
  Alpha_Parietal_sdChannel_logPower <- Alpha_Parietal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))
  Alpha_Parietal_sdChannel_Power <- Alpha_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = sd(Alpha_Trial_Power, na.rm = T))
  Alpha_Parietal_sdChannel_Duration <- Alpha_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Parietal_sdChannel <- merge(Alpha_Parietal_sdChannel_Events, Alpha_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Parietal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Theta
  Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
  Theta_Parietal <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "P"))
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Theta_Parietal_avgChannel_Events <- Theta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))
  Theta_Parietal_avgChannel_logPower <- Theta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))
  Theta_Parietal_avgChannel_Power <- Theta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = mean(Theta_Trial_Power, na.rm = T))
  Theta_Parietal_avgChannel_Duration <- Theta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))
  
  Theta_Parietal_avgChannel <- merge(Theta_Parietal_avgChannel_Events, Theta_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Parietal_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Parietal channels and THEN take SD of the trials 
  Theta_Parietal_sdChannel_Events <- Theta_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))
  Theta_Parietal_sdChannel_logPower <- Theta_Parietal_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))
  Theta_Parietal_sdChannel_Power <- Theta_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = sd(Theta_Trial_Power, na.rm = T))
  Theta_Parietal_sdChannel_Duration <- Theta_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))
  
  Theta_Parietal_sdChannel <- merge(Theta_Parietal_sdChannel_Events, Theta_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Theta_Parietal_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  allParietal <- merge(Gamma_Parietal_avgChannel, Beta_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))
  
  allParietal_sd <- merge(Gamma_Parietal_sdChannel, Beta_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
  
  alldf$allParietal <- allParietal
  alldf$allParietal_sd <- allParietal_sd
  
  
  #Occipital 
  #Gamma
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  Gamma_Occipital <- filter(GammaDelay_Age_Channel, str_detect(GammaDelay_Age_Channel$Label, "O"))
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Gamma_Occipital_avgChannel_Events <- Gamma_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  Gamma_Occipital_avgChannel_logPower <- Gamma_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  Gamma_Occipital_avgChannel_Power <- Gamma_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = mean(Gamma_Trial_Power, na.rm = T))
  
  Gamma_Occipital_avgChannel_Duration <- Gamma_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Occipital_avgChannel <- merge(Gamma_Occipital_avgChannel_Events, Gamma_Occipital_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Gamma_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Gamma_Occipital_sdChannel_Events <- Gamma_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  Gamma_Occipital_sdChannel_logPower <- Gamma_Occipital_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  Gamma_Occipital_sdChannel_Power <- Gamma_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = sd(Gamma_Trial_Power, na.rm = T))
  Gamma_Occipital_sdChannel_Duration <- Gamma_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Occipital_sdChannel <- merge(Gamma_Occipital_sdChannel_Events, Gamma_Occipital_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Gamma_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  #Beta
  Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
  Beta_Occipital <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "O"))
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Beta_Occipital_avgChannel_Events <- Beta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))
  Beta_Occipital_avgChannel_logPower <- Beta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))
  Beta_Occipital_avgChannel_Power <- Beta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = mean(Beta_Trial_Power, na.rm = T))
  Beta_Occipital_avgChannel_Duration <- Beta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))
  
  Beta_Occipital_avgChannel <- merge(Beta_Occipital_avgChannel_Events, Beta_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Beta_Occipital_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Beta_Occipital_sdChannel_Events <- Beta_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))
  Beta_Occipital_sdChannel_logPower <- Beta_Occipital_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))
  Beta_Occipital_sdChannel_Power <- Beta_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = sd(Beta_Trial_Power, na.rm = T))
  Beta_Occipital_sdChannel_Duration <- Beta_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))
  
  Beta_Occipital_sdChannel <- merge(Beta_Occipital_sdChannel_Events, Beta_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Beta_Occipital_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Alpha
  Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
  Alpha_Occipital <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "O"))
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Alpha_Occipital_avgChannel_Events <- Alpha_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
  Alpha_Occipital_avgChannel_logPower <- Alpha_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))
  Alpha_Occipital_avgChannel_Power <- Alpha_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = mean(Alpha_Trial_Power, na.rm = T))
  Alpha_Occipital_avgChannel_Duration <- Alpha_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Occipital_avgChannel <- merge(Alpha_Occipital_avgChannel_Events, Alpha_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Occipital_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Occipital channels and THEN take SD of the trials 
  Alpha_Occipital_sdChannel_Events <- Alpha_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))
  Alpha_Occipital_sdChannel_logPower <- Alpha_Occipital_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))
  Alpha_Occipital_sdChannel_Power <- Alpha_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = sd(Alpha_Trial_Power, na.rm = T))
  Alpha_Occipital_sdChannel_Duration <- Alpha_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Occipital_sdChannel <- merge(Alpha_Occipital_sdChannel_Events, Alpha_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Occipital_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Theta
  Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
  Theta_Occipital <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "O"))
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Theta_Occipital_avgChannel_Events <- Theta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))
  Theta_Occipital_avgChannel_logPower <- Theta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))
  Theta_Occipital_avgChannel_Power <- Theta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = mean(Theta_Trial_Power, na.rm = T))
  Theta_Occipital_avgChannel_Duration <- Theta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))
  
  Theta_Occipital_avgChannel <- merge(Theta_Occipital_avgChannel_Events, Theta_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Occipital_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Occipital channels and THEN take SD of the trials 
  Theta_Occipital_sdChannel_Events <- Theta_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))
  Theta_Occipital_sdChannel_logPower <- Theta_Occipital_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))
  Theta_Occipital_sdChannel_Power <- Theta_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = sd(Theta_Trial_Power, na.rm = T))
  Theta_Occipital_sdChannel_Duration <- Theta_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))
  
  Theta_Occipital_sdChannel <- merge(Theta_Occipital_sdChannel_Events, Theta_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Theta_Occipital_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  allOccipital <- merge(Gamma_Occipital_avgChannel, Beta_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))
  
  allOccipital_sd <- merge(Gamma_Occipital_sdChannel, Beta_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
  
  alldf$allOccipital <- allOccipital
  alldf$allOccipital_sd <- allOccipital_sd
  
  
  
  #Whole Brain 
  #Gamma
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel

  #keep the data on the trial level, but aggregate all the WholeBrain channels together
  Gamma_WholeBrain_avgChannel_Events <- GammaDelay_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  Gamma_WholeBrain_avgChannel_logPower <- GammaDelay_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  Gamma_WholeBrain_avgChannel_Power <- GammaDelay_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = mean(Gamma_Trial_Power, na.rm = T))
  
  Gamma_WholeBrain_avgChannel_Duration <- GammaDelay_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_WholeBrain_avgChannel <- merge(Gamma_WholeBrain_avgChannel_Events, Gamma_WholeBrain_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_WholeBrain_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Gamma_WholeBrain_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Gamma_WholeBrain_sdChannel_Events <- Gamma_WholeBrain_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  Gamma_WholeBrain_sdChannel_logPower <- Gamma_WholeBrain_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  Gamma_WholeBrain_sdChannel_Power <- Gamma_WholeBrain_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Trial_Power = sd(Gamma_Trial_Power, na.rm = T))
  Gamma_WholeBrain_sdChannel_Duration <- Gamma_WholeBrain_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_WholeBrain_sdChannel <- merge(Gamma_WholeBrain_sdChannel_Events, Gamma_WholeBrain_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_WholeBrain_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Gamma_WholeBrain_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  #Beta
  Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new

  #keep the data on the trial level, but aggregate all the WholeBrain channels together
  Beta_WholeBrain_avgChannel_Events <- Beta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))
  Beta_WholeBrain_avgChannel_logPower <- Beta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))
  Beta_WholeBrain_avgChannel_Power <- Beta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = mean(Beta_Trial_Power, na.rm = T))
  Beta_WholeBrain_avgChannel_Duration <- Beta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))
  
  Beta_WholeBrain_avgChannel <- merge(Beta_WholeBrain_avgChannel_Events, Beta_WholeBrain_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_WholeBrain_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))%>% merge(.,Beta_WholeBrain_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Beta_WholeBrain_sdChannel_Events <- Beta_WholeBrain_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))
  Beta_WholeBrain_sdChannel_logPower <- Beta_WholeBrain_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))
  Beta_WholeBrain_sdChannel_Power <- Beta_WholeBrain_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Trial_Power = sd(Beta_Trial_Power, na.rm = T))
  Beta_WholeBrain_sdChannel_Duration <- Beta_WholeBrain_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))
  
  Beta_WholeBrain_sdChannel <- merge(Beta_WholeBrain_sdChannel_Events, Beta_WholeBrain_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_WholeBrain_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Beta_WholeBrain_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Alpha
  Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new

  #keep the data on the trial level, but aggregate all the WholeBrain channels together
  Alpha_WholeBrain_avgChannel_Events <- Alpha_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
  Alpha_WholeBrain_avgChannel_logPower <- Alpha_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))
  Alpha_WholeBrain_avgChannel_Power <- Alpha_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = mean(Alpha_Trial_Power, na.rm = T))
  Alpha_WholeBrain_avgChannel_Duration <- Alpha_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
  
  Alpha_WholeBrain_avgChannel <- merge(Alpha_WholeBrain_avgChannel_Events, Alpha_WholeBrain_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_WholeBrain_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_WholeBrain_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the WholeBrain channels and THEN take SD of the trials 
  Alpha_WholeBrain_sdChannel_Events <- Alpha_WholeBrain_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))
  Alpha_WholeBrain_sdChannel_logPower <- Alpha_WholeBrain_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))
  Alpha_WholeBrain_sdChannel_Power <- Alpha_WholeBrain_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Trial_Power = sd(Alpha_Trial_Power, na.rm = T))
  Alpha_WholeBrain_sdChannel_Duration <- Alpha_WholeBrain_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))
  
  Alpha_WholeBrain_sdChannel <- merge(Alpha_WholeBrain_sdChannel_Events, Alpha_WholeBrain_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_WholeBrain_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_WholeBrain_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Theta
  Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new

  #keep the data on the trial level, but aggregate all the WholeBrain channels together
  Theta_WholeBrain_avgChannel_Events <- Theta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))
  Theta_WholeBrain_avgChannel_logPower <- Theta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))
  Theta_WholeBrain_avgChannel_Power <- Theta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = mean(Theta_Trial_Power, na.rm = T))
  Theta_WholeBrain_avgChannel_Duration <- Theta_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))
  
  Theta_WholeBrain_avgChannel <- merge(Theta_WholeBrain_avgChannel_Events, Theta_WholeBrain_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_WholeBrain_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_WholeBrain_avgChannel_logPower, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the WholeBrain channels and THEN take SD of the trials 
  Theta_WholeBrain_sdChannel_Events <- Theta_WholeBrain_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))
  Theta_WholeBrain_sdChannel_logPower <- Theta_WholeBrain_avgChannel_logPower %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))
  Theta_WholeBrain_sdChannel_Power <- Theta_WholeBrain_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Trial_Power = sd(Theta_Trial_Power, na.rm = T))
  Theta_WholeBrain_sdChannel_Duration <- Theta_WholeBrain_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))
  
  Theta_WholeBrain_sdChannel <- merge(Theta_WholeBrain_sdChannel_Events, Theta_WholeBrain_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_WholeBrain_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))%>% merge(.,Theta_WholeBrain_sdChannel_logPower, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  allWholeBrain <- merge(Gamma_WholeBrain_avgChannel, Beta_WholeBrain_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_WholeBrain_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_WholeBrain_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))
  
  allWholeBrain_sd <- merge(Gamma_WholeBrain_sdChannel, Beta_WholeBrain_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_WholeBrain_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_WholeBrain_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
  
  alldf$allWholeBrain <- allWholeBrain
  alldf$allWholeBrain_sd <- allWholeBrain_sd
  
  
  
  
  
  
  return(alldf)
}

RegionstoSpectralEventData_SubjectLevel <- function () {
  
  # load in gamma
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20210204.csv')
  FrontalGamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Frontal_trialLevel.csv') 
  ParietalGamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Parietal_trialLevel.csv') 
  OccipitalGamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Occipital_trialLevel.csv') 
  dlpfcGamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_dlpfc_trialLevel.csv') 
  
  colnames(OccipitalGamma)[3:5] <- paste(colnames(OccipitalGamma)[3:5],"Occipital",sep="_")
  colnames(dlpfcGamma)[3:5] <- paste(colnames(dlpfcGamma)[3:5],"dlpfc",sep="_")
  
  allGamma <- merge(FrontalGamma, ParietalGamma, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalGamma, by = c("Subject", "Trial"))%>% merge(., dlpfcGamma, by = c("Subject", "Trial"))
  
  allGamma$log_Gamma_Power_Frontal <- log1p(allGamma$Gamma_Trial_Power_Frontal)
  allGamma$log_Gamma_Power_Parietal <- log1p(allGamma$Gamma_Trial_Power_Parietal)
  allGamma$log_Gamma_Power_Occipital <- log1p(allGamma$Gamma_Trial_Power_Occipital)
  allGamma$log_Gamma_Power_dlpfc <- log1p(allGamma$Gamma_Trial_Power_dlpfc)
  
  
  #load in beta 
  FrontalBeta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Frontal_trialLevel.csv') 
  ParietalBeta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Parietal_trialLevel.csv') 
  OccipitalBeta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Occipital_trialLevel.csv') 
  dlpfcBeta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_dlpfc_trialLevel.csv') 
  
  colnames(OccipitalBeta)[3:5] <- paste(colnames(OccipitalBeta)[3:5],"Occipital",sep="_")
  colnames(dlpfcBeta)[3:5] <- paste(colnames(dlpfcBeta)[3:5],"dlpfc",sep="_")
  
  allBeta <- merge(FrontalBeta, ParietalBeta, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalBeta, by = c("Subject", "Trial"))%>% merge(., dlpfcBeta, by = c("Subject", "Trial"))
  
  allBeta$log_Beta_Power_Frontal <- log1p(allBeta$Beta_Trial_Power_Frontal)
  allBeta$log_Beta_Power_Parietal <- log1p(allBeta$Beta_Trial_Power_Parietal)
  allBeta$log_Beta_Power_Occipital <- log1p(allBeta$Beta_Trial_Power_Occipital)
  allBeta$log_Beta_Power_dlpfc <- log1p(allBeta$Beta_Trial_Power_dlpfc)
  
  #load in Alpha 
  FrontalAlpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_Frontal_data_trialLevel.csv') 
  ParietalAlpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_Parietal_data_trialLevel.csv') 
  OccipitalAlpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_Occipital_data_trialLevel.csv') 
  dlpfcAlpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_dlpfc_data_trialLevel.csv') 
  
  colnames(OccipitalAlpha)[3:5] <- paste(colnames(OccipitalAlpha)[3:5],"Occipital",sep="_")
  colnames(dlpfcAlpha)[3:5] <- paste(colnames(dlpfcAlpha)[3:5],"dlpfc",sep="_")
  
  
  allAlpha <- merge(FrontalAlpha, ParietalAlpha, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalAlpha, by = c("Subject", "Trial"))%>% merge(., dlpfcAlpha, by = c("Subject", "Trial"))
  
  allAlpha$log_Alpha_Power_Frontal <- log1p(allAlpha$Alpha_Trial_Power_Frontal)
  allAlpha$log_Alpha_Power_Parietal <- log1p(allAlpha$Alpha_Trial_Power_Parietal)
  allAlpha$log_Alpha_Power_Occipital <- log1p(allAlpha$Alpha_Trial_Power_Occipital)
  allAlpha$log_Alpha_Power_dlpfc <- log1p(allAlpha$Alpha_Trial_Power_dlpfc)
  
  #load in Theta 
  FrontalTheta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Frontal_trialLevel.csv') 
  ParietalTheta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Parietal_trialLevel.csv') 
  OccipitalTheta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Occipital_trialLevel.csv') 
  dlpfcTheta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_dlpfc_trialLevel.csv') 
  
  
  colnames(OccipitalTheta)[3:5] <- paste(colnames(OccipitalTheta)[3:5],"Occipital",sep="_")
  colnames(dlpfcTheta)[3:5] <- paste(colnames(dlpfcTheta)[3:5],"dlpfc",sep="_")
  
  
  allTheta <- merge(FrontalTheta, ParietalTheta, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalTheta, by = c("Subject", "Trial"))%>% merge(., dlpfcTheta, by = c("Subject", "Trial"))
  
  allTheta$log_Theta_Power_Frontal <- log1p(allTheta$Theta_Trial_Power_Frontal)
  allTheta$log_Theta_Power_Parietal <- log1p(allTheta$Theta_Trial_Power_Parietal)
  allTheta$log_Theta_Power_Occipital <- log1p(allTheta$Theta_Trial_Power_Occipital)
  allTheta$log_Theta_Power_dlpfc <- log1p(allTheta$Theta_Trial_Power_dlpfc)
  
  
  #combine all into one dataset 
  alldata_trialLevel <- merge(allGamma, allBeta, by = c("Subject", "Trial")) %>% merge(., allAlpha, by = c("Subject", "Trial")) %>% merge(., allTheta, by = c("Subject", "Trial")) 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_trialLevel_new <- alldata_trialLevel
  
  cols = names(alldata_trialLevel[c(3:50)])
  for (col in cols) {
    
    alldata_trialLevel_group <- alldata_trialLevel %>% group_by(Subject)
    indx <- outliers(alldata_trialLevel_group[[col]])
    
    alldata_trialLevel_new[[col]] <- Map(replace, alldata_trialLevel_new[[col]], indx, NA)
    alldata_trialLevel_new[[col]] <- as.numeric(alldata_trialLevel_new[[col]])
    
  }  
  

  
  #create subject level data
  avgdata <- aggregate(.~ Subject, data = alldata_trialLevel_new, mean)
  sddata <- aggregate(.~Subject, data = alldata_trialLevel_new, sd)
  
  alldata_sublevel <- merge(avgdata, sddata, by = "Subject", suffixes = c("", "_sd")) %>% merge(., agefile, by = "Subject")
  
  alldata_sublevel_new <- alldata_sublevel
  
  cols = names(alldata_sublevel[c(3:99)])
  for (col in cols) {
    
    indx <- outliers(alldata_sublevel[[col]])
    
    alldata_sublevel_new[[col]] <- Map(replace, alldata_sublevel_new[[col]], indx, NA)
    alldata_sublevel_new[[col]] <- as.numeric(alldata_sublevel_new[[col]])
    
  }  
  
  
  return(alldata_sublevel_new)
}


DelayMinusRest_IndividualChannels_SubLevel <- function () {
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  
  channelLocations <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')
  
  # Gamma
  
  GammaResting <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Resting.csv')  
  
  GammaDelay <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
  
  colnames(GammaResting)[1] <- "idvalues"
  colnames(GammaDelay)[1] <- "idvalues"
  
  GammaResting_Sublevel <- aggregate(.~Channel, GammaResting, mean)
  GammaResting_Sublevel$log_Gamma_Power <- log1p(GammaResting_Sublevel$Gamma_Trial_Power)
  GammaResting_SubLevel_Channel <- merge(GammaResting_Sublevel, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)

  GammaDelay_Sublevel <- aggregate(.~Channel, GammaDelay, mean)
  GammaDelay_Sublevel$log_Gamma_Power <- log1p(GammaDelay_Sublevel$Gamma_Trial_Power)
  GammaDelay_Sublevel_Channel <- merge(GammaDelay_Sublevel, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  GammaDelayRest <- merge(GammaResting_SubLevel_Channel[1:7], GammaDelay_Sublevel_Channel, by = "Channel", suffixes = c("_Rest", "_Delay"))
  
  GammaDelayRest$logGammaPower <- GammaDelayRest$log_Gamma_Power_Delay - GammaDelayRest$log_Gamma_Power_Rest
  GammaDelayRest$gammaEventNumber <- GammaDelayRest$Gamma_Event_Number_Delay - GammaDelayRest$Gamma_Event_Number_Rest
  
}

DelayOnly_peakFreq <- function (){
  
  alldf <- list()
  
  #gamma 
  GammaDelayPeaks3_4 <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_PeakFreq_Power3_4.csv')
  
  GammaDelayPeaks3_4_mean <- setnames(aggregate(. ~ Subject, data = GammaDelayPeaks3_4, mean), c("Subject", "Trial",  "Peak_Frequency", "Peak_Power"))
  
  GammaDelayPeaks3_4_var <- setnames(aggregate(. ~ Subject, data = GammaDelayPeaks3_4, sd), c("Subject", "Trial", "Peak_Frequency", "Peak_Power"))
  
  allGammaDelayPeaks_3_4 <- merge(GammaDelayPeaks3_4_mean, GammaDelayPeaks3_4_var, by = "Subject", suffixes = c("", "_Variability"))
  
  index <- names(allGammaDelayPeaks_3_4[,c(3,4,6,7)])
  names(allGammaDelayPeaks_3_4)[names(allGammaDelayPeaks_3_4) %in% index] <- paste("Gamma",names(allGammaDelayPeaks_3_4[, index]),sep=".")
  gammavars <- grep("Gamma",names(allGammaDelayPeaks_3_4),value=TRUE)
  
  
  write.csv(allGammaDelayPeaks_3_4,'H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/delayOnlyGammaPeaks_3_4_20212403.csv')
  
  allGammaDelayPeaks_3_4 <- merge(allGammaDelayPeaks_3_4, agefile, by = "Subject" )
  allGammaDelayPeaks_3_4 <- allGammaDelayPeaks_3_4[allGammaDelayPeaks_3_4$visitno < 2,]
  
  alldf$delayGammaPeaks <- allGammaDelayPeaks_3_4
  return(alldf)
  
  
}
