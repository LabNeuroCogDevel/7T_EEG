
Behavior_Sublevel_Maria <- function() {
  
  library("dplyr")
  # load behavior
  z_thresh = 2
  data <- merge(
    read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20220819.csv') %>% mutate(Subject = paste0(LunaID, '_', ScanDate)), 
    read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv'), 
    by='Subject')
  #filter out entire trial if it is not to be trusted
  #MP - added best error and mgs latency filtering
  # SM - added VGS latency
  
  Behavior <- data %>%
    mutate(absPositionError = abs(PositionError)) %>%
    mutate(absBestError = abs(BestError)) %>%
    mutate(
      absPositionError = ifelse(absPositionError > 23, NA, absPositionError),
      absBestError = ifelse(absBestError > 23, NA, absBestError))
  
  # filter(absPositionError < 23) %>% 
  #   filter(absBestError < 23) %>%
  #   filter(mgsLatency > .1) %>%
  #   filter(vgsLatency > .1)
  # Compute per-trial z-scores, and remove abs(z)>2 for absPosErr & mgsLat
  # MP - added best error, fixed position error typo
  # SM - added VGS latency
  
  behaviorClean <- Behavior %>%
    mutate(
      mgsLatency = ifelse(mgsLatency < .1, NA, mgsLatency), 
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      vgsLatency = ifelse(vgsLatency < .1, NA, vgsLatency), 
      vgsLatency = ifelse(abs(scale(vgsLatency)[,1])<z_thresh, vgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA)
    )
  # Merge with age file
  # MP - added best error
  # SM - added VGS latency
  
  m_mean <- behaviorClean %>%
    group_by(Subject, age, visitno) %>%
    summarize_at(vars(absPositionError, absBestError, mgsLatency, vgsLatency), mean, na.rm=T) %>%  mutate(inverseAge = 1/age)
  
  m_sd <- behaviorClean %>%
    group_by(Subject, age, visitno) %>%
    summarize_at(vars(absPositionError, absBestError, mgsLatency, vgsLatency), sd, na.rm=T) %>%  mutate(inverseAge = 1/age)
  
  m <- merge(m_mean, m_sd, by = "Subject", suffixes = c("","_sd"))
  
  # Check for subject-level outliers
  # MP - added best error
  # SM - added VGS latency
  m_clean <- m %>% ungroup() %>%
    dplyr::mutate(
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA),
      vgsLatency = ifelse(abs(scale(vgsLatency)[,1])<z_thresh, vgsLatency, NA),
      mgsLatency_sd = ifelse(abs(scale(mgsLatency_sd)[,1])<z_thresh, mgsLatency_sd, NA),
      vgsLatency_sd = ifelse(abs(scale(vgsLatency_sd)[,1])<z_thresh, vgsLatency_sd, NA),
      absPositionError_sd = ifelse(abs(scale(absPositionError_sd)[,1])<z_thresh, absPositionError_sd, NA), 
      absBestError_sd = ifelse(abs(scale(absBestError_sd)[,1])<z_thresh, absBestError_sd, NA)
    )
  summary(lm(absPositionError ~ inverseAge, data=m_clean))
  ggplot(data=m_clean, aes(x=age, y=absPositionError)) + geom_point() + stat_smooth(method='lm', formula='y~I(1/x)')
  
  return(m_clean)
  
  
  
}

Behavior_TrialLevel_Maria <- function() {
  
  library("dplyr")
  # load behavior
  z_thresh = 2
  data <- merge(
    read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200929.csv') %>% mutate(Subject = paste0(LunaID, '_', ScanDate)), 
    read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220914.csv'), 
    by='Subject')
  #filter out entire trial if it is not to be trusted
  #MP - added best error and mgs latency filtering
  Behavior <- data %>%
    filter(visitno == 1) %>%
    mutate(absPositionError = abs(PositionError)) %>%
    mutate(absBestError = abs(BestError)) %>%
    filter(absPositionError < 23) %>% 
    filter(absBestError < 23) %>%
    filter(mgsLatency > .1) %>%
    filter(vgsLatency > .1)
  # Compute per-trial z-scores, and remove abs(z)>2 for absPosErr & mgsLat
  #MP - added best error, fixed position error typo
  behaviorClean <- Behavior %>%
    mutate(
      mgsLatency = ifelse(mgsLatency < .1, NA, mgsLatency), 
      mgsLatency = ifelse(abs(scale(mgsLatency)[,1])<z_thresh, mgsLatency, NA),
      vgsLatency = ifelse(vgsLatency < .1, NA, vgsLatency), 
      vgsLatency = ifelse(abs(scale(vgsLatency)[,1])<z_thresh, vgsLatency, NA),
      absPositionError = ifelse(abs(scale(absPositionError)[,1])<z_thresh, absPositionError, NA), 
      absBestError = ifelse(abs(scale(absBestError)[,1])<z_thresh, absBestError, NA)
    )
  
}

Combine_New_Subs <- function () {
  
  # Gamma Frequency Band 
  newsubsGammaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_newSubs_1_2.csv')
  newsubsGammaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_newSubs_2_3.csv')
  newsubsGammaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_newSubs_3_4.csv')
  newsubsGammaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_newSubs_4_5.csv')
  newsubsGammaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_newSubs_5_6.csv')
  
  GammaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_data_1_2.csv')
  GammaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_all_data_2_3.csv')
  GammaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_all_data_3_4.csv')
  GammaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_all_data_4_5.csv')
  GammaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_all_data_5_6.csv')
  
  GammaDelay1_2 <- rbind(GammaDelay1_2, newsubsGammaDelay1_2)
  GammaDelay2_3 <- rbind(GammaDelay2_3, newsubsGammaDelay2_3)
  GammaDelay3_4 <- rbind(GammaDelay3_4, newsubsGammaDelay3_4)
  GammaDelay4_5 <- rbind(GammaDelay4_5, newsubsGammaDelay4_5)
  GammaDelay5_6 <- rbind(GammaDelay5_6, newsubsGammaDelay5_6)
  
  write.csv(GammaDelay1_2, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_data_1_2_20200810.csv')
  write.csv(GammaDelay2_3, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay2_3/Gamma_all_data_2_3_20200810.csv')
  write.csv(GammaDelay3_4, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_all_data_3_4_20200810.csv')
  write.csv(GammaDelay4_5, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay4_5/Gamma_all_data_4_5_20200810.csv')
  write.csv(GammaDelay5_6, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_all_data_5_6_20200810.csv')
  
  
  #alpha frequency band
  newsubsAlphaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_newSubs_data1_2.csv')
  newsubsAlphaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_newSubs_data2_3.csv')
  newsubsAlphaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_newSubs_data3_4.csv')
  newsubsAlphaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_newSubs_data4_5.csv')
  newsubsAlphaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_newSubs_data5_6.csv')
  
  AlphaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_all_data1_2.csv')
  AlphaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_all_data2_3.csv')
  AlphaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_all_data3_4.csv')
  AlphaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_all_data4_5.csv')
  AlphaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_all_data5_6.csv')
  
  AlphaDelay1_2 <- rbind(AlphaDelay1_2, newsubsAlphaDelay1_2)
  AlphaDelay2_3 <- rbind(AlphaDelay2_3, newsubsAlphaDelay2_3)
  AlphaDelay3_4 <- rbind(AlphaDelay3_4, newsubsAlphaDelay3_4)
  AlphaDelay4_5 <- rbind(AlphaDelay4_5, newsubsAlphaDelay4_5)
  AlphaDelay5_6 <- rbind(AlphaDelay5_6, newsubsAlphaDelay5_6)
  
  write.csv(AlphaDelay1_2, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay1_2/Alpha_all_data_1_2_20200810.csv')
  write.csv(AlphaDelay2_3, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay2_3/Alpha_all_data_2_3_20200810.csv')
  write.csv(AlphaDelay3_4, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_all_data_3_4_20200810.csv')
  write.csv(AlphaDelay4_5, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay4_5/Alpha_all_data_4_5_20200810.csv')
  write.csv(AlphaDelay5_6, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_all_data_5_6_20200810.csv')
  
  
  
  #theta frequency band 
  newsubsThetaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_newSubs_data1_2.csv')
  newsubsThetaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_newSubs_data2_3.csv')
  newsubsThetaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_newSubs_data3_4.csv')
  newsubsThetaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_newSubs_data4_5.csv')
  newsubsThetaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_newSubs_data5_6.csv')
  
  ThetaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_all_data1_2.csv')
  ThetaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_all_data2_3.csv')
  ThetaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_all_data3_4.csv')
  ThetaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_all_data4_5.csv')
  ThetaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_all_data5_6.csv')
  
  ThetaDelay1_2 <- rbind(ThetaDelay1_2, newsubsThetaDelay1_2)
  ThetaDelay2_3 <- rbind(ThetaDelay2_3, newsubsThetaDelay2_3)
  ThetaDelay3_4 <- rbind(ThetaDelay3_4, newsubsThetaDelay3_4)
  ThetaDelay4_5 <- rbind(ThetaDelay4_5, newsubsThetaDelay4_5)
  ThetaDelay5_6 <- rbind(ThetaDelay5_6, newsubsThetaDelay5_6)
  
  
  write.csv(ThetaDelay1_2, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay1_2/Theta_all_data_1_2_20200810.csv')
  write.csv(ThetaDelay2_3, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay2_3/Theta_all_data_2_3_20200810.csv')
  write.csv(ThetaDelay3_4, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_all_data_3_4_20200810.csv')
  write.csv(ThetaDelay4_5, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay4_5/Theta_all_data_4_5_20200810.csv')
  write.csv(ThetaDelay5_6, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_all_data_5_6_20200810.csv')
  
  
  #beta frequency band
  
  newsubsBetaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_newSubs_data1_2.csv')
  newsubsBetaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_newSubs_data2_3.csv')
  newsubsBetaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_newSubs_data3_4.csv')
  newsubsBetaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_newSubs_data4_5.csv')
  newsubsBetaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_newSubs_data5_6.csv')
  
  BetaDelay1_2 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_data1_2.csv')
  BetaDelay2_3 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_all_data2_3.csv')
  BetaDelay3_4 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_all_data3_4.csv')
  BetaDelay4_5 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_all_data4_5.csv')
  BetaDelay5_6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_all_data5_6.csv')
  
  BetaDelay1_2 <- rbind(BetaDelay1_2, newsubsBetaDelay1_2)
  BetaDelay2_3 <- rbind(BetaDelay2_3, newsubsBetaDelay2_3)
  BetaDelay3_4 <- rbind(BetaDelay3_4, newsubsBetaDelay3_4)
  BetaDelay4_5 <- rbind(BetaDelay4_5, newsubsBetaDelay4_5)
  BetaDelay5_6 <- rbind(BetaDelay5_6, newsubsBetaDelay5_6)
  
  
  write.csv(BetaDelay1_2, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_data_1_2_20200810.csv')
  write.csv(BetaDelay2_3, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay2_3/Beta_all_data_2_3_20200810.csv')
  write.csv(BetaDelay3_4, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_all_data_3_4_20200810.csv')
  write.csv(BetaDelay4_5, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay4_5/Beta_all_data_4_5_20200810.csv')
  write.csv(BetaDelay5_6, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_all_data_5_6_20200810.csv')
  
  
  
  # FIXATION PERIOD
  newGammaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_newSubs_FIX.csv')
  newAlphaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_newSubs_datafix.csv')
  newBetaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_newSUbs_datafix.csv')
  newThetaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_newSubs_datafix.csv')
  
  GammaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_all_data_fix.csv')
  AlphaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_all_datafix.csv')
  BetaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_all_datafix.csv')
  ThetaT  <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_all_datafix.csv')
  
  GammaT <- rbind(GammaT, newGammaT)
  AlphaT <- rbind(AlphaT, newAlphaT)
  BetaT <- rbind(BetaT, newBetaT)
  ThetaT <- rbind(ThetaT, newThetaT)
  
  
  write.csv(GammaT, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_all_data_fix_20200810.csv')
  write.csv(AlphaT, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_all_datafix_20200810.csv')
  write.csv(BetaT, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_all_datafix_20200810.csv')
  write.csv(ThetaT, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_all_datafix_20200810.csv')
  
}

Only_Take_One_Delay_Bin_TrialLevel <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_data3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_data_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_data3_4.csv' )
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
  #alldata 
  alldata_TrialLevel <- merge(alpha[,c("Subject","Trial",alphavars)],beta[,c("Subject","Trial",betavars)],by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,theta[,c("Subject","Trial",thetavars)],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,gamma[,c("Subject","Trial",gammavars)],by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  alldata_TrialLevel$Gamma.log1p_Trial_Power <- log1p(alldata_TrialLevel$Gamma.Trial_Power)
  alldata_TrialLevel$Beta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Beta.Trial_Power)
  alldata_TrialLevel$Alpha.log1p_Trial_Power <- log1p(alldata_TrialLevel$Alpha.Trial_Power)
  alldata_TrialLevel$Theta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Theta.Trial_Power)
  
  
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  agefile$inverseAge <- 1/agefile$age
  
  alldata_TrialLevel <- merge(alldata_TrialLevel, agefile, by = c("Subject"))
  alldata_TrialLevel <- alldata_TrialLevel[alldata_TrialLevel$visitno < 2,]
  alldata_TrialLevel <- na.omit(alldata_TrialLevel)
  
  
  
  write.csv(alldata_TrialLevel,"/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/alldata_TrialLevel_20210528.csv")
  
  
  
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

DelayOnly_Sublevel_5_6 <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_Subs_data5_6.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_Subs_data_5_6.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_Subs_data5_6.csv' )
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_Subs_5_6.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
  #alldata
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
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  
  
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

DelayOnly_Sublevel_3_4 <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_data3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_data_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_data3_4.csv' )
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
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
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220914.csv')
  
  
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
  alldf$alldata_delayOnly_subOutliersNotDone <- alldata_delayOnly_age
  
  return(alldf)
  
}

DelayOnly_Triallevel_3_4 <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Subs_data3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Subs_data_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Subs_data3_4.csv' )
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Subs_3_4.csv')
  analysisvars <- grep("Power|Duration|Variability|Delay|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
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
  
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  
  
  alldata_delayOnly_age <- merge(alldata_TrialLevel_new, agefile, by = "Subject")
  
  alldata_delayOnly_age <- subset(alldata_delayOnly_age, visitno == 1)
  
  alldata_delayOnly_age$inverseAge <- 1/alldata_delayOnly_age$age
  
  
  
  return(alldata_delayOnly_age)
  
}


DelayOnly_IndividualChannels_TrialLevel <- function () {
  
  individualChannelDF <- list()
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220914.csv')
  
  channelLocations <-  read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')
  
  # Gamma
  
  GammaDelay <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
  
  
  GammaDelay_Age <- merge(GammaDelay, agefile, by = "Subject")
  GammaDelay_Age$log_Gamma_Power <- log1p(GammaDelay_Age$Gamma_Trial_Power)
  GammaDelay_Age$inverseAge <- 1/GammaDelay_Age$age
  GammaDelay_Age_Channel <- merge(GammaDelay_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  GammaDelay_Age_Channel$Task <- 'Delay'
  

  # outlier detection 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  GammaDelay_Age_Channel_new <- GammaDelay_Age_Channel
  
  cols = names(GammaDelay_Age_Channel[c(4:6,11)])
  for (col in cols) {
    
    GammaDelay_Age_Channel_group <- GammaDelay_Age_Channel %>% group_by(Subject) %>% group_by(Channel)
    indx <- outliers(GammaDelay_Age_Channel_group[[col]])
    
    GammaDelay_Age_Channel_new[[col]] <- Map(replace, GammaDelay_Age_Channel_new[[col]], indx, NA)
    GammaDelay_Age_Channel_new[[col]] <- as.numeric(GammaDelay_Age_Channel_new[[col]])
    
  }  
  
  individualChannelDF$GammaDelay_Age_Channel_new <- GammaDelay_Age_Channel_new
  
  
  # Beta
  Beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  
  
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
  Alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  
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
  Theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  
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

FixOnly_Sublevel <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Subs_dataFix.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Subs_datafix.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Subs_dataFix.csv' )
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Subs_FIX.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
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
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220914.csv')
  
  
  alldata_FixOnly_age <- merge(alldata_subLevel, agefile, by = "Subject")
  
  alldata_FixOnly_age <- subset(alldata_FixOnly_age, visitno == 1)
  
  alldata_FixOnly_age$inverseAge <- 1/alldata_FixOnly_age$age
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_FixOnly_age_new <- alldata_FixOnly_age
  
  cols = names(alldata_FixOnly_age[3:35])
  for (col in cols) {
    
    indx <- outliers(alldata_FixOnly_age[[col]])
    
    alldata_FixOnly_age_new[[col]] <- Map(replace, alldata_FixOnly_age_new[[col]], indx, NA)
    alldata_FixOnly_age_new[[col]] <- as.numeric(alldata_FixOnly_age_new[[col]])
  }
  
  
  
  alldf$alldata_FixOnly <- alldata_FixOnly_age_new
  return(alldf)
  
}

Only_Take_Fix_TrialLevel <- function () {
  
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Subs_dataFix.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  #alpha frequency band 
  alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Subs_datafix.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|X|Group",analysisvars)]
  names(alpha) <- gsub('Alpha_','Alpha.',names(alpha))
  alphavars <- grep("Alpha",names(alpha),value=TRUE)
  
  allAlpha <- alpha[,c("Subject",alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Subs_dataFix.csv' )
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|X|Age|Group",analysisvars)]
  names(beta) <- gsub('Beta_','Beta.',names(beta))
  betavars <- grep("Beta",names(beta),value=TRUE)
  
  allBeta <- beta[,c("Subject",betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Subs_FIX.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(gamma),value=TRUE) ###matches names
  analysisvars <- analysisvars[!grepl("id|age|X",analysisvars)] ### removes "id" "age" names etc
  names(gamma) <- gsub('Gamma_','Gamma.',names(gamma))
  
  
  gammavars <- grep("Gamma",names(gamma),value=TRUE)
  
  allGamma <- gamma[,c("Subject",gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
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
  
  alldata_TrialLevel_new <- merge(alldata_TrialLevel_new, agefile, by = c("Subject"))
  
  alldata_TrialLevel_new <- subset(alldata_TrialLevel_new, visitno == 1)
  
  alldata_TrialLevel_new$inverseAge <- 1/alldata_TrialLevel_new$age
  
  alldf$alldata_Fix_TrialLevel <- alldata_TrialLevel_new
  
  return(alldf)
  
  
}

Resting_State_Data_TrialLevel <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  
  # Setting up the dataframes 
  #theta frequency band 
  Theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_trialLevel_data_RS.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(theta) <- gsub('Theta_','Theta.',names(theta))
  thetavars <- grep("Theta",names(theta),value=TRUE)
  
  allTheta <- theta[,c("Subject", "Trial",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  
  #alpha frequency band 
  Alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_trialLevel_data_RS.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(Alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(Alpha) <- gsub('Alpha_','Alpha.',names(Alpha))
  Alphavars <- grep("Alpha",names(Alpha),value=TRUE)
  
  allAlpha <- Alpha[,c("Subject", "Trial",Alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  Beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_trialLevel_data_RS.csv' )
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(Beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(Beta) <- gsub('Beta_','Beta.',names(Beta))
  Betavars <- grep("Beta",names(Beta),value=TRUE)
  
  allBeta <- Beta[,c("Subject", "Trial",Betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  Gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_trialLevel_data_RS.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(Gamma),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(Gamma) <- gsub('Gamma_','Gamma.',names(Gamma))
  Gammavars <- grep("Gamma",names(Gamma),value=TRUE)
  
  allGamma <- Gamma[,c("Subject", "Trial",Gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
  alldata_TrialLevel <- merge(Alpha,Beta,by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% 
    merge(.,Theta, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,Gamma, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  agefile$Subject <- agefile$idvalues
  agefile$inverseAge <- 1/agefile$age
  
  alldata_TrialLevel <- merge(alldata_TrialLevel, agefile, by = "Subject",  all.x = TRUE, all.y = TRUE)
  
  alldf_RS$alldata_RS_TrialLevel <- alldata_TrialLevel
  
  return(alldf_RS)
  
  
  
  
  
  
}

Resting_State_Data_SubjectLevel <- function () {
  
  alldf <- list()
  
  # Setting up the dataframes 
  #theta frequency band 
  
  # Setting up the dataframes 
  #theta frequency band 
  Theta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_trialLevel_data_RS.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(Theta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(Theta) <- gsub('Theta_','Theta.',names(Theta))
  thetavars <- grep("Theta",names(Theta),value=TRUE)
  
  allTheta <- Theta[,c("Subject", "Trial",thetavars)]
  
  allTheta <- allTheta[complete.cases(allTheta),]
  allTheta <- allTheta[!duplicated(allTheta),]
  
  alldf$allTheta <- allTheta
  
  
  #alpha frequency band 
  Alpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_trialLevel_data_RS.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(Alpha),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(Alpha) <- gsub('Alpha_','Alpha.',names(Alpha))
  Alphavars <- grep("Alpha",names(Alpha),value=TRUE)
  
  allAlpha <- Alpha[,c("Subject", "Trial",Alphavars)]
  
  allAlpha <- allAlpha[complete.cases(allAlpha),]
  allAlpha <- allAlpha[!duplicated(allAlpha),]
  
  alldf$allAlpha <- allAlpha
  
  #beta frequency band 
  Beta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_trialLevel_data_RS.csv' )
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(Beta),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(Beta) <- gsub('Beta_','Beta.',names(Beta))
  Betavars <- grep("Beta",names(Beta),value=TRUE)
  
  allBeta <- Beta[,c("Subject", "Trial",Betavars)]
  
  allBeta <- allBeta[complete.cases(allBeta),]
  allBeta <- allBeta[!duplicated(allBeta),]
  
  alldf$allBeta <- allBeta
  
  #gamma frequency band 
  Gamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_trialLevel_data_RS.csv')
  analysisvars <- grep("Power|Duration|Variability|Fix|Number|visit",names(Gamma),value=TRUE) 
  analysisvars <- analysisvars[!grepl("id|age|Age|Group|X",analysisvars)]
  names(Gamma) <- gsub('Gamma_','Gamma.',names(Gamma))
  Gammavars <- grep("Gamma",names(Gamma),value=TRUE)
  
  allGamma <- Gamma[,c("Subject", "Trial",Gammavars)]
  
  allGamma <- allGamma[complete.cases(allGamma),]
  allGamma <- allGamma[!duplicated(allGamma),]
  
  alldf$allGamma <- allGamma
  
  alldata_TrialLevel <- merge(Alpha,Beta,by= c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% 
    merge(.,Theta, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE) %>% merge(.,Gamma, by=c("Subject","Trial"), all.x = TRUE, all.y = TRUE)
  
  
  
  alldata_TrialLevel$Gamma.log1p_Trial_Power <- log1p(alldata_TrialLevel$Gamma.Trial_Power)
  alldata_TrialLevel$Beta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Beta.Trial_Power)
  alldata_TrialLevel$Alpha.log1p_Trial_Power <- log1p(alldata_TrialLevel$Alpha.Trial_Power)
  alldata_TrialLevel$Theta.log1p_Trial_Power <- log1p(alldata_TrialLevel$Theta.Trial_Power)
  
  # Trying to remove outliers from each variable but keep the trial for all the other variables if the value isnt an outlier 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_TrialLevel_new <- alldata_TrialLevel
  
  cols = names(alldata_TrialLevel[c(3:18)])
  for (col in cols) {
    
    alldata_TrialLevel_grouped <- alldata_TrialLevel %>% group_by(Subject)
    
    indx <- outliers(alldata_TrialLevel_grouped[[col]])
    
    alldata_TrialLevel_new[[col]] <- Map(replace, alldata_TrialLevel_new[[col]], indx, NA)
    alldata_TrialLevel_new[[col]] <- as.numeric(alldata_TrialLevel_new[[col]])
    
  }  
  
  alldata_subLevel_avg <- aggregate(.~ Subject, alldata_TrialLevel_new , mean)
  alldata_subLevel_sd <- aggregate(.~ Subject, alldata_TrialLevel_new , sd)
  alldata_subLevel <- merge(alldata_subLevel_avg, alldata_subLevel_sd, by = ("Subject"), suffixes = c("", "_Variability"))
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220914.csv')
  
  
  alldata_RestOnly_age <- merge(alldata_subLevel, agefile, by = c("Subject"))
  
  alldata_RestOnly_age <- subset(alldata_RestOnly_age, visitno == 1)
  
  alldata_RestOnly_age$inverseAge <- 1/alldata_RestOnly_age$age
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_RestOnly_age_new <- alldata_RestOnly_age
  
  cols = names(alldata_RestOnly_age[3:35])
  for (col in cols) {
    
    indx <- outliers(alldata_RestOnly_age[[col]])
    
    alldata_RestOnly_age_new[[col]] <- Map(replace, alldata_RestOnly_age_new[[col]], indx, NA)
    alldata_RestOnly_age_new[[col]] <- as.numeric(alldata_RestOnly_age_new[[col]])
  }
  
  
  
  alldf$alldata_RestOnly <- alldata_RestOnly_age_new
  return(alldf)
  
  
  
  
  
}

RestDelay_IndividualChannels_TrialLevel <- function () {
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  
  channelLocations <-  read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')
  
  # Gamma
  
  GammaResting <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Resting.csv')  
  
  GammaDelay <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv') 
  
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

RegionstoSpectralEventData_SubjectLevel <- function () {
  
  # load in gamma
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  FrontalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Frontal_trialLevel.csv') 
  ParietalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Parietal_trialLevel.csv') 
  OccipitalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_Occipital_trialLevel.csv') 
  dlpfcGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_dlpfc_trialLevel.csv') 
  selectiveFrontalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay3_4/Gamma_selectiveFrontal_trialLevel.csv') 
  
  
  colnames(OccipitalGamma)[3:5] <- paste(colnames(OccipitalGamma)[3:5],"Occipital",sep="_")
  colnames(dlpfcGamma)[3:5] <- paste(colnames(dlpfcGamma)[3:5],"dlpfc",sep="_")
  colnames(selectiveFrontalGamma)[3:5] <- paste(colnames(selectiveFrontalGamma)[3:5],"selectiveFrontal",sep="_")
  
  
  allGamma <- merge(FrontalGamma, ParietalGamma, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalGamma, by = c("Subject", "Trial"))%>% merge(., dlpfcGamma, by = c("Subject", "Trial"))%>% merge(., selectiveFrontalGamma, by = c("Subject", "Trial"))
  
  allGamma$log1p_Gamma_Power_Frontal <- log1p(allGamma$Gamma_Trial_Power_Frontal)
  allGamma$log1p_Gamma_Power_Parietal <- log1p(allGamma$Gamma_Trial_Power_Parietal)
  allGamma$log1p_Gamma_Power_Occipital <- log1p(allGamma$Gamma_Trial_Power_Occipital)
  allGamma$log1p_Gamma_Power_dlpfc <- log1p(allGamma$Gamma_Trial_Power_dlpfc)
  allGamma$log1p_Gamma_Power_selectiveFrontal <- log1p(allGamma$Gamma_Trial_Power_selectiveFrontal)
  
  
  #load in beta 
  FrontalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Frontal_trialLevel.csv') 
  ParietalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Parietal_trialLevel.csv') 
  OccipitalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_Occipital_trialLevel.csv') 
  dlpfcBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_dlpfc_trialLevel.csv') 
  selectiveFrontalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay3_4/Beta_selectiveFrontal_trialLevel.csv') 
  
  
  colnames(OccipitalBeta)[3:5] <- paste(colnames(OccipitalBeta)[3:5],"Occipital",sep="_")
  colnames(dlpfcBeta)[3:5] <- paste(colnames(dlpfcBeta)[3:5],"dlpfc",sep="_")
  colnames(selectiveFrontalBeta)[3:5] <- paste(colnames(selectiveFrontalBeta)[3:5],"selectiveFrontal",sep="_")
  
  
  allBeta <- merge(FrontalBeta, ParietalBeta, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalBeta, by = c("Subject", "Trial"))%>% merge(., dlpfcBeta, by = c("Subject", "Trial"))%>% merge(., selectiveFrontalBeta, by = c("Subject", "Trial"))
  
  allBeta$log1p_Beta_Power_Frontal <- log1p(allBeta$Beta_Trial_Power_Frontal)
  allBeta$log1p_Beta_Power_Parietal <- log1p(allBeta$Beta_Trial_Power_Parietal)
  allBeta$log1p_Beta_Power_Occipital <- log1p(allBeta$Beta_Trial_Power_Occipital)
  allBeta$log1p_Beta_Power_dlpfc <- log1p(allBeta$Beta_Trial_Power_dlpfc)
  allBeta$log1p_Beta_Power_selectiveFrontal <- log1p(allBeta$Beta_Trial_Power_selectiveFrontal)
  
  
  #load in Alpha 
  FrontalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Frontal_data_trialLevel.csv') 
  ParietalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Parietal_data_trialLevel.csv') 
  OccipitalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_Occipital_data_trialLevel.csv') 
  dlpfcAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_dlpfc_data_trialLevel.csv') 
  selectiveFrontalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay3_4/Alpha_selectiveFrontal_data_trialLevel.csv') 
  
  colnames(OccipitalAlpha)[3:5] <- paste(colnames(OccipitalAlpha)[3:5],"Occipital",sep="_")
  colnames(dlpfcAlpha)[3:5] <- paste(colnames(dlpfcAlpha)[3:5],"dlpfc",sep="_")
  colnames(selectiveFrontalAlpha)[3:5] <- paste(colnames(selectiveFrontalAlpha)[3:5],"selectiveFrontal",sep="_")
  
  
  allAlpha <- merge(FrontalAlpha, ParietalAlpha, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalAlpha, by = c("Subject", "Trial"))%>% merge(., dlpfcAlpha, by = c("Subject", "Trial"))%>% merge(., selectiveFrontalAlpha, by = c("Subject", "Trial"))
  
  allAlpha$log1p_Alpha_Power_Frontal <- log1p(allAlpha$Alpha_Trial_Power_Frontal)
  allAlpha$log1p_Alpha_Power_Parietal <- log1p(allAlpha$Alpha_Trial_Power_Parietal)
  allAlpha$log1p_Alpha_Power_Occipital <- log1p(allAlpha$Alpha_Trial_Power_Occipital)
  allAlpha$log1p_Alpha_Power_dlpfc <- log1p(allAlpha$Alpha_Trial_Power_dlpfc)
  allAlpha$log1p_Alpha_Power_selectiveFrontal <- log1p(allAlpha$Alpha_Trial_Power_selectiveFrontal)
  
  
  #load in Theta 
  FrontalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Frontal_trialLevel.csv') 
  ParietalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Parietal_trialLevel.csv') 
  OccipitalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_Occipital_trialLevel.csv') 
  dlpfcTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_dlpfc_trialLevel.csv') 
  selectiveFrontalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay3_4/Theta_selectiveFrontal_trialLevel.csv') 
  
  
  colnames(OccipitalTheta)[3:5] <- paste(colnames(OccipitalTheta)[3:5],"Occipital",sep="_")
  colnames(dlpfcTheta)[3:5] <- paste(colnames(dlpfcTheta)[3:5],"dlpfc",sep="_")
  colnames(selectiveFrontalTheta)[3:5] <- paste(colnames(selectiveFrontalTheta)[3:5],"selectiveFrontal",sep="_")
  
  
  allTheta <- merge(FrontalTheta, ParietalTheta, by = c("Subject", "Trial"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., OccipitalTheta, by = c("Subject", "Trial"))%>% merge(., dlpfcTheta, by = c("Subject", "Trial"))%>% merge(., selectiveFrontalTheta, by = c("Subject", "Trial"))
  
  allTheta$log1p_Theta_Power_Frontal <- log1p(allTheta$Theta_Trial_Power_Frontal)
  allTheta$log1p_Theta_Power_Parietal <- log1p(allTheta$Theta_Trial_Power_Parietal)
  allTheta$log1p_Theta_Power_Occipital <- log1p(allTheta$Theta_Trial_Power_Occipital)
  allTheta$log1p_Theta_Power_dlpfc <- log1p(allTheta$Theta_Trial_Power_dlpfc)
  allTheta$log1p_Theta_Power_selectiveFrontal <- log1p(allTheta$Theta_Trial_Power_selectiveFrontal)
  
  
  #combine all into one dataset 
  alldata_trialLevel <- merge(allGamma, allBeta, by = c("Subject", "Trial")) %>% merge(., allAlpha, by = c("Subject", "Trial")) %>% merge(., allTheta, by = c("Subject", "Trial")) 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_trialLevel_new <- alldata_trialLevel
  
  cols = names(alldata_trialLevel[c(3:82)])
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
  
  cols = names(alldata_sublevel[c(3:161)])
  for (col in cols) {
    
    indx <- outliers(alldata_sublevel[[col]])
    
    alldata_sublevel_new[[col]] <- Map(replace, alldata_sublevel_new[[col]], indx, NA)
    alldata_sublevel_new[[col]] <- as.numeric(alldata_sublevel_new[[col]])
    
  }  
  alldata_sublevel_new$inverseAge <- 1/alldata_sublevel_new$age
  
  return(alldata_sublevel_new)
}


RegionstoSpectralEventData_SubjectLevel5_6 <- function () {
  
  # load in gamma
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  ParietalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_Parietal_trialLevel.csv') 
  OccipitalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_Occipital_trialLevel.csv') 
  dlpfcGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_dlpfc_trialLevel.csv') 
  selectiveFrontalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay5_6/Gamma_selectiveFrontal_trialLevel.csv') 
  
  
  colnames(OccipitalGamma)[3:5] <- paste(colnames(OccipitalGamma)[3:5],"Occipital",sep="_")
  colnames(dlpfcGamma)[3:5] <- paste(colnames(dlpfcGamma)[3:5],"dlpfc",sep="_")

  
  allGamma <- merge(selectiveFrontalGamma, ParietalGamma, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalGamma, by = c("Subject", "Trial"))%>% merge(., dlpfcGamma, by = c("Subject", "Trial"))
  
  allGamma$log1p_Gamma_Power_Parietal <- log1p(allGamma$Gamma_Trial_Power_Parietal)
  allGamma$log1p_Gamma_Power_Occipital <- log1p(allGamma$Gamma_Trial_Power_Occipital)
  allGamma$log1p_Gamma_Power_dlpfc <- log1p(allGamma$Gamma_Trial_Power_dlpfc)
  allGamma$log1p_Gamma_Power_selectiveFrontal <- log1p(allGamma$Gamma_Trial_Power_selectiveFrontal)
  
  
  #load in beta 
  ParietalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_Parietal_trialLevel.csv') 
  OccipitalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_Occipital_trialLevel.csv') 
  dlpfcBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_dlpfc_trialLevel.csv') 
  selectiveFrontalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay5_6/Beta_selectiveFrontal_trialLevel.csv') 
  
  
  colnames(OccipitalBeta)[3:5] <- paste(colnames(OccipitalBeta)[3:5],"Occipital",sep="_")
  colnames(dlpfcBeta)[3:5] <- paste(colnames(dlpfcBeta)[3:5],"dlpfc",sep="_")

  
  allBeta <- merge(selectiveFrontalBeta, ParietalBeta, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalBeta, by = c("Subject", "Trial"))%>% merge(., dlpfcBeta, by = c("Subject", "Trial"))
  
  allBeta$log1p_Beta_Power_Parietal <- log1p(allBeta$Beta_Trial_Power_Parietal)
  allBeta$log1p_Beta_Power_Occipital <- log1p(allBeta$Beta_Trial_Power_Occipital)
  allBeta$log1p_Beta_Power_dlpfc <- log1p(allBeta$Beta_Trial_Power_dlpfc)
  allBeta$log1p_Beta_Power_selectiveFrontal <- log1p(allBeta$Beta_Trial_Power_selectiveFrontal)
  
  
  #load in Alpha 
  ParietalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_Parietal_data_trialLevel.csv') 
  OccipitalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_Occipital_data_trialLevel.csv') 
  dlpfcAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_dlpfc_data_trialLevel.csv') 
  selectiveFrontalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Delay5_6/Alpha_selectiveFrontal_data_trialLevel.csv') 
  
  colnames(OccipitalAlpha)[3:5] <- paste(colnames(OccipitalAlpha)[3:5],"Occipital",sep="_")
  colnames(dlpfcAlpha)[3:5] <- paste(colnames(dlpfcAlpha)[3:5],"dlpfc",sep="_")

  
  allAlpha <- merge(selectiveFrontalAlpha, ParietalAlpha, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalAlpha, by = c("Subject", "Trial"))%>% merge(., dlpfcAlpha, by = c("Subject", "Trial"))
  
  allAlpha$log1p_Alpha_Power_Parietal <- log1p(allAlpha$Alpha_Trial_Power_Parietal)
  allAlpha$log1p_Alpha_Power_Occipital <- log1p(allAlpha$Alpha_Trial_Power_Occipital)
  allAlpha$log1p_Alpha_Power_dlpfc <- log1p(allAlpha$Alpha_Trial_Power_dlpfc)
  allAlpha$log1p_Alpha_Power_selectiveFrontal <- log1p(allAlpha$Alpha_Trial_Power_selectiveFrontal)
  
  
  #load in Theta 
  ParietalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_Parietal_trialLevel.csv') 
  OccipitalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_Occipital_trialLevel.csv') 
  dlpfcTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_dlpfc_trialLevel.csv') 
  selectiveFrontalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Delay5_6/Theta_selectiveFrontal_trialLevel.csv') 
  
  
  colnames(OccipitalTheta)[3:5] <- paste(colnames(OccipitalTheta)[3:5],"Occipital",sep="_")
  colnames(dlpfcTheta)[3:5] <- paste(colnames(dlpfcTheta)[3:5],"dlpfc",sep="_")

  
  allTheta <- merge(selectiveFrontalTheta, ParietalTheta, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalTheta, by = c("Subject", "Trial"))%>% merge(., dlpfcTheta, by = c("Subject", "Trial"))
  
  allTheta$log1p_Theta_Power_Parietal <- log1p(allTheta$Theta_Trial_Power_Parietal)
  allTheta$log1p_Theta_Power_Occipital <- log1p(allTheta$Theta_Trial_Power_Occipital)
  allTheta$log1p_Theta_Power_dlpfc <- log1p(allTheta$Theta_Trial_Power_dlpfc)
  allTheta$log1p_Theta_Power_selectiveFrontal <- log1p(allTheta$Theta_Trial_Power_selectiveFrontal)
  
  
  #combine all into one dataset 
  alldata_trialLevel <- merge(allGamma, allBeta, by = c("Subject", "Trial")) %>% merge(., allAlpha, by = c("Subject", "Trial")) %>% merge(., allTheta, by = c("Subject", "Trial")) 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_trialLevel_new <- alldata_trialLevel
  
  cols = names(alldata_trialLevel[c(3:66)])
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
  
  cols = names(alldata_sublevel[c(3:131)])
  for (col in cols) {
    
    indx <- outliers(alldata_sublevel[[col]])
    
    alldata_sublevel_new[[col]] <- Map(replace, alldata_sublevel_new[[col]], indx, NA)
    alldata_sublevel_new[[col]] <- as.numeric(alldata_sublevel_new[[col]])
    write.csv(alldata_sublevel_new, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/alldata_regional.csv')
  }  
  
  
  return(alldata_sublevel_new)
}

RegionstoSpectralEventData_SubjectLevelRS <- function () {
  
  # load in gamma
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  ParietalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_Parietal_trialLevel_RS.csv') 
  OccipitalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_Occipital_trialLevel_RS.csv') 
  dlpfcGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_dlpfc_trialLevel_RS.csv') 
  selectiveFrontalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_selectiveFrontal_trialLevel_RS.csv')   
  FrontalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_Frontal_trialLevel_RS.csv') 

  
  
  colnames(OccipitalGamma)[3:5] <- paste(colnames(OccipitalGamma)[3:5],"Occipital",sep="_")
  colnames(dlpfcGamma)[3:5] <- paste(colnames(dlpfcGamma)[3:5],"dlpfc",sep="_")
  colnames(FrontalGamma)[3:5] <- paste(colnames(FrontalGamma)[3:5],"Frontal",sep="_")
  
  
  allGamma <- merge(selectiveFrontalGamma, ParietalGamma, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalGamma, by = c("Subject", "Trial"))%>% merge(., dlpfcGamma, by = c("Subject", "Trial"))%>% merge(., FrontalGamma, by = c("Subject", "Trial"))
  
  allGamma$log1p_Gamma_Power_Parietal <- log1p(allGamma$Gamma_Trial_Power_Parietal)
  allGamma$log1p_Gamma_Power_Occipital <- log1p(allGamma$Gamma_Trial_Power_Occipital)
  allGamma$log1p_Gamma_Power_dlpfc <- log1p(allGamma$Gamma_Trial_Power_dlpfc)
  allGamma$log1p_Gamma_Power_selectiveFrontal <- log1p(allGamma$Gamma_Trial_Power_selectiveFrontal)
  allGamma$log1p_Gamma_Power_Frontal <- log1p(allGamma$Gamma_Trial_Power_Frontal)
  
  
  #load in beta 
  ParietalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_Parietal_trialLevel_RS.csv') 
  OccipitalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_Occipital_trialLevel_RS.csv') 
  dlpfcBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_dlpfc_trialLevel_RS.csv') 
  selectiveFrontalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_selectiveFrontal_trialLevel_RS.csv') 
  FrontalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_Frontal_trialLevel_RS.csv') 
  
  
  colnames(OccipitalBeta)[3:5] <- paste(colnames(OccipitalBeta)[3:5],"Occipital",sep="_")
  colnames(dlpfcBeta)[3:5] <- paste(colnames(dlpfcBeta)[3:5],"dlpfc",sep="_")
  colnames(FrontalBeta)[3:5] <- paste(colnames(FrontalBeta)[3:5],"Frontal",sep="_")
  
  
  allBeta <- merge(selectiveFrontalBeta, ParietalBeta, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalBeta, by = c("Subject", "Trial"))%>% merge(., dlpfcBeta, by = c("Subject", "Trial"))%>% merge(., FrontalBeta, by = c("Subject", "Trial"))
  
  allBeta$log1p_Beta_Power_Parietal <- log1p(allBeta$Beta_Trial_Power_Parietal)
  allBeta$log1p_Beta_Power_Occipital <- log1p(allBeta$Beta_Trial_Power_Occipital)
  allBeta$log1p_Beta_Power_dlpfc <- log1p(allBeta$Beta_Trial_Power_dlpfc)
  allBeta$log1p_Beta_Power_selectiveFrontal <- log1p(allBeta$Beta_Trial_Power_selectiveFrontal)
  allBeta$log1p_Beta_Power_Frontal <- log1p(allBeta$Beta_Trial_Power_Frontal)
  
  
  #load in Alpha 
  ParietalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_Parietal_data_trialLevel.csv') 
  OccipitalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_Occipital_data_trialLevel.csv') 
  dlpfcAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_dlpfc_data_trialLevel.csv') 
  selectiveFrontalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_selectiveFrontal_data_trialLevel.csv') 
  FrontalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_Frontal_data_trialLevel.csv') 
  
  
  colnames(OccipitalAlpha)[3:5] <- paste(colnames(OccipitalAlpha)[3:5],"Occipital",sep="_")
  colnames(dlpfcAlpha)[3:5] <- paste(colnames(dlpfcAlpha)[3:5],"dlpfc",sep="_")
  colnames(FrontalAlpha)[3:5] <- paste(colnames(FrontalAlpha)[3:5],"Frontal",sep="_")
  
  
  allAlpha <- merge(selectiveFrontalAlpha, ParietalAlpha, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalAlpha, by = c("Subject", "Trial"))%>% merge(., dlpfcAlpha, by = c("Subject", "Trial"))%>% merge(., FrontalAlpha, by = c("Subject", "Trial"))
  
  allAlpha$log1p_Alpha_Power_Parietal <- log1p(allAlpha$Alpha_Trial_Power_Parietal)
  allAlpha$log1p_Alpha_Power_Occipital <- log1p(allAlpha$Alpha_Trial_Power_Occipital)
  allAlpha$log1p_Alpha_Power_dlpfc <- log1p(allAlpha$Alpha_Trial_Power_dlpfc)
  allAlpha$log1p_Alpha_Power_selectiveFrontal <- log1p(allAlpha$Alpha_Trial_Power_selectiveFrontal)
  allAlpha$log1p_Alpha_Power_Frontal <- log1p(allAlpha$Alpha_Trial_Power_Frontal)
  
  
  #load in Theta 
  ParietalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_Parietal_trialLevel_RS.csv') 
  OccipitalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_Occipital_trialLevel_RS.csv') 
  dlpfcTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_dlpfc_trialLevel_RS.csv') 
  selectiveFrontalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_selectiveFrontal_trialLevel_RS.csv') 
  FrontalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_Frontal_trialLevel_RS.csv') 
  
  
  colnames(OccipitalTheta)[3:5] <- paste(colnames(OccipitalTheta)[3:5],"Occipital",sep="_")
  colnames(dlpfcTheta)[3:5] <- paste(colnames(dlpfcTheta)[3:5],"dlpfc",sep="_")
  colnames(FrontalTheta)[3:5] <- paste(colnames(FrontalTheta)[3:5],"Frontal",sep="_")
  
  
  
  allTheta <- merge(selectiveFrontalTheta, ParietalTheta, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalTheta, by = c("Subject", "Trial"))%>% merge(., dlpfcTheta, by = c("Subject", "Trial"))%>% merge(., FrontalTheta, by = c("Subject", "Trial"))
  
  allTheta$log1p_Theta_Power_Parietal <- log1p(allTheta$Theta_Trial_Power_Parietal)
  allTheta$log1p_Theta_Power_Occipital <- log1p(allTheta$Theta_Trial_Power_Occipital)
  allTheta$log1p_Theta_Power_dlpfc <- log1p(allTheta$Theta_Trial_Power_dlpfc)
  allTheta$log1p_Theta_Power_selectiveFrontal <- log1p(allTheta$Theta_Trial_Power_selectiveFrontal)
  allTheta$log1p_Theta_Power_Frontal <- log1p(allTheta$Theta_Trial_Power_Frontal)
  
  
  #combine all into one dataset 
  alldata_trialLevel <- merge(allGamma, allBeta, by = c("Subject", "Trial")) %>% merge(., allAlpha, by = c("Subject", "Trial")) %>% merge(., allTheta, by = c("Subject", "Trial")) 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_trialLevel_new <- alldata_trialLevel
  
  cols = names(alldata_trialLevel[c(3:82)])
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
  
  cols = names(alldata_sublevel[c(3:163)])
  for (col in cols) {
    
    indx <- outliers(alldata_sublevel[[col]])
    
    alldata_sublevel_new[[col]] <- Map(replace, alldata_sublevel_new[[col]], indx, NA)
    alldata_sublevel_new[[col]] <- as.numeric(alldata_sublevel_new[[col]])
    
  }  
  
  alldata_sublevel_new$inverseAge <- 1/alldata_sublevel_new$age
  return(alldata_sublevel_new)
}

RegionstoSpectralEventData_SubjectLevel_fix <- function () {
  
  # load in gamma
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  ParietalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Parietal_trialLevel_Fix.csv') 
  OccipitalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Occipital_trialLevel_Fix.csv') 
  dlpfcGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_dlpfc_trialLevel_Fix.csv') 
  selectiveFrontalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_selectiveFrontal_trialLevel_Fix.csv') 
  FrontalGamma <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Frontal_trialLevel_Fix.csv') 
  
  
  colnames(OccipitalGamma)[3:5] <- paste(colnames(OccipitalGamma)[3:5],"Occipital",sep="_")
  colnames(dlpfcGamma)[3:5] <- paste(colnames(dlpfcGamma)[3:5],"dlpfc",sep="_")
  colnames(FrontalGamma)[3:5] <- paste(colnames(FrontalGamma)[3:5],"Frontal",sep="_")
  
  
  allGamma <- merge(selectiveFrontalGamma, ParietalGamma, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalGamma, by = c("Subject", "Trial"))%>% merge(., dlpfcGamma, by = c("Subject", "Trial"))%>% merge(., FrontalGamma, by = c("Subject", "Trial"))
  
  allGamma$log1p_Gamma_Power_Parietal <- log1p(allGamma$Gamma_Trial_Power_Parietal)
  allGamma$log1p_Gamma_Power_Occipital <- log1p(allGamma$Gamma_Trial_Power_Occipital)
  allGamma$log1p_Gamma_Power_dlpfc <- log1p(allGamma$Gamma_Trial_Power_dlpfc)
  allGamma$log1p_Gamma_Power_selectiveFrontal <- log1p(allGamma$Gamma_Trial_Power_selectiveFrontal)
  allGamma$log1p_Gamma_Power_Frontal <- log1p(allGamma$Gamma_Trial_Power_Frontal)
  
  
  #load in beta 
  ParietalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Parietal_trialLevel_Fix.csv') 
  OccipitalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Occipital_trialLevel_Fix.csv') 
  dlpfcBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_dlpfc_trialLevel_Fix.csv') 
  selectiveFrontalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_selectiveFrontal_trialLevel_Fix.csv') 
  FrontalBeta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Frontal_trialLevel_Fix.csv')
  
  
  colnames(OccipitalBeta)[3:5] <- paste(colnames(OccipitalBeta)[3:5],"Occipital",sep="_")
  colnames(dlpfcBeta)[3:5] <- paste(colnames(dlpfcBeta)[3:5],"dlpfc",sep="_")
  colnames(FrontalBeta)[3:5] <- paste(colnames(FrontalBeta)[3:5],"Frontal",sep="_")
  
  
  allBeta <- merge(selectiveFrontalBeta, ParietalBeta, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalBeta, by = c("Subject", "Trial"))%>% merge(., dlpfcBeta, by = c("Subject", "Trial"))%>% merge(., FrontalBeta, by = c("Subject", "Trial"))
  
  allBeta$log1p_Beta_Power_Parietal <- log1p(allBeta$Beta_Trial_Power_Parietal)
  allBeta$log1p_Beta_Power_Occipital <- log1p(allBeta$Beta_Trial_Power_Occipital)
  allBeta$log1p_Beta_Power_dlpfc <- log1p(allBeta$Beta_Trial_Power_dlpfc)
  allBeta$log1p_Beta_Power_selectiveFrontal <- log1p(allBeta$Beta_Trial_Power_selectiveFrontal)
  allBeta$log1p_Beta_Power_Frontal <- log1p(allBeta$Beta_Trial_Power_Frontal)
  
  
  #load in Alpha 
  ParietalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Parietal_trialLevel_Fix.csv') 
  OccipitalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Occipital_trialLevel_Fix.csv') 
  dlpfcAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_dlpfc_trialLevel_Fix.csv') 
  selectiveFrontalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_selectiveFrontal_trialLevel_Fix.csv') 
  FrontalAlpha <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Frontal_trialLevel_Fix.csv') 
  
  colnames(OccipitalAlpha)[3:5] <- paste(colnames(OccipitalAlpha)[3:5],"Occipital",sep="_")
  colnames(dlpfcAlpha)[3:5] <- paste(colnames(dlpfcAlpha)[3:5],"dlpfc",sep="_")
  colnames(FrontalAlpha)[3:5] <- paste(colnames(FrontalAlpha)[3:5],"Frontal",sep="_")
  
  
  allAlpha <- merge(selectiveFrontalAlpha, ParietalAlpha, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalAlpha, by = c("Subject", "Trial"))%>% merge(., dlpfcAlpha, by = c("Subject", "Trial"))%>% merge(., FrontalAlpha, by = c("Subject", "Trial"))
  
  allAlpha$log1p_Alpha_Power_Parietal <- log1p(allAlpha$Alpha_Trial_Power_Parietal)
  allAlpha$log1p_Alpha_Power_Occipital <- log1p(allAlpha$Alpha_Trial_Power_Occipital)
  allAlpha$log1p_Alpha_Power_dlpfc <- log1p(allAlpha$Alpha_Trial_Power_dlpfc)
  allAlpha$log1p_Alpha_Power_selectiveFrontal <- log1p(allAlpha$Alpha_Trial_Power_selectiveFrontal)
  allAlpha$log1p_Alpha_Power_Frontal <- log1p(allAlpha$Alpha_Trial_Power_Frontal)
  
  
  #load in Theta 
  ParietalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Parietal_trialLevel_Fix.csv') 
  OccipitalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Occipital_trialLevel_Fix.csv') 
  dlpfcTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_dlpfc_trialLevel_Fix.csv') 
  selectiveFrontalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_selectiveFrontal_trialLevel_Fix.csv') 
  FrontalTheta <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Frontal_trialLevel_Fix.csv') 
  
  
  colnames(OccipitalTheta)[3:5] <- paste(colnames(OccipitalTheta)[3:5],"Occipital",sep="_")
  colnames(dlpfcTheta)[3:5] <- paste(colnames(dlpfcTheta)[3:5],"dlpfc",sep="_")
  colnames(FrontalTheta)[3:5] <- paste(colnames(FrontalTheta)[3:5],"Frontal",sep="_")
  
  
  
  allTheta <- merge(selectiveFrontalTheta, ParietalTheta, by = c("Subject", "Trial"), suffixes = c("_selectiveFrontal", "_Parietal")) %>% merge(., OccipitalTheta, by = c("Subject", "Trial"))%>% merge(., dlpfcTheta, by = c("Subject", "Trial"))%>% merge(., FrontalTheta, by = c("Subject", "Trial"))
  
  allTheta$log1p_Theta_Power_Parietal <- log1p(allTheta$Theta_Trial_Power_Parietal)
  allTheta$log1p_Theta_Power_Occipital <- log1p(allTheta$Theta_Trial_Power_Occipital)
  allTheta$log1p_Theta_Power_dlpfc <- log1p(allTheta$Theta_Trial_Power_dlpfc)
  allTheta$log1p_Theta_Power_selectiveFrontal <- log1p(allTheta$Theta_Trial_Power_selectiveFrontal)
  allTheta$log1p_Theta_Power_Frontal <- log1p(allTheta$Theta_Trial_Power_Frontal)
  
  
  #combine all into one dataset 
  alldata_trialLevel <- merge(allGamma, allBeta, by = c("Subject", "Trial")) %>% merge(., allAlpha, by = c("Subject", "Trial")) %>% merge(., allTheta, by = c("Subject", "Trial")) 
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  alldata_trialLevel_new <- alldata_trialLevel
  
  cols = names(alldata_trialLevel[c(3:82)])
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
  
  cols = names(alldata_sublevel[c(3:163)])
  for (col in cols) {
    
    indx <- outliers(alldata_sublevel[[col]])
    
    alldata_sublevel_new[[col]] <- Map(replace, alldata_sublevel_new[[col]], indx, NA)
    alldata_sublevel_new[[col]] <- as.numeric(alldata_sublevel_new[[col]])
    
  }  
  
  alldata_sublevel_new$inverseAge <- 1/alldata_sublevel_new$age
  
  return(alldata_sublevel_new)
}




traditionalEEG <- function() {
  
  # Delay Period
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile_20220204.csv')
  agefile$inverseAge <- 1/agefile$age
  
  DelaygammaPower <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/traditionalEEGanalysis/gammaPower_traditionalMethod_Delay.csv')
  DelaygammaPower$logPower <- 20*log10(DelaygammaPower$Power)
  
  DelaygammaPower_age <- merge(agefile, DelaygammaPower, by = "Subject")
  DelaygammaPower_age <- subset(DelaygammaPower_age, visitno == 1)
  
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  DelaygammaPower_age_new <- DelaygammaPower_age
  
  cols = names(DelaygammaPower_age_new[7:8])
  for (col in cols) {
    
    
    indx <- outliers(DelaygammaPower_age_new[[col]])
    
    DelaygammaPower_age_new[[col]] <- Map(replace, DelaygammaPower_age_new[[col]], indx, NA)
    DelaygammaPower_age_new[[col]] <- as.numeric(DelaygammaPower_age_new[[col]])
    
  }  
  
  # Fixation 
  
  FixgammaPower <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/traditionalEEGanalysis/gammaPower_traditionalMethod_Fix.csv')
  FixgammaPower$logPower <- 20*log10(FixgammaPower$Power)
  
  FixgammaPower_age <- merge(agefile, FixgammaPower, by = "Subject")
  FixgammaPower_age <- subset(FixgammaPower_age, visitno == 1)
  
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  FixgammaPower_age_new <- FixgammaPower_age
  
  cols = names(FixgammaPower_age_new[7:8])
  for (col in cols) {
    
    
    indx <- outliers(FixgammaPower_age_new[[col]])
    
    FixgammaPower_age_new[[col]] <- Map(replace, FixgammaPower_age_new[[col]], indx, NA)
    FixgammaPower_age_new[[col]] <- as.numeric(FixgammaPower_age_new[[col]])
    
  } 
  
  
  # Resting State 
  
  RestgammaPower <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/traditionalEEGanalysis/gammaPower_traditionalMethod_Rest.csv')
  RestgammaPower$logPower <- 20*log10(RestgammaPower$Power)
  
  RestgammaPower_age <- merge(agefile, RestgammaPower, by = "Subject")
  RestgammaPower_age <- subset(RestgammaPower_age, visitno == 1)
  
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  RestgammaPower_age_new <- RestgammaPower_age
  
  cols = names(RestgammaPower_age_new[7:8])
  for (col in cols) {
    
    
    indx <- outliers(RestgammaPower_age_new[[col]])
    
    RestgammaPower_age_new[[col]] <- Map(replace, RestgammaPower_age_new[[col]], indx, NA)
    RestgammaPower_age_new[[col]] <- as.numeric(RestgammaPower_age_new[[col]])
    
  } 
  
  

  return(DelaygammaPower_age_new)
  
}

traditionalEEG_IndividualChannels <- function() {
  
  agefile <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/agefile_20220204.csv')
  agefile$inverseAge <- 1/agefile$age
  
  DelaygammaPower <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/traditionalEEGanalysis/gammaPower_traditionalMethod_allChannels_Delay_newMethod.csv')
  
  
  DelaygammaPower_age <- merge(agefile, DelaygammaPower, by = "Subject")
  
  DelaygammaPower_age <- subset(DelaygammaPower_age, visitno  == 1)

  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  DelaygammaPower_age_new <- DelaygammaPower_age
  
  cols = names(DelaygammaPower_age_new[8])
  for (col in cols) {
    
    DelayGammaPower_Age_Channel_group <- DelaygammaPower_age_new %>% group_by(Subject) %>% group_by(Channel)
    
    indx <- outliers(DelayGammaPower_Age_Channel_group[[col]])
    
    DelaygammaPower_age_new[[col]] <- Map(replace, DelaygammaPower_age_new[[col]], indx, NA)
    DelaygammaPower_age_new[[col]] <- as.numeric(DelaygammaPower_age_new[[col]])
    
  }  
  DelaygammaPower_age_new$Epoch <- 'Delay'
  
  
  # Fixation
  

  FixgammaPower <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/traditionalEEGanalysis/gammaPower_traditionalMethod_allChannels_Fix_newMethod.csv')
  
  

  FixgammaPower_age <- merge(agefile, FixgammaPower, by = "Subject")
  
  FixgammaPower_age <- subset(FixgammaPower_age, visitno  == 1)
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  FixgammaPower_age_new <- FixgammaPower_age
  
  cols = names(FixgammaPower_age_new[8])
  for (col in cols) {
    
    FixGammaPower_Age_Channel_group <- FixgammaPower_age_new %>% group_by(Subject) %>% group_by(Channel)
    
    indx <- outliers(FixGammaPower_Age_Channel_group[[col]])
    
    FixgammaPower_age_new[[col]] <- Map(replace, FixgammaPower_age_new[[col]], indx, NA)
    FixgammaPower_age_new[[col]] <- as.numeric(FixgammaPower_age_new[[col]])
    
  }  
  FixgammaPower_age_new$Epoch <- 'Fix'
  
  
  
  # Resting State
 
  RestgammaPower <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/traditionalEEGanalysis/gammaPower_traditionalMethod_allChannels_Rest_newMethod.csv')
  
  

  RestgammaPower_age <- merge(agefile, RestgammaPower, by = "Subject")
  
  RestgammaPower_age <- subset(RestgammaPower_age, visitno  == 1)
  
  outliers <- function(x) {
    
    (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
    
    
  }
  
  RestgammaPower_age_new <- RestgammaPower_age
  
  cols = names(RestgammaPower_age_new[8])
  for (col in cols) {
    
    RestGammaPower_Age_Channel_group <- RestgammaPower_age_new %>% group_by(Subject) %>% group_by(Channel)
    
    indx <- outliers(RestGammaPower_Age_Channel_group[[col]])
    
    RestgammaPower_age_new[[col]] <- Map(replace, RestgammaPower_age_new[[col]], indx, NA)
    RestgammaPower_age_new[[col]] <- as.numeric(RestgammaPower_age_new[[col]])
    
  }  
  RestgammaPower_age_new$Epoch <- 'Rest'
  
  allepochs <- rbind(DelaygammaPower_age_new, FixgammaPower_age_new) %>% rbind(., RestgammaPower_age_new)
  
return(allepochs)
  
}

