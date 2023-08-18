
TopographyAcrossAllParticipants <- function() {
  
  library(eegUtils)
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  
  #Gamma

  Gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.csv')
  Gamma_avgChannels <- aggregate(. ~ Gamma$Channel,Gamma, mean)
  
  channelLocations <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')


  Gamma_avgChannels <- merge(Gamma_avgChannels, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  Gamma_avgChannels$log_Gamma_power <- log1p(Gamma_avgChannels$Gamma_Trial_Power)


  ggplot(Gamma_avgChannels, aes(x = Gamma_avgChannels$X, y = Gamma_avgChannels$Y, fill = Gamma_avgChannels$log_Gamma_power)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Gamma_avgChannels, aes(x = Gamma_avgChannels$X, y = Gamma_avgChannels$Y, fill = Gamma_avgChannels$Gamma_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu")

  ggplot(Gamma_avgChannels, aes(x = Gamma_avgChannels$X, y = Gamma_avgChannels$Y, fill = Gamma_avgChannels$Gamma_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu")



  #Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')
  Beta_avgChannels <- aggregate(. ~ Beta$Channel, Beta, mean)
  
  Beta_avgChannels <-  merge(Beta_avgChannels, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  Beta_avgChannels$log_Beta_Power <- log1p(Beta_avgChannels$Beta_Trial_Power)
  
  ggplot(Beta_avgChannels, aes(x = Beta_avgChannels$X, y = Beta_avgChannels$Y, fill = Beta_avgChannels$log_Beta_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu")

  ggplot(Beta_avgChannels, aes(x = Beta_avgChannels$X, y = Beta_avgChannels$Y, fill = Beta_avgChannels$Beta_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Beta_avgChannels, aes(x = Beta_avgChannels$X, y = Beta_avgChannels$Y, fill = Beta_avgChannels$Beta_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  

  #Theta
  
  Theta <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')
  Theta_avgChannels <- aggregate(. ~ Theta$Channel, Theta, mean)
  
  
  Theta_avgChannels <-  merge(Theta_avgChannels, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  Theta_avgChannels$log_Theta_Power <- log1p(Theta_avgChannels$Theta_Trial_Power)
  
  ggplot(Theta_avgChannels, aes(x = Theta_avgChannels$X, y = Theta_avgChannels$Y, fill = Theta_avgChannels$log_Theta_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Theta_avgChannels, aes(x = Theta_avgChannels$X, y = Theta_avgChannels$Y, fill = Theta_avgChannels$Theta_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu")

  ggplot(Theta_avgChannels, aes(x = Theta_avgChannels$X, y = Theta_avgChannels$Y, fill = Theta_avgChannels$Theta_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  
  
  #Alpha
  
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')
  Alpha_avgChannels <- aggregate(. ~ Alpha$Channel, Alpha, mean)
  
  Alpha_avgChannels <- merge(Alpha_avgChannels, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  Alpha_avgChannels$log_Alpha_Power <- log1p(Alpha_avgChannels$Alpha_Trial_Power)
  
  ggplot(Alpha_avgChannels, aes(x = Alpha_avgChannels$X, y = Alpha_avgChannels$Y, fill = Alpha_avgChannels$log_Alpha_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Alpha_avgChannels, aes(x = Alpha_avgChannels$X, y = Alpha_avgChannels$Y, fill = Alpha_avgChannels$Alpha_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Alpha_avgChannels, aes(x = Alpha_avgChannels$X, y = Alpha_avgChannels$Y, fill = Alpha_avgChannels$Alpha_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  
}

TopographyinAgeGroups <- function () {
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  
  Gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel.csv')  
  colnames(Gamma)[1] <- "idvalues"

  Gamma_avgSubjects <- merge(Gamma, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Gamma_avgSubjects$log_Gamma_Power <- log1p(Gamma_avgSubjects$Gamma_Trial_Power)
  
  #Trial Power
  Gamma_avgSubjects <- aggregate(Gamma_avgSubjects$log_Gamma_Power, by = list(Gamma_avgSubjects$Channel, Gamma_avgSubjects$Group), mean)
  colnames(Gamma_avgSubjects)[1] <- "Channel"
  colnames(Gamma_avgSubjects)[2] <- "Group"
  colnames(Gamma_avgSubjects)[3] <- "log_Gamma_Power"
  
  Gamma_avgSubjects <- merge(Gamma_avgSubjects, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  Gamma_avgSubjects_Group1 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 1,]
  Gamma_avgSubjects_Group2 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 2,]
  Gamma_avgSubjects_Group3 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 3,]
  Gamma_avgSubjects_Group4 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 4,]
  

  TP_group1 <- ggplot(Gamma_avgSubjects_Group1, aes(x = Gamma_avgSubjects_Group1$X, y = Gamma_avgSubjects_Group1$Y, fill = Gamma_avgSubjects_Group1$log_Gamma_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 10-15") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group2 <- ggplot(Gamma_avgSubjects_Group2, aes(x = Gamma_avgSubjects_Group2$X, y = Gamma_avgSubjects_Group2$Y, fill = Gamma_avgSubjects_Group2$log_Gamma_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 16-20") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group3 <- ggplot(Gamma_avgSubjects_Group3, aes(x = Gamma_avgSubjects_Group3$X, y = Gamma_avgSubjects_Group3$Y, fill = Gamma_avgSubjects_Group3$log_Gamma_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 21-25") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group4 <- ggplot(Gamma_avgSubjects_Group4, aes(x = Gamma_avgSubjects_Group4$X, y = Gamma_avgSubjects_Group4$Y, fill = Gamma_avgSubjects_Group4$log_Gamma_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 26+") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  
  cowplot::plot_grid(TP_group1, TP_group2, TP_group3, TP_group4, labels = c("A", "B", "C", "D"), nrow = 2, ncol = 2, label_size = 30)

  
  
  #Number of Events 
  Gamma_avgSubjects <- merge(Gamma, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Gamma_avgSubjects <- aggregate(Gamma_avgSubjects$Gamma_Event_Number, by = list(Gamma_avgSubjects$Channel, Gamma_avgSubjects$Group), mean)
  colnames(Gamma_avgSubjects)[1] <- "Channel"
  colnames(Gamma_avgSubjects)[2] <- "Group"
  colnames(Gamma_avgSubjects)[3] <- "Gamma_Event_Number"
  
  Gamma_avgSubjects <- merge(Gamma_avgSubjects, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  Gamma_avgSubjects_Group1 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 1,]
  Gamma_avgSubjects_Group2 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 2,]
  Gamma_avgSubjects_Group3 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 3,]
  Gamma_avgSubjects_Group4 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 4,]
  
  
  TP_group1 <- ggplot(Gamma_avgSubjects_Group1, aes(x = Gamma_avgSubjects_Group1$X, y = Gamma_avgSubjects_Group1$Y, fill = Gamma_avgSubjects_Group1$Gamma_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 10-15") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group2 <- ggplot(Gamma_avgSubjects_Group2, aes(x = Gamma_avgSubjects_Group2$X, y = Gamma_avgSubjects_Group2$Y, fill = Gamma_avgSubjects_Group2$Gamma_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 16-20") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group3 <- ggplot(Gamma_avgSubjects_Group3, aes(x = Gamma_avgSubjects_Group3$X, y = Gamma_avgSubjects_Group3$Y, fill = Gamma_avgSubjects_Group3$Gamma_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 21-25") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group4 <- ggplot(Gamma_avgSubjects_Group4, aes(x = Gamma_avgSubjects_Group4$X, y = Gamma_avgSubjects_Group4$Y, fill = Gamma_avgSubjects_Group4$Gamma_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 26+") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  
  cowplot::plot_grid(TP_group1, TP_group2, TP_group3, TP_group4, labels = c("A", "B", "C", "D"), nrow = 2, ncol = 2, label_size = 30)
  
  
  #Duration of Events 
  Gamma_avgSubjects <- merge(Gamma, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Gamma_avgSubjects <- aggregate(Gamma_avgSubjects$Gamma_Event_Duration, by = list(Gamma_avgSubjects$Channel, Gamma_avgSubjects$Group), mean)
  colnames(Gamma_avgSubjects)[1] <- "Channel"
  colnames(Gamma_avgSubjects)[2] <- "Group"
  colnames(Gamma_avgSubjects)[3] <- "Gamma_Event_Duration"
  
  Gamma_avgSubjects <- merge(Gamma_avgSubjects, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  Gamma_avgSubjects_Group1 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 1,]
  Gamma_avgSubjects_Group2 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 2,]
  Gamma_avgSubjects_Group3 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 3,]
  Gamma_avgSubjects_Group4 <- Gamma_avgSubjects[Gamma_avgSubjects$Group == 4,]
  
  
  TP_group1 <- ggplot(Gamma_avgSubjects_Group1, aes(x = Gamma_avgSubjects_Group1$X, y = Gamma_avgSubjects_Group1$Y, fill = Gamma_avgSubjects_Group1$Gamma_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 10-15") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group2 <- ggplot(Gamma_avgSubjects_Group2, aes(x = Gamma_avgSubjects_Group2$X, y = Gamma_avgSubjects_Group2$Y, fill = Gamma_avgSubjects_Group2$Gamma_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 16-20") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group3 <- ggplot(Gamma_avgSubjects_Group3, aes(x = Gamma_avgSubjects_Group3$X, y = Gamma_avgSubjects_Group3$Y, fill = Gamma_avgSubjects_Group3$Gamma_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 21-25") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  TP_group4 <- ggplot(Gamma_avgSubjects_Group4, aes(x = Gamma_avgSubjects_Group4$X, y = Gamma_avgSubjects_Group4$Y, fill = Gamma_avgSubjects_Group4$Gamma_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu") + ggtitle("Ages 26+") + theme(plot.title = element_text(face = "bold", size = (30), hjust = 0.5))
  
  
  cowplot::plot_grid(TP_group1, TP_group2, TP_group3, TP_group4, labels = c("A", "B", "C", "D"), nrow = 2, ncol = 2, label_size = 30)
  
  
  
  
}

TopographyAcrossAllParticipants_DelayMinusRest_SubLevel <- function() {
  
  library(eegUtils)
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  
  DelayMinusRest_IndividualChannels_SubLevel() #output of which should be GammaDelayRest
  
  #Gamma
  
  ggplot(GammaDelayRest, aes(x = GammaDelayRest$X, y = GammaDelayRest$Y, fill = GammaDelayRest$logGammaPower)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(GammaDelayRest, aes(x = GammaDelayRest$X, y = GammaDelayRest$Y, fill = GammaDelayRest$gammaEventNumber)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  
  
  
  
  #Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')
  Beta_avgChannels <- aggregate(. ~ Beta$Channel, Beta, mean)
  
  Beta_avgChannels <-  merge(Beta_avgChannels, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  Beta_avgChannels$log_Beta_Power <- log1p(Beta_avgChannels$Beta_Trial_Power)
  
  ggplot(Beta_avgChannels, aes(x = Beta_avgChannels$X, y = Beta_avgChannels$Y, fill = Beta_avgChannels$log_Beta_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Beta_avgChannels, aes(x = Beta_avgChannels$X, y = Beta_avgChannels$Y, fill = Beta_avgChannels$Beta_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Beta_avgChannels, aes(x = Beta_avgChannels$X, y = Beta_avgChannels$Y, fill = Beta_avgChannels$Beta_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  
  
  #Theta
  
  Theta <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')
  Theta_avgChannels <- aggregate(. ~ Theta$Channel, Theta, mean)
  
  
  Theta_avgChannels <-  merge(Theta_avgChannels, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  Theta_avgChannels$log_Theta_Power <- log1p(Theta_avgChannels$Theta_Trial_Power)
  
  ggplot(Theta_avgChannels, aes(x = Theta_avgChannels$X, y = Theta_avgChannels$Y, fill = Theta_avgChannels$log_Theta_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Theta_avgChannels, aes(x = Theta_avgChannels$X, y = Theta_avgChannels$Y, fill = Theta_avgChannels$Theta_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Theta_avgChannels, aes(x = Theta_avgChannels$X, y = Theta_avgChannels$Y, fill = Theta_avgChannels$Theta_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  
  
  #Alpha
  
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')
  Alpha_avgChannels <- aggregate(. ~ Alpha$Channel, Alpha, mean)
  
  Alpha_avgChannels <- merge(Alpha_avgChannels, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  Alpha_avgChannels$log_Alpha_Power <- log1p(Alpha_avgChannels$Alpha_Trial_Power)
  
  ggplot(Alpha_avgChannels, aes(x = Alpha_avgChannels$X, y = Alpha_avgChannels$Y, fill = Alpha_avgChannels$log_Alpha_Power)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Alpha_avgChannels, aes(x = Alpha_avgChannels$X, y = Alpha_avgChannels$Y, fill = Alpha_avgChannels$Alpha_Event_Number)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  ggplot(Alpha_avgChannels, aes(x = Alpha_avgChannels$X, y = Alpha_avgChannels$Y, fill = Alpha_avgChannels$Alpha_Event_Duration)) + geom_topo() + scale_fill_distiller(palette = "RdBu")
  
  
}


DLPFC_analysis <- function () {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  # Gamma

  GammaDelay_Age_Channel <- DelayOnly_IndividualChannels_TrialLevel()

  Gamma_DLPFC <- GammaDelay_Age_Channel %>% filter(Label == "'F3'")
  Gamma_DRPFC <- GammaDelay_Age_Channel %>% filter(Label == "'F4'")
 
  Gamma_DLPFC_all <- merge(Gamma_DLPFC, Gamma_DRPFC, by = c("Subject", "Trial", "age"), suffixes = c("_Left","_Right"))
  
  #create avg gamma number between the left and right DLPFC
  Gamma_DLPFC_all$avg_Gamma_Event_Number <- (Gamma_DLPFC_all$Gamma_Event_Number_Left + Gamma_DLPFC_all$Gamma_Event_Number_Right)/2
  
  #create avg gamma power between the left and right DLPFC
  Gamma_DLPFC_all$avg_log_Gamma_Trial_Power <- (Gamma_DLPFC_all$log_Gamma_Power_Left + Gamma_DLPFC_all$log_Gamma_Power_Right)/2
  
  #clean data using events
  Gamma_DLPFC_all_cleanEvents <- Gamma_DLPFC_all %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2))
  
  #clean data using power
  Gamma_DLPFC_all_cleanPower <- Gamma_DLPFC_all %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2))
  
  
  # PLOTTING EEG MEASURE BY AGE On TRIAL LEVEL
  ## Trial Power
 lunaize(ggplot(Gamma_DLPFC_all_cleanPower, aes(x = age, y = avg_log_Gamma_Trial_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Gamma Power in the DLPFC (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_DLPFC_all_cleanPower, avg_log_Gamma_Trial_Power ~ inverseAge_Left + (1|Subject))
  print(car::Anova(lm.model))
   
  ## Number of events
  
  lunaize(ggplot(Gamma_DLPFC_all_cleanEvents, aes(x = age_Left, y = avg_Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the DLPFC (Delay Only)") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_DLPFC_all_cleanEvents, avg_Gamma_Event_Number ~ inverseAge_Left + (1|Subject))
  print(car::Anova(lm.model))
  
 
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  
  ## Trial Power vs accuracy
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("avg_log_Gamma_Trial_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all, other_vars)
  
  # Mediation on separate age groups 
  Gamma_DLPFC_all_youngest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 1,]
  Gamma_DLPFC_all_secondyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 2,]
  Gamma_DLPFC_all_thirdyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 3,]
  Gamma_DLPFC_all_oldest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 4,]
  
  other_vars <-c("avg_log_Gamma_Trial_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all_oldest, other_vars)
  
  
  
  #Trial power vs  latency
  other_vars <-c("avg_log_Gamma_Trial_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all, other_vars)
  
  
  # Mediation on separate age groups 
  Gamma_DLPFC_all_youngest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 1,]
  Gamma_DLPFC_all_secondyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 2,]
  Gamma_DLPFC_all_thirdyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 3,]
  Gamma_DLPFC_all_oldest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 4,]
  
  other_vars <-c("avg_log_Gamma_Trial_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all_oldest, other_vars)
  
  

  ## Number of Events vs accuracy
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("avg_Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all, other_vars)
  
  # Mediation on separate age groups 
  Gamma_DLPFC_all_youngest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 1,]
  Gamma_DLPFC_all_secondyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 2,]
  Gamma_DLPFC_all_thirdyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 3,]
  Gamma_DLPFC_all_oldest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 4,]
  
  other_vars <-c("avg_Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all_oldest, other_vars)
  
  ## Number of Events vs Latency
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("avg_Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all, other_vars)
  
  # Mediation on separate age groups 
  Gamma_DLPFC_all_youngest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 1,]
  Gamma_DLPFC_all_secondyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 2,]
  Gamma_DLPFC_all_thirdyoungest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 3,]
  Gamma_DLPFC_all_oldest <- Gamma_DLPFC_all[Gamma_DLPFC_all$Group_Left == 4,]
  
  other_vars <-c("avg_Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_DLPFC_all_oldest, other_vars)
  
  # Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  
  colnames(Beta)[1] <- "idvalues"
  
  Beta_Age <- merge(Beta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Beta_Age$log_Beta_Power <- log1p(Beta_Age$Beta_Trial_Power)
  Beta_Age$inverseAge <- 1/Beta_Age$age
  
  Beta_Age_Channel <- merge(Beta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Beta_DLPFC <- Beta_Age_Channel %>% filter(Label == "'F3'")
  Beta_DRPFC <- Beta_Age_Channel %>% filter(Label == "'F4'")
  
  Beta_DLPFC_AvgSubjects <- aggregate(. ~ idvalues, Beta_DLPFC, mean)
  
  # Trial Power
  print(ggplot(Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Beta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the DLPFC") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 ,], log_Beta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_DLPFC_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_DLPFC_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 & Beta_DLPFC_AvgSubjects$Beta_Event_Number < BetaCutoff ,], aes(x = age, y = Beta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the DLPFC") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2  & Beta_DLPFC_AvgSubjects$Beta_Event_Number < BetaCutoff,], Beta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_DLPFC_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_DLPFC_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 & Beta_DLPFC_AvgSubjects$Beta_Event_Duration < BetaCutoff ,], aes(x = age, y = Beta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the DLPFC") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2  & Beta_DLPFC_AvgSubjects$Beta_Event_Duration < BetaCutoff,], Beta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Theta
  
  Theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  colnames(Theta)[1] <- "idvalues"
  
  Theta_Age <- merge(Theta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Theta_Age$log_Theta_Power <- log1p(Theta_Age$Theta_Trial_Power)
  Theta_Age$inverseAge <- 1/Theta_Age$age
  
  Theta_Age_Channel <- merge(Theta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Theta_DLPFC <- Theta_Age_Channel %>% filter(Label == "'F3'")
  Theta_DRPFC <- Theta_Age_Channel %>% filter(Label == "'F4'")
  
  Theta_DLPFC_AvgSubjects <- aggregate(. ~ idvalues, Theta_DLPFC, mean)
  
  # Trial Power
  print(ggplot(Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Theta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the DLPFC") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 ,], log_Theta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_DLPFC_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_DLPFC_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
 print(ggplot(Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 & Theta_DLPFC_AvgSubjects$Theta_Event_Number < ThetaCutoff ,], aes(x = age, y = Theta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the DLPFC") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2  & Theta_DLPFC_AvgSubjects$Theta_Event_Number < ThetaCutoff,], Theta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_DLPFC_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_DLPFC_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 & Theta_DLPFC_AvgSubjects$Theta_Event_Duration < ThetaCutoff ,], aes(x = age, y = Theta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the DLPFC") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2  & Theta_DLPFC_AvgSubjects$Theta_Event_Duration < ThetaCutoff,], Theta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  #Alpha
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  colnames(Alpha)[1] <- "idvalues"
  
  Alpha_Age <- merge(Alpha, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Alpha_Age$log_Alpha_Power <- log1p(Alpha_Age$Alpha_Trial_Power)
  Alpha_Age$inverseAge <- 1/Alpha_Age$age
  
  Alpha_Age_Channel <- merge(Alpha_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Alpha_DLPFC <- Alpha_Age_Channel %>% filter(Label == "'F3'")
  Alpha_DRPFC <- Alpha_Age_Channel %>% filter(Label == "'F4'")
  
  Alpha_DLPFC_AvgSubjects <- aggregate(. ~ idvalues, Alpha_DLPFC, mean)
  
  # Trial Power
  print(ggplot(Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Alpha_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the DLPFC") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 ,], log_Alpha_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_DLPFC_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_DLPFC_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 & Alpha_DLPFC_AvgSubjects$Alpha_Event_Number < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the DLPFC") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2  & Alpha_DLPFC_AvgSubjects$Alpha_Event_Number < AlphaCutoff,], Alpha_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_DLPFC_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_DLPFC_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 & Alpha_DLPFC_AvgSubjects$Alpha_Event_Duration < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the DLPFC") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2  & Alpha_DLPFC_AvgSubjects$Alpha_Event_Duration < AlphaCutoff,], Alpha_Event_Duration ~ inverseAge)
  print(anova(lm.model))
}

FrontalCortex_Analysis <- function () { 
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  # Gamma
  individualChannelDF <- DelayOnly_IndividualChannels_TrialLevel()
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  
  Gamma_Frontal <- filter(GammaDelay_Age_Channel, str_detect(GammaDelay_Age_Channel$Label, "F"))
  

  
  #keep the data on the trial level, but aggregate all the frontal channels together
  Gamma_Frontal_avgChannel_Events <- Gamma_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
    
  Gamma_Frontal_avgChannel_Power <- Gamma_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  
  Gamma_Frontal_avgChannel_Duration <- Gamma_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Frontal_avgChannel <- merge(Gamma_Frontal_avgChannel_Events, Gamma_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Gamma_Frontal_sdChannel_Events <- Gamma_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  
  Gamma_Frontal_sdChannel_Power <- Gamma_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  
  Gamma_Frontal_sdChannel_Duration <- Gamma_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Frontal_sdChannel <- merge(Gamma_Frontal_sdChannel_Events, Gamma_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  

  # PLOTTING EEG MEASURE BY AGE On TRIAL LEVEL
  ## Trial Power
  lunaize(ggplot(Gamma_Frontal_avgChannel, aes(x = age, y = log_Gamma_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Gamma Power in the Frontal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Frontal_avgChannel, log_Gamma_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Gamma_Frontal_sdChannel, aes(x = age, y = log_Gamma_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Gamma Power in the Frontal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Frontal_sdChannel, log_Gamma_Power ~ inverseAge )
  print(car::Anova(lm.model))
  
  ## Number of events
  
  lunaize(ggplot(Gamma_Frontal_avgChannel, aes(x = age, y = Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the Frontal (Delay Only)") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Frontal_avgChannel, Gamma_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Number of events var
  
  lunaize(ggplot(Gamma_Frontal_sdChannel, aes(x = age, y = Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Gamma Bursts in the Frontal (Delay Only)") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Frontal_sdChannel, Gamma_Event_Number ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  ## Duration of events
  
  lunaize(ggplot(Gamma_Frontal_avgChannel, aes(x = age, y = Gamma_Event_Duration)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the Frontal (Delay Only)") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Frontal_avgChannel, Gamma_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Duration of events var
  
  lunaize(ggplot(Gamma_Frontal_sdChannel, aes(x = age, y = Gamma_Event_Duration)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Duration of Gamma Bursts in the Frontal (Delay Only)") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Frontal_sdChannel, Gamma_Event_Duration~ inverseAge )
  print(car::Anova(lm.model))
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  
  ## Trial Power vs accuracy
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("log_Gamma_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Frontal_avgChannel_youngest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 1,]
  Gamma_Frontal_avgChannel_secondyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 2,]
  Gamma_Frontal_avgChannel_thirdyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 3,]
  Gamma_Frontal_avgChannel_oldest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 4,]
  
  other_vars <-c("log_Gamma_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel_oldest, other_vars)
  
  
  #Trial power vs  latency
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Frontal_avgChannel_youngest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 1,]
  Gamma_Frontal_avgChannel_secondyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 2,]
  Gamma_Frontal_avgChannel_thirdyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 3,]
  Gamma_Frontal_avgChannel_oldest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 4,]
  
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel_oldest, other_vars)
  
  
  
  
  ## Number of Events vs accuracy
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Frontal_avgChannel_youngest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 1,]
  Gamma_Frontal_avgChannel_secondyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 2,]
  Gamma_Frontal_avgChannel_thirdyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 3,]
  Gamma_Frontal_avgChannel_oldest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 4,]
  
  other_vars <-c("Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel_thirdyoungest, other_vars)
  
  
  ## Number of Events vs Latency
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Frontal_avgChannel_youngest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 1,]
  Gamma_Frontal_avgChannel_secondyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 2,]
  Gamma_Frontal_avgChannel_thirdyoungest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 3,]
  Gamma_Frontal_avgChannel_oldest <- Gamma_Frontal_avgChannel[Gamma_Frontal_avgChannel$Group == 4,]
  
  other_vars <-c("Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Frontal_avgChannel_oldest, other_vars)
  


# Beta
Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new

Beta_Frontal <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "F"))

#keep the data on the trial level, but aggregate all the frontal channels together
Beta_Frontal_avgChannel_Events <- Beta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))

Beta_Frontal_avgChannel_Power <- Beta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))

Beta_Frontal_avgChannel_Duration <- Beta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))

Beta_Frontal_avgChannel <- merge(Beta_Frontal_avgChannel_Events, Beta_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))


#variability
Beta_Frontal_sdChannel_Events <- Beta_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))

Beta_Frontal_sdChannel_Power <- Beta_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))

Beta_Frontal_sdChannel_Duration <- Beta_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))

Beta_Frontal_sdChannel <- merge(Beta_Frontal_sdChannel_Events, Beta_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))



# Trial Power
print(lunaize(ggplot(Beta_Frontal_avgChannel, aes(x = age, y = log_Beta_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the Frontal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Beta_Frontal_avgChannel, log_Beta_Power ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

## Trial Power variability 
lunaize(ggplot(Beta_Frontal_sdChannel, aes(x = age, y = log_Beta_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Beta Power in the Frontal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))

lm.model <- lm(data = Beta_Frontal_sdChannel, log_Beta_Power ~ inverseAge)
print(car::Anova(lm.model))

# Number of events

print(lunaize(ggplot(Beta_Frontal_avgChannel, aes(x = age, y = Beta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Beta_Frontal_avgChannel, Beta_Event_Number ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

# Number of events Var

print(lunaize(ggplot(Beta_Frontal_sdChannel, aes(x = age, y = Beta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Beta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lm(data = Beta_Frontal_sdChannel, Beta_Event_Number ~ inverseAge)
print(car::Anova(lm.model))

# Duration of Events
print(lunaize(ggplot(Beta_Frontal_avgChannel, aes(x = age, y = Beta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Beta_Frontal_avgChannel, Beta_Event_Duration ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

# Duration of Events var
print(lunaize(ggplot(Beta_Frontal_sdChannel, aes(x = age, y = Beta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Beta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lm(data = Beta_Frontal_sdChannel, Beta_Event_Duration ~ inverseAge )
print(car::Anova(lm.model))


# MEDITATION ANALYSIS 
source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
# the dataframe becomes merged with behavior in the function

other_vars <-c("log_Beta_Power", "mgsLatency") #define what you want in the mediation analysis 

MediationAnalysis(BetaFrontal_final, other_vars)

# Alpha
Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new

Alpha_Frontal <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "F"))

#keep the data on the trial level, but aggregate all the frontal channels together
Alpha_Frontal_avgChannel_Events <- Alpha_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))

Alpha_Frontal_avgChannel_Power <- Alpha_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))

Alpha_Frontal_avgChannel_Duration <- Alpha_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))

Alpha_Frontal_avgChannel <- merge(Alpha_Frontal_avgChannel_Events, Alpha_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))


#variability
#has to be subject level, so average the frontal channels and THEN take SD of the trials 
Alpha_Frontal_sdChannel_Events <- Alpha_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))

Alpha_Frontal_sdChannel_Power <- Alpha_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))

Alpha_Frontal_sdChannel_Duration <- Alpha_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))

Alpha_Frontal_sdChannel <- merge(Alpha_Frontal_sdChannel_Events, Alpha_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))



# Trial Power
print(lunaize(ggplot(Alpha_Frontal_avgChannel, aes(x = age, y = log_Alpha_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the Frontal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Alpha_Frontal_avgChannel, log_Alpha_Power ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

## Trial Power variability 
lunaize(ggplot(Alpha_Frontal_sdChannel, aes(x = age, y = log_Alpha_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Alpha Power in the Frontal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))

lm.model <- lm(data = Alpha_Frontal_sdChannel, log_Alpha_Power ~ inverseAge )
print(car::Anova(lm.model))

# Number of events

print(lunaize(ggplot(Alpha_Frontal_avgChannel, aes(x = age, y = Alpha_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Alpha_Frontal_avgChannel, Alpha_Event_Number ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

# Number of events Var

print(lunaize(ggplot(Alpha_Frontal_sdChannel, aes(x = age, y = Alpha_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Alpha Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lm(data = Alpha_Frontal_sdChannel, Alpha_Event_Number ~ inverseAge )
print(car::Anova(lm.model))

# Duration of Events
print(lunaize(ggplot(Alpha_Frontal_avgChannel, aes(x = age, y = Alpha_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Alpha_Frontal_avgChannel, Alpha_Event_Duration ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

# Duration of Events var
print(lunaize(ggplot(Alpha_Frontal_sdChannel, aes(x = age, y = Alpha_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Alpha Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lm(data = Alpha_Frontal_sdChannel, Alpha_Event_Duration ~ inverseAge )
print(car::Anova(lm.model))


# MEDITATION ANALYSIS 
source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
# the dataframe becomes merged with behavior in the function

other_vars <-c("log_Alpha_Power", "mgsLatency") #define what you want in the mediation analysis 

MediationAnalysis_SubjectLevel(Alpha_Frontal_sdChannel, other_vars)


# Theta

Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new

Theta_Frontal <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "F"))

#keep the data on the trial level, but aggregate all the frontal channels together
Theta_Frontal_avgChannel_Events <- Theta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))

Theta_Frontal_avgChannel_Power <- Theta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))

Theta_Frontal_avgChannel_Duration <- Theta_Frontal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))

Theta_Frontal_avgChannel <- merge(Theta_Frontal_avgChannel_Events, Theta_Frontal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Frontal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))


#variability
#has to be subject level, so average the frontal channels and THEN take SD of the trials 
Theta_Frontal_sdChannel_Events <- Theta_Frontal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))

Theta_Frontal_sdChannel_Power <- Theta_Frontal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))

Theta_Frontal_sdChannel_Duration <- Theta_Frontal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))

Theta_Frontal_sdChannel <- merge(Theta_Frontal_sdChannel_Events, Theta_Frontal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_Frontal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))



# Trial Power
print(lunaize(ggplot(Theta_Frontal_avgChannel, aes(x = age, y = log_Theta_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the Frontal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Theta_Frontal_avgChannel, log_Theta_Power ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

## Trial Power variability 
lunaize(ggplot(Theta_Frontal_sdChannel, aes(x = age, y = log_Theta_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Theta Power in the Frontal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))

lm.model <- lm(data = Theta_Frontal_sdChannel, log_Theta_Power ~ inverseAge )
print(car::Anova(lm.model))

# Number of events

print(lunaize(ggplot(Theta_Frontal_avgChannel, aes(x = age, y = Theta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Theta_Frontal_avgChannel, Theta_Event_Number ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

# Number of events Var

print(lunaize(ggplot(Theta_Frontal_sdChannel, aes(x = age, y = Theta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Theta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lm(data = Theta_Frontal_sdChannel, Theta_Event_Number ~ inverseAge )
print(car::Anova(lm.model))

# Duration of Events
print(lunaize(ggplot(Theta_Frontal_avgChannel, aes(x = age, y = Theta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lmer(data = Theta_Frontal_avgChannel, Theta_Event_Duration ~ inverseAge + (1|Subject))
print(car::Anova(lm.model))

# Duration of Events var
print(lunaize(ggplot(Theta_Frontal_sdChannel, aes(x = age, y = Theta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Theta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))

lm.model <- lm(data = Theta_Frontal_sdChannel, Theta_Event_Duration ~ inverseAge )
print(car::Anova(lm.model))


# MEDITATION ANALYSIS 
source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
# the dataframe becomes merged with behavior in the function

other_vars <-c("log_Theta_Power", "mgsLatency") #define what you want in the mediation analysis 

MediationAnalysis(Theta_Frontal_avgChannel, other_vars)


other_vars <-c("log_Theta_Power", "absPositionError") #define what you want in the mediation analysis 

MediationAnalysis_SubjectLevel(Theta_Frontal_sdChannel, other_vars)
}

ParietalCortex_Analysis <- function () { 
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  # Gamma
  individualChannelDF <- DelayOnly_IndividualChannels_TrialLevel()
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  
  Gamma_Parietal <- filter(GammaDelay_Age_Channel, str_detect(GammaDelay_Age_Channel$Label, "P"))
  
  
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Gamma_Parietal_avgChannel_Events <- Gamma_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  
  Gamma_Parietal_avgChannel_Power <- Gamma_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  
  Gamma_Parietal_avgChannel_Duration <- Gamma_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Parietal_avgChannel <- merge(Gamma_Parietal_avgChannel_Events, Gamma_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Gamma_Parietal_sdChannel_Events <- Gamma_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  
  Gamma_Parietal_sdChannel_Power <- Gamma_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  
  Gamma_Parietal_sdChannel_Duration <- Gamma_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Parietal_sdChannel <- merge(Gamma_Parietal_sdChannel_Events, Gamma_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  # PLOTTING EEG MEASURE BY AGE On TRIAL LEVEL
  ## Trial Power
  lunaize(ggplot(Gamma_Parietal_avgChannel, aes(x = age, y = log_Gamma_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Gamma Power in the Parietal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Parietal_avgChannel, log_Gamma_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Gamma_Parietal_sdChannel, aes(x = age, y = log_Gamma_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Gamma Power in the Parietal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Parietal_sdChannel, log_Gamma_Power ~ inverseAge )
  print(car::Anova(lm.model))
  
  ## Number of events
  
  lunaize(ggplot(Gamma_Parietal_avgChannel, aes(x = age, y = Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the Parietal (Delay Only)") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Parietal_avgChannel, Gamma_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Number of events var
  
  lunaize(ggplot(Gamma_Parietal_sdChannel, aes(x = age, y = Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Gamma Bursts in the Parietal (Delay Only)") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Parietal_sdChannel, Gamma_Event_Number ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  ## Duration of events
  
  lunaize(ggplot(Gamma_Parietal_avgChannel, aes(x = age, y = Gamma_Event_Duration)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the Parietal (Delay Only)") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Parietal_avgChannel, Gamma_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Duration of events var
  
  lunaize(ggplot(Gamma_Parietal_sdChannel, aes(x = age, y = Gamma_Event_Duration)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Duration of Gamma Bursts in the Parietal (Delay Only)") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Parietal_sdChannel, Gamma_Event_Duration~ inverseAge )
  print(car::Anova(lm.model))
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Parietal_avgChannel_youngest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 1,]
  Gamma_Parietal_avgChannel_secondyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 2,]
  Gamma_Parietal_avgChannel_thirdyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 3,]
  Gamma_Parietal_avgChannel_oldest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 4,]
  
  other_vars <-c("log_Gamma_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel_oldest, other_vars)
  
  
  #Trial power vs  latency
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Parietal_avgChannel_youngest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 1,]
  Gamma_Parietal_avgChannel_secondyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 2,]
  Gamma_Parietal_avgChannel_thirdyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 3,]
  Gamma_Parietal_avgChannel_oldest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 4,]
  
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel_oldest, other_vars)
  
  
  
  
  ## Number of Events vs accuracy
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Parietal_avgChannel_youngest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 1,]
  Gamma_Parietal_avgChannel_secondyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 2,]
  Gamma_Parietal_avgChannel_thirdyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 3,]
  Gamma_Parietal_avgChannel_oldest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 4,]
  
  other_vars <-c("Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel_thirdyoungest, other_vars)
  
  
  ## Number of Events vs Latency
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Parietal_avgChannel_youngest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 1,]
  Gamma_Parietal_avgChannel_secondyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 2,]
  Gamma_Parietal_avgChannel_thirdyoungest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 3,]
  Gamma_Parietal_avgChannel_oldest <- Gamma_Parietal_avgChannel[Gamma_Parietal_avgChannel$Group == 4,]
  
  other_vars <-c("Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Parietal_avgChannel_oldest, other_vars)
  
  
  
  # Beta
  Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
  
  Beta_Parietal <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Beta_Parietal_avgChannel_Events <- Beta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))
  
  Beta_Parietal_avgChannel_Power <- Beta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))
  
  Beta_Parietal_avgChannel_Duration <- Beta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))
  
  Beta_Parietal_avgChannel <- merge(Beta_Parietal_avgChannel_Events, Beta_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Beta_Parietal_sdChannel_Events <- Beta_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))
  
  Beta_Parietal_sdChannel_Power <- Beta_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))
  
  Beta_Parietal_sdChannel_Duration <- Beta_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))
  
  Beta_Parietal_sdChannel <- merge(Beta_Parietal_sdChannel_Events, Beta_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Trial Power
  print(lunaize(ggplot(Beta_Parietal_avgChannel, aes(x = age, y = log_Beta_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the Parietal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Beta_Parietal_avgChannel, log_Beta_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Beta_Parietal_sdChannel, aes(x = age, y = log_Beta_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Beta Power in the Parietal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Parietal_sdChannel, log_Beta_Power ~ inverseAge)
  print(car::Anova(lm.model))
  
  # Number of events
  
  print(lunaize(ggplot(Beta_Parietal_avgChannel, aes(x = age, y = Beta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Beta_Parietal_avgChannel, Beta_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Number of events Var
  
  print(lunaize(ggplot(Beta_Parietal_sdChannel, aes(x = age, y = Beta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Beta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Beta_Parietal_sdChannel, Beta_Event_Number ~ inverseAge)
  print(car::Anova(lm.model))
  
  # Duration of Events
  print(lunaize(ggplot(Beta_Parietal_avgChannel, aes(x = age, y = Beta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Beta_Parietal_avgChannel, Beta_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Duration of Events var
  print(lunaize(ggplot(Beta_Parietal_sdChannel, aes(x = age, y = Beta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Beta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Beta_Parietal_sdChannel, Beta_Event_Duration ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("log_Beta_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Beta_Parietal_avgChannel, other_vars)
  
  # Alpha
  Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
  
  Alpha_Parietal <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Alpha_Parietal_avgChannel_Events <- Alpha_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
  
  Alpha_Parietal_avgChannel_Power <- Alpha_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))
  
  Alpha_Parietal_avgChannel_Duration <- Alpha_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Parietal_avgChannel <- merge(Alpha_Parietal_avgChannel_Events, Alpha_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Parietal channels and THEN take SD of the trials 
  Alpha_Parietal_sdChannel_Events <- Alpha_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))
  
  Alpha_Parietal_sdChannel_Power <- Alpha_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))
  
  Alpha_Parietal_sdChannel_Duration <- Alpha_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Parietal_sdChannel <- merge(Alpha_Parietal_sdChannel_Events, Alpha_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Trial Power
  print(lunaize(ggplot(Alpha_Parietal_avgChannel, aes(x = age, y = log_Alpha_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the Parietal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Alpha_Parietal_avgChannel, log_Alpha_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Alpha_Parietal_sdChannel, aes(x = age, y = log_Alpha_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Alpha Power in the Parietal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Parietal_sdChannel, log_Alpha_Power ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Number of events
  
  print(lunaize(ggplot(Alpha_Parietal_avgChannel, aes(x = age, y = Alpha_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the Parietal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Alpha_Parietal_avgChannel, Alpha_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Number of events Var
  
  print(lunaize(ggplot(Alpha_Parietal_sdChannel, aes(x = age, y = Alpha_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Alpha Bursts in the Parietal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Alpha_Parietal_sdChannel, Alpha_Event_Number ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Duration of Events
  print(lunaize(ggplot(Alpha_Parietal_avgChannel, aes(x = age, y = Alpha_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the Parietal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Alpha_Parietal_avgChannel, Alpha_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Duration of Events var
  print(lunaize(ggplot(Alpha_Parietal_sdChannel, aes(x = age, y = Alpha_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Alpha Bursts in the Parietal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Alpha_Parietal_sdChannel, Alpha_Event_Duration ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("log_Alpha_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis_SubjectLevel(Alpha_Parietal_sdChannel, other_vars)
  
  
  # Theta
  
  Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
  
  Theta_Parietal <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the Parietal channels together
  Theta_Parietal_avgChannel_Events <- Theta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))
  
  Theta_Parietal_avgChannel_Power <- Theta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))
  
  Theta_Parietal_avgChannel_Duration <- Theta_Parietal %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))
  
  Theta_Parietal_avgChannel <- merge(Theta_Parietal_avgChannel_Events, Theta_Parietal_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Parietal_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Parietal channels and THEN take SD of the trials 
  Theta_Parietal_sdChannel_Events <- Theta_Parietal_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))
  
  Theta_Parietal_sdChannel_Power <- Theta_Parietal_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))
  
  Theta_Parietal_sdChannel_Duration <- Theta_Parietal_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))
  
  Theta_Parietal_sdChannel <- merge(Theta_Parietal_sdChannel_Events, Theta_Parietal_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_Parietal_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Trial Power
  print(lunaize(ggplot(Theta_Parietal_avgChannel, aes(x = age, y = log_Theta_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the Parietal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Theta_Parietal_avgChannel, log_Theta_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Theta_Parietal_sdChannel, aes(x = age, y = log_Theta_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Theta Power in the Parietal (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Parietal_sdChannel, log_Theta_Power ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Number of events
  
  print(lunaize(ggplot(Theta_Parietal_avgChannel, aes(x = age, y = Theta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Theta_Parietal_avgChannel, Theta_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Number of events Var
  
  print(lunaize(ggplot(Theta_Parietal_sdChannel, aes(x = age, y = Theta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Theta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Theta_Parietal_sdChannel, Theta_Event_Number ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Duration of Events
  print(lunaize(ggplot(Theta_Parietal_avgChannel, aes(x = age, y = Theta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Theta_Parietal_avgChannel, Theta_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Duration of Events var
  print(lunaize(ggplot(Theta_Parietal_sdChannel, aes(x = age, y = Theta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Theta Bursts in the Parietal Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Theta_Parietal_sdChannel, Theta_Event_Duration ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Theta_Event_Duration", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Theta_Parietal_avgChannel, other_vars)
  
  
  other_vars <-c("Theta_Event_Duration", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis_SubjectLevel(Theta_Parietal_sdChannel, other_vars)
  
}

OccipitalCortex_Analysis <- function () { 
  
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  # Gamma
  individualChannelDF <- DelayOnly_IndividualChannels_TrialLevel()
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  
  Gamma_Occipital <- filter(GammaDelay_Age_Channel, str_detect(GammaDelay_Age_Channel$Label, "O"))
  
  
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Gamma_Occipital_avgChannel_Events <- Gamma_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  
  Gamma_Occipital_avgChannel_Power <- Gamma_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  
  Gamma_Occipital_avgChannel_Duration <- Gamma_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Occipital_avgChannel <- merge(Gamma_Occipital_avgChannel_Events, Gamma_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Gamma_Occipital_sdChannel_Events <- Gamma_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  
  Gamma_Occipital_sdChannel_Power <- Gamma_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  
  Gamma_Occipital_sdChannel_Duration <- Gamma_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_Occipital_sdChannel <- merge(Gamma_Occipital_sdChannel_Events, Gamma_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  # PLOTTING EEG MEASURE BY AGE On TRIAL LEVEL
  ## Trial Power
  lunaize(ggplot(Gamma_Occipital_avgChannel, aes(x = age, y = log_Gamma_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Gamma Power in the Occipital (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Occipital_avgChannel, log_Gamma_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Gamma_Occipital_sdChannel, aes(x = age, y = log_Gamma_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Gamma Power in the Occipital (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Occipital_sdChannel, log_Gamma_Power ~ inverseAge )
  print(car::Anova(lm.model))
  
  ## Number of events
  
  lunaize(ggplot(Gamma_Occipital_avgChannel, aes(x = age, y = Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the Occipital (Delay Only)") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Occipital_avgChannel, Gamma_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Number of events var
  
  lunaize(ggplot(Gamma_Occipital_sdChannel, aes(x = age, y = Gamma_Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Gamma Bursts in the Occipital (Delay Only)") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Occipital_sdChannel, Gamma_Event_Number ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  ## Duration of events
  
  lunaize(ggplot(Gamma_Occipital_avgChannel, aes(x = age, y = Gamma_Event_Duration)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the Occipital (Delay Only)") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lmer(data = Gamma_Occipital_avgChannel, Gamma_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Duration of events var
  
  lunaize(ggplot(Gamma_Occipital_sdChannel, aes(x = age, y = Gamma_Event_Duration)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Duration of Gamma Bursts in the Occipital (Delay Only)") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Occipital_sdChannel, Gamma_Event_Duration~ inverseAge )
  print(car::Anova(lm.model))
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Occipital_avgChannel_youngest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 1,]
  Gamma_Occipital_avgChannel_secondyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 2,]
  Gamma_Occipital_avgChannel_thirdyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 3,]
  Gamma_Occipital_avgChannel_oldest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 4,]
  
  other_vars <-c("log_Gamma_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel_oldest, other_vars)
  
  
  #Trial power vs  latency
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Occipital_avgChannel_youngest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 1,]
  Gamma_Occipital_avgChannel_secondyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 2,]
  Gamma_Occipital_avgChannel_thirdyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 3,]
  Gamma_Occipital_avgChannel_oldest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 4,]
  
  other_vars <-c("log_Gamma_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel_oldest, other_vars)
  
  
  
  
  ## Number of Events vs accuracy
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Occipital_avgChannel_youngest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 1,]
  Gamma_Occipital_avgChannel_secondyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 2,]
  Gamma_Occipital_avgChannel_thirdyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 3,]
  Gamma_Occipital_avgChannel_oldest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 4,]
  
  other_vars <-c("Gamma_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel_thirdyoungest, other_vars)
  
  
  ## Number of Events vs Latency
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel, other_vars)
  
  # mediation on age groups 
  Gamma_Occipital_avgChannel_youngest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 1,]
  Gamma_Occipital_avgChannel_secondyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 2,]
  Gamma_Occipital_avgChannel_thirdyoungest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 3,]
  Gamma_Occipital_avgChannel_oldest <- Gamma_Occipital_avgChannel[Gamma_Occipital_avgChannel$Group == 4,]
  
  other_vars <-c("Gamma_Event_Number", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Gamma_Occipital_avgChannel_oldest, other_vars)
  
  
  
  # Beta
  Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
  
  Beta_Occipital <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Beta_Occipital_avgChannel_Events <- Beta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Number = mean(Beta_Event_Number, na.rm = T))
  
  Beta_Occipital_avgChannel_Power <- Beta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Beta_Power = mean(log_Beta_Power, na.rm = T))
  
  Beta_Occipital_avgChannel_Duration <- Beta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = mean(Beta_Event_Duration, na.rm = T))
  
  Beta_Occipital_avgChannel <- merge(Beta_Occipital_avgChannel_Events, Beta_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Beta_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  Beta_Occipital_sdChannel_Events <- Beta_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Number = sd(Beta_Event_Number, na.rm = T))
  
  Beta_Occipital_sdChannel_Power <- Beta_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Beta_Power = sd(log_Beta_Power, na.rm = T))
  
  Beta_Occipital_sdChannel_Duration <- Beta_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Beta_Event_Duration = sd(Beta_Event_Duration, na.rm = T))
  
  Beta_Occipital_sdChannel <- merge(Beta_Occipital_sdChannel_Events, Beta_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Beta_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Trial Power
  print(lunaize(ggplot(Beta_Occipital_avgChannel, aes(x = age, y = log_Beta_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the Occipital Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Beta_Occipital_avgChannel, log_Beta_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Beta_Occipital_sdChannel, aes(x = age, y = log_Beta_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Beta Power in the Occipital (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Occipital_sdChannel, log_Beta_Power ~ inverseAge)
  print(car::Anova(lm.model))
  
  # Number of events
  
  print(lunaize(ggplot(Beta_Occipital_avgChannel, aes(x = age, y = Beta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Beta_Occipital_avgChannel, Beta_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Number of events Var
  
  print(lunaize(ggplot(Beta_Occipital_sdChannel, aes(x = age, y = Beta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Beta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Beta_Occipital_sdChannel, Beta_Event_Number ~ inverseAge)
  print(car::Anova(lm.model))
  
  # Duration of Events
  print(lunaize(ggplot(Beta_Occipital_avgChannel, aes(x = age, y = Beta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Beta_Occipital_avgChannel, Beta_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Duration of Events var
  print(lunaize(ggplot(Beta_Occipital_sdChannel, aes(x = age, y = Beta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Beta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Beta_Occipital_sdChannel, Beta_Event_Duration ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("log_Beta_Power", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis(Beta_Occipital_avgChannel, other_vars)
  
  # Alpha
  Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
  
  Alpha_Occipital <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Alpha_Occipital_avgChannel_Events <- Alpha_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = mean(Alpha_Event_Number, na.rm = T))
  
  Alpha_Occipital_avgChannel_Power <- Alpha_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Alpha_Power = mean(log_Alpha_Power, na.rm = T))
  
  Alpha_Occipital_avgChannel_Duration <- Alpha_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = mean(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Occipital_avgChannel <- merge(Alpha_Occipital_avgChannel_Events, Alpha_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Occipital channels and THEN take SD of the trials 
  Alpha_Occipital_sdChannel_Events <- Alpha_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Number = sd(Alpha_Event_Number, na.rm = T))
  
  Alpha_Occipital_sdChannel_Power <- Alpha_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Alpha_Power = sd(log_Alpha_Power, na.rm = T))
  
  Alpha_Occipital_sdChannel_Duration <- Alpha_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Alpha_Event_Duration = sd(Alpha_Event_Duration, na.rm = T))
  
  Alpha_Occipital_sdChannel <- merge(Alpha_Occipital_sdChannel_Events, Alpha_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Alpha_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Trial Power
  print(lunaize(ggplot(Alpha_Occipital_avgChannel, aes(x = age, y = log_Alpha_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the Occipital Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Alpha_Occipital_avgChannel, log_Alpha_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Alpha_Occipital_sdChannel, aes(x = age, y = log_Alpha_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Alpha Power in the Occipital (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Occipital_sdChannel, log_Alpha_Power ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Number of events
  
  print(lunaize(ggplot(Alpha_Occipital_avgChannel, aes(x = age, y = Alpha_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Alpha_Occipital_avgChannel, Alpha_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Number of events Var
  
  print(lunaize(ggplot(Alpha_Occipital_sdChannel, aes(x = age, y = Alpha_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Alpha Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Alpha_Occipital_sdChannel, Alpha_Event_Number ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Duration of Events
  print(lunaize(ggplot(Alpha_Occipital_avgChannel, aes(x = age, y = Alpha_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Alpha_Occipital_avgChannel, Alpha_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Duration of Events var
  print(lunaize(ggplot(Alpha_Occipital_sdChannel, aes(x = age, y = Alpha_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Alpha Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Alpha_Occipital_sdChannel, Alpha_Event_Duration ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("log_Alpha_Power", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis_SubjectLevel(Alpha_Occipital_sdChannel, other_vars)
  
  
  # Theta
  
  Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
  
  Theta_Occipital <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "F"))
  
  #keep the data on the trial level, but aggregate all the Occipital channels together
  Theta_Occipital_avgChannel_Events <- Theta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Number = mean(Theta_Event_Number, na.rm = T))
  
  Theta_Occipital_avgChannel_Power <- Theta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Theta_Power = mean(log_Theta_Power, na.rm = T))
  
  Theta_Occipital_avgChannel_Duration <- Theta_Occipital %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = mean(Theta_Event_Duration, na.rm = T))
  
  Theta_Occipital_avgChannel <- merge(Theta_Occipital_avgChannel_Events, Theta_Occipital_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Theta_Occipital_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability
  #has to be subject level, so average the Occipital channels and THEN take SD of the trials 
  Theta_Occipital_sdChannel_Events <- Theta_Occipital_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Number = sd(Theta_Event_Number, na.rm = T))
  
  Theta_Occipital_sdChannel_Power <- Theta_Occipital_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Theta_Power = sd(log_Theta_Power, na.rm = T))
  
  Theta_Occipital_sdChannel_Duration <- Theta_Occipital_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Theta_Event_Duration = sd(Theta_Event_Duration, na.rm = T))
  
  Theta_Occipital_sdChannel <- merge(Theta_Occipital_sdChannel_Events, Theta_Occipital_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Theta_Occipital_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  # Trial Power
  print(lunaize(ggplot(Theta_Occipital_avgChannel, aes(x = age, y = log_Theta_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the Occipital Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Theta_Occipital_avgChannel, log_Theta_Power ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  ## Trial Power variability 
  lunaize(ggplot(Theta_Occipital_sdChannel, aes(x = age, y = log_Theta_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Theta Power in the Occipital (Delay Only)") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Occipital_sdChannel, log_Theta_Power ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Number of events
  
  print(lunaize(ggplot(Theta_Occipital_avgChannel, aes(x = age, y = Theta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Theta_Occipital_avgChannel, Theta_Event_Number ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Number of events Var
  
  print(lunaize(ggplot(Theta_Occipital_sdChannel, aes(x = age, y = Theta_Event_Number))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Number of Theta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Theta_Occipital_sdChannel, Theta_Event_Number ~ inverseAge )
  print(car::Anova(lm.model))
  
  # Duration of Events
  print(lunaize(ggplot(Theta_Occipital_avgChannel, aes(x = age, y = Theta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lmer(data = Theta_Occipital_avgChannel, Theta_Event_Duration ~ inverseAge + (1|Subject))
  print(car::Anova(lm.model))
  
  # Duration of Events var
  print(lunaize(ggplot(Theta_Occipital_sdChannel, aes(x = age, y = Theta_Event_Duration))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Var Duration of Theta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5))))
  
  lm.model <- lm(data = Theta_Occipital_sdChannel, Theta_Event_Duration ~ inverseAge )
  print(car::Anova(lm.model))
  
  
  # MEDITATION ANALYSIS 
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")
  # the dataframe becomes merged with behavior in the function
  
  other_vars <-c("Theta_Event_Duration", "mgsLatency") #define what you want in the mediation analysis 
  
  MediationAnalysis(Theta_Occipital_avgChannel, other_vars)
  
  
  other_vars <-c("Theta_Event_Number", "absPositionError") #define what you want in the mediation analysis 
  
  MediationAnalysis_SubjectLevel(Theta_Occipital_sdChannel, other_vars)
  
}

CentralCortex_Analysis <- function () { 
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  channelLocations <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv')
  
  # Gamma
  
  Gamma <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel.csv')  
  colnames(Gamma)[1] <- "idvalues"
  
  Gamma_Age <- merge(Gamma, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Gamma_Age$log_Gamma_Power <- log1p(Gamma_Age$Gamma_Trial_Power)
  Gamma_Age$inverseAge <- 1/Gamma_Age$age
  
  Gamma_Age_Channel <- merge(Gamma_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Gamma_Central <- filter(Gamma_Age_Channel, str_detect(Gamma_Age_Channel$Label, "C"))
  
  Gamma_Central_AvgSubjects <- aggregate(. ~ idvalues, Gamma_Central, mean)
  
  # Trial Power
  print(ggplot(Gamma_Central_AvgSubjects[Gamma_Central_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Gamma_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Gamma Power in the Central Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Central_AvgSubjects[Gamma_Central_AvgSubjects$visitno < 2 ,], log_Gamma_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Central_AvgSubjects, mean)
  sdGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Central_AvgSubjects, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  print(ggplot(Gamma_Central_AvgSubjects[Gamma_Central_AvgSubjects$visitno < 2 & Gamma_Central_AvgSubjects$Gamma_Event_Number < gammaCutoff ,], aes(x = age, y = Gamma_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the Central Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Central_AvgSubjects[Gamma_Central_AvgSubjects$visitno < 2  & Gamma_Central_AvgSubjects$Gamma_Event_Number < gammaCutoff,], Gamma_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Central_AvgSubjects, mean)
  sdGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Central_AvgSubjects, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  print(ggplot(Gamma_Central_AvgSubjects[Gamma_Central_AvgSubjects$visitno < 2 & Gamma_Central_AvgSubjects$Gamma_Event_Duration < gammaCutoff ,], aes(x = age, y = Gamma_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the Central Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Central_AvgSubjects[Gamma_Central_AvgSubjects$visitno < 2  & Gamma_Central_AvgSubjects$Gamma_Event_Duration < gammaCutoff,], Gamma_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  # Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  
  colnames(Beta)[1] <- "idvalues"
  
  Beta_Age <- merge(Beta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Beta_Age$log_Beta_Power <- log1p(Beta_Age$Beta_Trial_Power)
  Beta_Age$inverseAge <- 1/Beta_Age$age
  
  Beta_Age_Channel <- merge(Beta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Beta_Central <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "C"))
  
  Beta_Central_AvgSubjects <- aggregate(. ~ idvalues, Beta_Central, mean)
  
  # Trial Power
  print(ggplot(Beta_Central_AvgSubjects[Beta_Central_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Beta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the Central Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Central_AvgSubjects[Beta_Central_AvgSubjects$visitno < 2 ,], log_Beta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Central_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Central_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Central_AvgSubjects[Beta_Central_AvgSubjects$visitno < 2 & Beta_Central_AvgSubjects$Beta_Event_Number < BetaCutoff ,], aes(x = age, y = Beta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the Central Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Central_AvgSubjects[Beta_Central_AvgSubjects$visitno < 2  & Beta_Central_AvgSubjects$Beta_Event_Number < BetaCutoff,], Beta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Central_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Central_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Central_AvgSubjects[Beta_Central_AvgSubjects$visitno < 2 & Beta_Central_AvgSubjects$Beta_Event_Duration < BetaCutoff ,], aes(x = age, y = Beta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the Central Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Central_AvgSubjects[Beta_Central_AvgSubjects$visitno < 2  & Beta_Central_AvgSubjects$Beta_Event_Duration < BetaCutoff,], Beta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Alpha
  
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  colnames(Alpha)[1] <- "idvalues"
  
  Alpha_Age <- merge(Alpha, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Alpha_Age$log_Alpha_Power <- log1p(Alpha_Age$Alpha_Trial_Power)
  Alpha_Age$inverseAge <- 1/Alpha_Age$age
  
  Alpha_Age_Channel <- merge(Alpha_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Alpha_Central <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "C"))
  
  Alpha_Central_AvgSubjects <- aggregate(. ~ idvalues, Alpha_Central, mean)
  
  # Trial Power
  print(ggplot(Alpha_Central_AvgSubjects[Alpha_Central_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Alpha_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the Central Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Central_AvgSubjects[Alpha_Central_AvgSubjects$visitno < 2 ,], log_Alpha_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Central_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Central_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Central_AvgSubjects[Alpha_Central_AvgSubjects$visitno < 2 & Alpha_Central_AvgSubjects$Alpha_Event_Number < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the Central Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Central_AvgSubjects[Alpha_Central_AvgSubjects$visitno < 2  & Alpha_Central_AvgSubjects$Alpha_Event_Number < AlphaCutoff,], Alpha_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Central_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Central_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Central_AvgSubjects[Alpha_Central_AvgSubjects$visitno < 2 & Alpha_Central_AvgSubjects$Alpha_Event_Duration < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the Central Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Central_AvgSubjects[Alpha_Central_AvgSubjects$visitno < 2  & Alpha_Central_AvgSubjects$Alpha_Event_Duration < AlphaCutoff,], Alpha_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Theta
  
  Theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  colnames(Theta)[1] <- "idvalues"
  
  Theta_Age <- merge(Theta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Theta_Age$log_Theta_Power <- log1p(Theta_Age$Theta_Trial_Power)
  Theta_Age$inverseAge <- 1/Theta_Age$age
  
  Theta_Age_Channel <- merge(Theta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Theta_Central <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "C"))
  
  Theta_Central_AvgSubjects <- aggregate(. ~ idvalues, Theta_Central, mean)
  
  # Trial Power
  print(ggplot(Theta_Central_AvgSubjects[Theta_Central_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Theta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the Central Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Central_AvgSubjects[Theta_Central_AvgSubjects$visitno < 2 ,], log_Theta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Central_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Central_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Central_AvgSubjects[Theta_Central_AvgSubjects$visitno < 2 & Theta_Central_AvgSubjects$Theta_Event_Number < ThetaCutoff ,], aes(x = age, y = Theta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the Central Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Central_AvgSubjects[Theta_Central_AvgSubjects$visitno < 2  & Theta_Central_AvgSubjects$Theta_Event_Number < ThetaCutoff,], Theta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Central_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Central_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Central_AvgSubjects[Theta_Central_AvgSubjects$visitno < 2 & Theta_Central_AvgSubjects$Theta_Event_Duration < ThetaCutoff ,], aes(x = age, y = Theta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the Central Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Central_AvgSubjects[Theta_Central_AvgSubjects$visitno < 2  & Theta_Central_AvgSubjects$Theta_Event_Duration < ThetaCutoff,], Theta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
}

WholeCortex_Analysis <- function () {
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  individualChannelDF <- DelayOnly_IndividualChannels_TrialLevel()
  GammaDelay_Age_Channel <- individualChannelDF$GammaDelay_Age_Channel
  
  Gamma_WholeBrain_avgChannel_Events <- GammaDelay_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = mean(Gamma_Event_Number, na.rm = T))
  
  Gamma_WholeBrain_avgChannel_Power <- GammaDelay_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(log_Gamma_Power = mean(log_Gamma_Power, na.rm = T))
  
  Gamma_WholeBrain_avgChannel_Duration <- GammaDelay_Age_Channel %>% group_by(Subject, Trial, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = mean(Gamma_Event_Duration, na.rm = T))
  
  Gamma_WholeBrain_avgChannel <- merge(Gamma_WholeBrain_avgChannel_Events, Gamma_WholeBrain_avgChannel_Power, by = c("Subject", "Trial", "age", "inverseAge", "Group")) %>% merge(.,Gamma_WholeBrain_avgChannel_Duration, by = c("Subject", "Trial", "age", "inverseAge", "Group"))
  
  
  #variability (subject level)
  Gamma_WholeBrain_sdChannel_Events <- Gamma_WholeBrain_avgChannel_Events %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Number = sd(Gamma_Event_Number, na.rm = T))
  
  Gamma_WholeBrain_sdChannel_Power <- Gamma_WholeBrain_avgChannel_Power %>% group_by(Subject, age, inverseAge, Group) %>% summarise(log_Gamma_Power = sd(log_Gamma_Power, na.rm = T))
  
  Gamma_WholeBrain_sdChannel_Duration <- Gamma_WholeBrain_avgChannel_Duration %>% group_by(Subject, age, inverseAge, Group) %>% summarise(Gamma_Event_Duration = sd(Gamma_Event_Duration, na.rm = T))
  
  Gamma_WholeBrain_sdChannel <- merge(Gamma_WholeBrain_sdChannel_Events, Gamma_WholeBrain_sdChannel_Power, by = c("Subject", "age", "inverseAge", "Group")) %>% merge(.,Gamma_WholeBrain_sdChannel_Duration, by = c("Subject", "age", "inverseAge", "Group"))
  
  
  
  
  # Trial Power
  print(ggplot(Gamma_WholeBrain_avgChannel[], aes(x = age, y = log_Gamma_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Gamma Power in the Whole Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Whole_AvgSubjects[Gamma_Whole_AvgSubjects$visitno < 2 ,], log_Gamma_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Whole_AvgSubjects, mean)
  sdGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Whole_AvgSubjects, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  print(ggplot(Gamma_Whole_AvgSubjects[Gamma_Whole_AvgSubjects$visitno < 2 & Gamma_Whole_AvgSubjects$Gamma_Event_Number < gammaCutoff ,], aes(x = age, y = Gamma_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the Whole Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Whole_AvgSubjects[Gamma_Whole_AvgSubjects$visitno < 2  & Gamma_Whole_AvgSubjects$Gamma_Event_Number < gammaCutoff,], Gamma_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Whole_AvgSubjects, mean)
  sdGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Whole_AvgSubjects, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  print(ggplot(Gamma_Whole_AvgSubjects[Gamma_Whole_AvgSubjects$visitno < 2 & Gamma_Whole_AvgSubjects$Gamma_Event_Duration < gammaCutoff ,], aes(x = age, y = Gamma_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the Whole Brain") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Whole_AvgSubjects[Gamma_Whole_AvgSubjects$visitno < 2  & Gamma_Whole_AvgSubjects$Gamma_Event_Duration < gammaCutoff,], Gamma_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  # Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  
  colnames(Beta)[1] <- "idvalues"
  
  Beta_Age <- merge(Beta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Beta_Age$log_Beta_Power <- log1p(Beta_Age$Beta_Trial_Power)
  Beta_Age$inverseAge <- 1/Beta_Age$age
  
  Beta_Age_Channel <- merge(Beta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Beta_Whole_AvgSubjects <- aggregate(. ~ idvalues, Beta_Age_Channel, mean)
  
  # Trial Power
  print(ggplot(Beta_Whole_AvgSubjects[Beta_Whole_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Beta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the Whole Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Whole_AvgSubjects[Beta_Whole_AvgSubjects$visitno < 2 ,], log_Beta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Whole_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Whole_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Whole_AvgSubjects[Beta_Whole_AvgSubjects$visitno < 2 & Beta_Whole_AvgSubjects$Beta_Event_Number < BetaCutoff ,], aes(x = age, y = Beta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the Whole Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Whole_AvgSubjects[Beta_Whole_AvgSubjects$visitno < 2  & Beta_Whole_AvgSubjects$Beta_Event_Number < BetaCutoff,], Beta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Whole_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Whole_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Whole_AvgSubjects[Beta_Whole_AvgSubjects$visitno < 2 & Beta_Whole_AvgSubjects$Beta_Event_Duration < BetaCutoff ,], aes(x = age, y = Beta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the Whole Brain") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Whole_AvgSubjects[Beta_Whole_AvgSubjects$visitno < 2  & Beta_Whole_AvgSubjects$Beta_Event_Duration < BetaCutoff,], Beta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Alpha
  
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  colnames(Alpha)[1] <- "idvalues"
  
  Alpha_Age <- merge(Alpha, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Alpha_Age$log_Alpha_Power <- log1p(Alpha_Age$Alpha_Trial_Power)
  Alpha_Age$inverseAge <- 1/Alpha_Age$age
  
  Alpha_Age_Channel <- merge(Alpha_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  

  Alpha_Whole_AvgSubjects <- aggregate(. ~ idvalues, Alpha_Age_Channel, mean)
  
  # Trial Power
  print(ggplot(Alpha_Whole_AvgSubjects[Alpha_Whole_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Alpha_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the Whole Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Whole_AvgSubjects[Alpha_Whole_AvgSubjects$visitno < 2 ,], log_Alpha_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Whole_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Whole_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Whole_AvgSubjects[Alpha_Whole_AvgSubjects$visitno < 2 & Alpha_Whole_AvgSubjects$Alpha_Event_Number < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the Whole Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Whole_AvgSubjects[Alpha_Whole_AvgSubjects$visitno < 2  & Alpha_Whole_AvgSubjects$Alpha_Event_Number < AlphaCutoff,], Alpha_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Whole_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Whole_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Whole_AvgSubjects[Alpha_Whole_AvgSubjects$visitno < 2 & Alpha_Whole_AvgSubjects$Alpha_Event_Duration < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the Whole Brain") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Whole_AvgSubjects[Alpha_Whole_AvgSubjects$visitno < 2  & Alpha_Whole_AvgSubjects$Alpha_Event_Duration < AlphaCutoff,], Alpha_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Theta
  
  Theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  colnames(Theta)[1] <- "idvalues"
  
  Theta_Age <- merge(Theta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Theta_Age$log_Theta_Power <- log1p(Theta_Age$Theta_Trial_Power)
  Theta_Age$inverseAge <- 1/Theta_Age$age
  
  Theta_Age_Channel <- merge(Theta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Theta_Whole_AvgSubjects <- aggregate(. ~ idvalues, Theta_Age_Channel, mean)
  
  # Trial Power
  print(ggplot(Theta_Whole_AvgSubjects[Theta_Whole_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Theta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the Whole Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Whole_AvgSubjects[Theta_Whole_AvgSubjects$visitno < 2 ,], log_Theta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Whole_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Whole_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Whole_AvgSubjects[Theta_Whole_AvgSubjects$visitno < 2 & Theta_Whole_AvgSubjects$Theta_Event_Number < ThetaCutoff ,], aes(x = age, y = Theta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the Whole Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Whole_AvgSubjects[Theta_Whole_AvgSubjects$visitno < 2  & Theta_Whole_AvgSubjects$Theta_Event_Number < ThetaCutoff,], Theta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Whole_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Whole_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Whole_AvgSubjects[Theta_Whole_AvgSubjects$visitno < 2 & Theta_Whole_AvgSubjects$Theta_Event_Duration < ThetaCutoff ,], aes(x = age, y = Theta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the Whole Brain") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Whole_AvgSubjects[Theta_Whole_AvgSubjects$visitno < 2  & Theta_Whole_AvgSubjects$Theta_Event_Duration < ThetaCutoff,], Theta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
}

DLPFC_analysis_DelayResting <- function () {
  
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
  
  
  Gamma_DLPFC <- Gamma_DelayRest %>% filter(Label == "'F3'")
  Gamma_DRPFC <- Gamma_DelayRest %>% filter(Label == "'F4'")
  
  Gamma_DLPFC_AvgSubjects <- aggregate(. ~ idvalues, Gamma_DLPFC, mean)
  
  # Trial Power
  avgGamma <- aggregate(log_Gamma_Power ~ visitno , Gamma_DLPFC, mean)
  sdGamma <- aggregate(log_Gamma_Power ~ visitno , Gamma_DLPFC, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_DLPFC[Gamma_DLPFC$visitno < 2 & Gamma_DLPFC$log_Gamma_Power < gammaCutoff,], aes(x = age, y = log_Gamma_Power, group = Task, color = Task)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Gamma Burst Power in the DLPFC") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_DLPFC[Gamma_DLPFC$visitno < 2 & Gamma_DLPFC$log_Gamma_Power < gammaCutoff & Gamma_DLPFC$Task == 'Delay' ,], log_Gamma_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_DLPFC, mean)
  sdGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_DLPFC, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_DLPFC[Gamma_DLPFC$visitno < 2 & Gamma_DLPFC$Gamma_Event_Number < gammaCutoff ,], aes(x = age, y = Gamma_Event_Number, group = Task, color = Task)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the DLPFC") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_DLPFC[Gamma_DLPFC$visitno < 2  & Gamma_DLPFC$Gamma_Event_Number < gammaCutoff & Gamma_DLPFC$Task == 'Delay',], Gamma_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_DLPFC, mean)
  sdGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_DLPFC, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_DLPFC[Gamma_DLPFC$visitno < 2 & Gamma_DLPFC$Gamma_Event_Duration < gammaCutoff ,], aes(x = age, y = Gamma_Event_Duration, group = Task, color = Task)) +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the DLPFC") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_DLPFC[Gamma_DLPFC$visitno < 2  & Gamma_DLPFC$Gamma_Event_Duration < gammaCutoff & Gamma_DLPFC$Task == 'Delay',], Gamma_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  
  colnames(Beta)[1] <- "idvalues"
  
  Beta_Age <- merge(Beta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Beta_Age$log_Beta_Power <- log1p(Beta_Age$Beta_Trial_Power)
  Beta_Age$inverseAge <- 1/Beta_Age$age
  
  Beta_Age_Channel <- merge(Beta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Beta_DLPFC <- Beta_Age_Channel %>% filter(Label == "'F3'")
  Beta_DRPFC <- Beta_Age_Channel %>% filter(Label == "'F4'")
  
  Beta_DLPFC_AvgSubjects <- aggregate(. ~ idvalues, Beta_DLPFC, mean)
  
  # Trial Power
  print(ggplot(Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Beta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the DLPFC") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 ,], log_Beta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_DLPFC_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_DLPFC_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 & Beta_DLPFC_AvgSubjects$Beta_Event_Number < BetaCutoff ,], aes(x = age, y = Beta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the DLPFC") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2  & Beta_DLPFC_AvgSubjects$Beta_Event_Number < BetaCutoff,], Beta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_DLPFC_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_DLPFC_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2 & Beta_DLPFC_AvgSubjects$Beta_Event_Duration < BetaCutoff ,], aes(x = age, y = Beta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the DLPFC") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_DLPFC_AvgSubjects[Beta_DLPFC_AvgSubjects$visitno < 2  & Beta_DLPFC_AvgSubjects$Beta_Event_Duration < BetaCutoff,], Beta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Theta
  
  Theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  colnames(Theta)[1] <- "idvalues"
  
  Theta_Age <- merge(Theta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Theta_Age$log_Theta_Power <- log1p(Theta_Age$Theta_Trial_Power)
  Theta_Age$inverseAge <- 1/Theta_Age$age
  
  Theta_Age_Channel <- merge(Theta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Theta_DLPFC <- Theta_Age_Channel %>% filter(Label == "'F3'")
  Theta_DRPFC <- Theta_Age_Channel %>% filter(Label == "'F4'")
  
  Theta_DLPFC_AvgSubjects <- aggregate(. ~ idvalues, Theta_DLPFC, mean)
  
  # Trial Power
  print(ggplot(Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Theta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the DLPFC") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 ,], log_Theta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_DLPFC_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_DLPFC_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 & Theta_DLPFC_AvgSubjects$Theta_Event_Number < ThetaCutoff ,], aes(x = age, y = Theta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the DLPFC") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2  & Theta_DLPFC_AvgSubjects$Theta_Event_Number < ThetaCutoff,], Theta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_DLPFC_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_DLPFC_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2 & Theta_DLPFC_AvgSubjects$Theta_Event_Duration < ThetaCutoff ,], aes(x = age, y = Theta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the DLPFC") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_DLPFC_AvgSubjects[Theta_DLPFC_AvgSubjects$visitno < 2  & Theta_DLPFC_AvgSubjects$Theta_Event_Duration < ThetaCutoff,], Theta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  #Alpha
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  colnames(Alpha)[1] <- "idvalues"
  
  Alpha_Age <- merge(Alpha, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Alpha_Age$log_Alpha_Power <- log1p(Alpha_Age$Alpha_Trial_Power)
  Alpha_Age$inverseAge <- 1/Alpha_Age$age
  
  Alpha_Age_Channel <- merge(Alpha_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Alpha_DLPFC <- Alpha_Age_Channel %>% filter(Label == "'F3'")
  Alpha_DRPFC <- Alpha_Age_Channel %>% filter(Label == "'F4'")
  
  Alpha_DLPFC_AvgSubjects <- aggregate(. ~ idvalues, Alpha_DLPFC, mean)
  
  # Trial Power
  print(ggplot(Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Alpha_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the DLPFC") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 ,], log_Alpha_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_DLPFC_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_DLPFC_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 & Alpha_DLPFC_AvgSubjects$Alpha_Event_Number < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the DLPFC") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2  & Alpha_DLPFC_AvgSubjects$Alpha_Event_Number < AlphaCutoff,], Alpha_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_DLPFC_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_DLPFC_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2 & Alpha_DLPFC_AvgSubjects$Alpha_Event_Duration < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the DLPFC") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_DLPFC_AvgSubjects[Alpha_DLPFC_AvgSubjects$visitno < 2  & Alpha_DLPFC_AvgSubjects$Alpha_Event_Duration < AlphaCutoff,], Alpha_Event_Duration ~ inverseAge)
  print(anova(lm.model))
}

OccipitalCortex_Analysis_DelayResting <- function () { 
  
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
  
  
  Gamma_Occipital <- filter(Gamma_DelayRest, str_detect(Gamma_DelayRest$Label, "O"))
  
  Gamma_Occipital_AvgSubjects <- aggregate(. ~ idvalues, Gamma_Occipital, mean)
  
  # Trial Power
  avgGamma <- aggregate(log_Gamma_Power ~ visitno , Gamma_Occipital, mean)
  sdGamma <- aggregate(log_Gamma_Power ~ visitno , Gamma_Occipital, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_Occipital[Gamma_Occipital$visitno < 2 & Gamma_Occipital$log_Gamma_Power < gammaCutoff,], aes(x = age, y = log_Gamma_Power, group = Task, color = Task))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Gamma Burst Power in the Occipital Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Occipital[Gamma_Occipital$visitno < 2 & Gamma_Occipital$log_Gamma_Power < gammaCutoff & Gamma_Occipital$Task == 'Rest' ,], log_Gamma_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Occipital, mean)
  sdGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Occipital, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_Occipital[Gamma_Occipital$visitno < 2 & Gamma_Occipital$Gamma_Event_Number < gammaCutoff ,], aes(x = age, y = Gamma_Event_Number, group = Task, color = Task)) +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Occipital[Gamma_Occipital$visitno < 2  & Gamma_Occipital$Gamma_Event_Number < gammaCutoff & Gamma_Occipital$Task == 'Rest',], Gamma_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Occipital, mean)
  sdGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Occipital, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_Occipital[Gamma_Occipital$visitno < 2 & Gamma_Occipital$Gamma_Event_Duration < gammaCutoff ,], aes(x = age, y = Gamma_Event_Duration, group = Task, color = Task)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Occipital[Gamma_Occipital$visitno < 2  & Gamma_Occipital$Gamma_Event_Duration < gammaCutoff & Gamma_Occipital$Task =='Rest',], Gamma_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  # Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  
  colnames(Beta)[1] <- "idvalues"
  
  Beta_Age <- merge(Beta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Beta_Age$log_Beta_Power <- log1p(Beta_Age$Beta_Trial_Power)
  Beta_Age$inverseAge <- 1/Beta_Age$age
  
  Beta_Age_Channel <- merge(Beta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Beta_Occipital <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "O"))
  
  Beta_Occipital_AvgSubjects <- aggregate(. ~ idvalues, Beta_Occipital, mean)
  
  # Trial Power
  print(ggplot(Beta_Occipital_AvgSubjects[Beta_Occipital_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Beta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the Occipital Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Occipital_AvgSubjects[Beta_Occipital_AvgSubjects$visitno < 2 ,], log_Beta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Occipital_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Occipital_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Occipital_AvgSubjects[Beta_Occipital_AvgSubjects$visitno < 2 & Beta_Occipital_AvgSubjects$Beta_Event_Number < BetaCutoff ,], aes(x = age, y = Beta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Occipital_AvgSubjects[Beta_Occipital_AvgSubjects$visitno < 2  & Beta_Occipital_AvgSubjects$Beta_Event_Number < BetaCutoff,], Beta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Occipital_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Occipital_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Occipital_AvgSubjects[Beta_Occipital_AvgSubjects$visitno < 2 & Beta_Occipital_AvgSubjects$Beta_Event_Duration < BetaCutoff ,], aes(x = age, y = Beta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Occipital_AvgSubjects[Beta_Occipital_AvgSubjects$visitno < 2  & Beta_Occipital_AvgSubjects$Beta_Event_Duration < BetaCutoff,], Beta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Alpha
  
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  colnames(Alpha)[1] <- "idvalues"
  
  Alpha_Age <- merge(Alpha, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Alpha_Age$log_Alpha_Power <- log1p(Alpha_Age$Alpha_Trial_Power)
  Alpha_Age$inverseAge <- 1/Alpha_Age$age
  
  Alpha_Age_Channel <- merge(Alpha_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Alpha_Occipital <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "O"))
  
  Alpha_Occipital_AvgSubjects <- aggregate(. ~ idvalues, Alpha_Occipital, mean)
  
  # Trial Power
  print(ggplot(Alpha_Occipital_AvgSubjects[Alpha_Occipital_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Alpha_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the Occipital Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Occipital_AvgSubjects[Alpha_Occipital_AvgSubjects$visitno < 2 ,], log_Alpha_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Occipital_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Occipital_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Occipital_AvgSubjects[Alpha_Occipital_AvgSubjects$visitno < 2 & Alpha_Occipital_AvgSubjects$Alpha_Event_Number < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Occipital_AvgSubjects[Alpha_Occipital_AvgSubjects$visitno < 2  & Alpha_Occipital_AvgSubjects$Alpha_Event_Number < AlphaCutoff,], Alpha_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Occipital_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Occipital_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Occipital_AvgSubjects[Alpha_Occipital_AvgSubjects$visitno < 2 & Alpha_Occipital_AvgSubjects$Alpha_Event_Duration < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Occipital_AvgSubjects[Alpha_Occipital_AvgSubjects$visitno < 2  & Alpha_Occipital_AvgSubjects$Alpha_Event_Duration < AlphaCutoff,], Alpha_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Theta
  
  Theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  colnames(Theta)[1] <- "idvalues"
  
  Theta_Age <- merge(Theta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Theta_Age$log_Theta_Power <- log1p(Theta_Age$Theta_Trial_Power)
  Theta_Age$inverseAge <- 1/Theta_Age$age
  
  Theta_Age_Channel <- merge(Theta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Theta_Occipital <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "O"))
  
  Theta_Occipital_AvgSubjects <- aggregate(. ~ idvalues, Theta_Occipital, mean)
  
  # Trial Power
  print(ggplot(Theta_Occipital_AvgSubjects[Theta_Occipital_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Theta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the Occipital Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Occipital_AvgSubjects[Theta_Occipital_AvgSubjects$visitno < 2 ,], log_Theta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Occipital_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Occipital_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Occipital_AvgSubjects[Theta_Occipital_AvgSubjects$visitno < 2 & Theta_Occipital_AvgSubjects$Theta_Event_Number < ThetaCutoff ,], aes(x = age, y = Theta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Occipital_AvgSubjects[Theta_Occipital_AvgSubjects$visitno < 2  & Theta_Occipital_AvgSubjects$Theta_Event_Number < ThetaCutoff,], Theta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Occipital_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Occipital_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Occipital_AvgSubjects[Theta_Occipital_AvgSubjects$visitno < 2 & Theta_Occipital_AvgSubjects$Theta_Event_Duration < ThetaCutoff ,], aes(x = age, y = Theta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the Occipital Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Occipital_AvgSubjects[Theta_Occipital_AvgSubjects$visitno < 2  & Theta_Occipital_AvgSubjects$Theta_Event_Duration < ThetaCutoff,], Theta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
}


FrontalCortex_Analysis_DelayResting <- function () { 
  
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
  
  
  Gamma_Frontal <- filter(Gamma_DelayRest, str_detect(Gamma_DelayRest$Label, "F"))
  
  Gamma_Frontal_AvgSubjects <- aggregate(. ~ idvalues, Gamma_Frontal, mean)
  
  # Trial Power
  avgGamma <- aggregate(log_Gamma_Power ~ visitno , Gamma_Frontal, mean)
  sdGamma <- aggregate(log_Gamma_Power ~ visitno , Gamma_Frontal, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_Frontal[Gamma_Frontal$visitno < 2 & Gamma_Frontal$log_Gamma_Power < gammaCutoff,], aes(x = age, y = log_Gamma_Power, group = Task, color = Task))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Gamma Burst Power in the Frontal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Frontal[Gamma_Frontal$visitno < 2 & Gamma_Frontal$log_Gamma_Power < gammaCutoff & Gamma_Frontal$Task == 'Rest' ,], log_Gamma_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Frontal, mean)
  sdGamma <- aggregate(Gamma_Event_Number ~ visitno , Gamma_Frontal, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_Frontal[Gamma_Frontal$visitno < 2 & Gamma_Frontal$Gamma_Event_Number < gammaCutoff ,], aes(x = age, y = Gamma_Event_Number, group = Task, color = Task)) +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Gamma Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Frontal[Gamma_Frontal$visitno < 2  & Gamma_Frontal$Gamma_Event_Number < gammaCutoff & Gamma_Frontal$Task == 'Delay',], Gamma_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Frontal, mean)
  sdGamma <- aggregate(Gamma_Event_Duration ~ visitno , Gamma_Frontal, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  lunaize(ggplot(Gamma_Frontal[Gamma_Frontal$visitno < 2 & Gamma_Frontal$Gamma_Event_Duration < gammaCutoff ,], aes(x = age, y = Gamma_Event_Duration, group = Task, color = Task)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Gamma Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Gamma_Frontal[Gamma_Frontal$visitno < 2  & Gamma_Frontal$Gamma_Event_Duration < gammaCutoff & Gamma_Frontal$Task =='Delay',], Gamma_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  # Beta
  
  Beta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allChannels_TrialLevel.csv')  
  colnames(Beta)[1] <- "idvalues"
  
  Beta_Age <- merge(Beta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Beta_Age$log_Beta_Power <- log1p(Beta_Age$Beta_Trial_Power)
  Beta_Age$inverseAge <- 1/Beta_Age$age
  
  Beta_Age_Channel <- merge(Beta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Beta_Frontal <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "O"))
  
  Beta_Frontal_AvgSubjects <- aggregate(. ~ idvalues, Beta_Frontal, mean)
  
  # Trial Power
  print(ggplot(Beta_Frontal_AvgSubjects[Beta_Frontal_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Beta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Beta Power in the Frontal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Frontal_AvgSubjects[Beta_Frontal_AvgSubjects$visitno < 2 ,], log_Beta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Frontal_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Number ~ visitno , Beta_Frontal_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Frontal_AvgSubjects[Beta_Frontal_AvgSubjects$visitno < 2 & Beta_Frontal_AvgSubjects$Beta_Event_Number < BetaCutoff ,], aes(x = age, y = Beta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Beta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Frontal_AvgSubjects[Beta_Frontal_AvgSubjects$visitno < 2  & Beta_Frontal_AvgSubjects$Beta_Event_Number < BetaCutoff,], Beta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Frontal_AvgSubjects, mean)
  sdBeta <- aggregate(Beta_Event_Duration ~ visitno , Beta_Frontal_AvgSubjects, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  print(ggplot(Beta_Frontal_AvgSubjects[Beta_Frontal_AvgSubjects$visitno < 2 & Beta_Frontal_AvgSubjects$Beta_Event_Duration < BetaCutoff ,], aes(x = age, y = Beta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Beta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Beta_Frontal_AvgSubjects[Beta_Frontal_AvgSubjects$visitno < 2  & Beta_Frontal_AvgSubjects$Beta_Event_Duration < BetaCutoff,], Beta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Alpha
  
  Alpha <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allChannels_TrialLevel.csv')  
  colnames(Alpha)[1] <- "idvalues"
  
  Alpha_Age <- merge(Alpha, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Alpha_Age$log_Alpha_Power <- log1p(Alpha_Age$Alpha_Trial_Power)
  Alpha_Age$inverseAge <- 1/Alpha_Age$age
  
  Alpha_Age_Channel <- merge(Alpha_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Alpha_Frontal <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "O"))
  
  Alpha_Frontal_AvgSubjects <- aggregate(. ~ idvalues, Alpha_Frontal, mean)
  
  # Trial Power
  print(ggplot(Alpha_Frontal_AvgSubjects[Alpha_Frontal_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Alpha_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Alpha Power in the Frontal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Frontal_AvgSubjects[Alpha_Frontal_AvgSubjects$visitno < 2 ,], log_Alpha_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Frontal_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Number ~ visitno , Alpha_Frontal_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Frontal_AvgSubjects[Alpha_Frontal_AvgSubjects$visitno < 2 & Alpha_Frontal_AvgSubjects$Alpha_Event_Number < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Alpha Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Frontal_AvgSubjects[Alpha_Frontal_AvgSubjects$visitno < 2  & Alpha_Frontal_AvgSubjects$Alpha_Event_Number < AlphaCutoff,], Alpha_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Frontal_AvgSubjects, mean)
  sdAlpha <- aggregate(Alpha_Event_Duration ~ visitno , Alpha_Frontal_AvgSubjects, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  print(ggplot(Alpha_Frontal_AvgSubjects[Alpha_Frontal_AvgSubjects$visitno < 2 & Alpha_Frontal_AvgSubjects$Alpha_Event_Duration < AlphaCutoff ,], aes(x = age, y = Alpha_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Alpha Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Alpha_Frontal_AvgSubjects[Alpha_Frontal_AvgSubjects$visitno < 2  & Alpha_Frontal_AvgSubjects$Alpha_Event_Duration < AlphaCutoff,], Alpha_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
  # Theta
  
  Theta <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allChannels_TrialLevel.csv')  
  colnames(Theta)[1] <- "idvalues"
  
  Theta_Age <- merge(Theta, agefile, by = "idvalues", all.x = TRUE, all.y = TRUE)
  Theta_Age$log_Theta_Power <- log1p(Theta_Age$Theta_Trial_Power)
  Theta_Age$inverseAge <- 1/Theta_Age$age
  
  Theta_Age_Channel <- merge(Theta_Age, channelLocations, by = "Channel", all.x = TRUE, all.y = TRUE)
  
  
  Theta_Frontal <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "O"))
  
  Theta_Frontal_AvgSubjects <- aggregate(. ~ idvalues, Theta_Frontal, mean)
  
  # Trial Power
  print(ggplot(Theta_Frontal_AvgSubjects[Theta_Frontal_AvgSubjects$visitno < 2 ,], aes(x = age, y = log_Theta_Power)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Trial Theta Power in the Frontal Cortex") + xlab("Age") + ylab("log(Power)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Frontal_AvgSubjects[Theta_Frontal_AvgSubjects$visitno < 2 ,], log_Theta_Power ~ inverseAge)
  print(anova(lm.model))
  
  # Number of events
  avgTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Frontal_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Number ~ visitno , Theta_Frontal_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Frontal_AvgSubjects[Theta_Frontal_AvgSubjects$visitno < 2 & Theta_Frontal_AvgSubjects$Theta_Event_Number < ThetaCutoff ,], aes(x = age, y = Theta_Event_Number)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Number of Theta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Number of Bursts") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Frontal_AvgSubjects[Theta_Frontal_AvgSubjects$visitno < 2  & Theta_Frontal_AvgSubjects$Theta_Event_Number < ThetaCutoff,], Theta_Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  # Duration of Events
  avgTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Frontal_AvgSubjects, mean)
  sdTheta <- aggregate(Theta_Event_Duration ~ visitno , Theta_Frontal_AvgSubjects, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  print(ggplot(Theta_Frontal_AvgSubjects[Theta_Frontal_AvgSubjects$visitno < 2 & Theta_Frontal_AvgSubjects$Theta_Event_Duration < ThetaCutoff ,], aes(x = age, y = Theta_Event_Duration)) + geom_point()  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Duration of Theta Bursts in the Frontal Cortex") + xlab("Age") + ylab("Duration of Bursts (s)") + theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20, hjust = 0.5)))
  
  lm.model <- lm(data = Theta_Frontal_AvgSubjects[Theta_Frontal_AvgSubjects$visitno < 2  & Theta_Frontal_AvgSubjects$Theta_Event_Duration < ThetaCutoff,], Theta_Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
}

