

delayRestCorrelation <- function() {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  #all data, delay and fix, at the subject level 
  alldf <- Only_Take_One_Delay_Bin_and_Fix()
  alldata <- alldf$alldata
  
  #only delay variables, trial level
  alldf_trialLevel <- Only_Take_One_Delay_Bin_TrialLevel()
  alldata_trialLevel <- alldf_trialLevel$alldata_TrialLevel
  alldata_trialLevel$Task <- "Delay"
  names(alldata_trialLevel) <- sub("Gamma_","",names(alldata_trialLevel))
  names(alldata_trialLevel) <- sub("Beta_","",names(alldata_trialLevel))
  names(alldata_trialLevel) <- sub("Theta_","",names(alldata_trialLevel))
  names(alldata_trialLevel) <- sub("Alpha_","",names(alldata_trialLevel))
  
  #only the delay variables at the subject level 
  delayOnlyVars <- grep("Delay|Subject|age|idvalues|Group|visitno|inverse", names(alldata), value = TRUE)
  delayOnlyVars <- delayOnlyVars[!grepl("Trial.x|Trial.y|Peak", delayOnlyVars)]
  delayOnly <- alldata[,delayOnlyVars]
  delayOnly <- subset(delayOnly, select = -c(Gamma.Trial_Power_Variability_Delay, Beta.Trial_Power_Variability_Delay, Alpha.Trial_Power_Variability_Delay,Theta.Trial_Power_Variability_Delay))
  
  
  delayOnly$Task <- "Delay"
  names(delayOnly) <- sub("_Delay","",names(delayOnly))
  
  #resting state at the trial level 
  alldf_RS <- Resting_State_Data_TrialLevel()
  alldata_RS <- alldf_RS$alldata_RS_TrialLevel
  alldata_RS$Task <- "Rest"
  
  #resting state at the subject level
  alldf_SubjectLevel_RS <- Resting_State_Data_SubjectLevel()
  alldata_RS_SubLevel <- alldf_SubjectLevel_RS$rest_SubLevel
  alldata_RS_SubLevel$Task <- "Rest"
  alldata_RS_SubLevel <- subset(alldata_RS_SubLevel, select = -Trial)
  
  #merging delay and rest, at the subject level 
  
  delayRest <- merge(delayOnly, alldata_RS_SubLevel[1:30], by = "Subject", suffixes = c("_Delay", "_Rest"))
  delayRest_long <- rbind(delayOnly, alldata_RS_SubLevel)
  
  # Correlation matrix 
  install.packages("corrplot")
  library(corrplot)
  
  res <- cor(delayRest[c(2:29,36:63)])
  round(res, 2)
  
  # Power and Duration Correlation for all variables 
  powerDurationVars <- grep("Duration|log",names(delayRest),value=TRUE) 
  powerDuration <- delayRest[,powerDurationVars]
  res <- cor(powerDuration)
  
  install.packages("Hmisc")
  library(Hmisc)
  res2 <- rcorr(as.matrix(powerDuration))
  pmat <- res2$P
  corrplot(res, type = "upper", tl.cex = 0.5, sig.level = 0.05)
  
  # Power and Duration Correlation for Gamma Only 
  powerDurationVars <- grep("Gamma",names(delayRest),value=TRUE)
  powerDurationVars <- powerDurationVars[!grepl("Number|a.Trial_Power", powerDurationVars)]
  
  powerDuration <- delayRest[,powerDurationVars]
  setnames(powerDuration, c("EDD", "TPD", "EDVD", "TPVD","EDR", "TPR", "TPVR", "EDVR"))
  
  res <- cor(powerDuration)
  
  res2 <- rcorr(as.matrix(powerDuration))
  pmat <- res2$P
  
  corrplot(res, type = "upper", tl.cex = 1, method = "number")
}

powerNumberCorrelation_ageGroups <- function() {
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  
  #only delay variables, trial level
  alldf_trialLevel <- Only_Take_One_Delay_Bin_TrialLevel()
  alldata_trialLevel <- alldf_trialLevel$alldata_TrialLevel
  alldata_trialLevel$Task <- "Delay"
  names(alldata_trialLevel) <- sub("Gamma_","",names(alldata_trialLevel))
  names(alldata_trialLevel) <- sub("Beta_","",names(alldata_trialLevel))
  names(alldata_trialLevel) <- sub("Theta_","",names(alldata_trialLevel))
  names(alldata_trialLevel) <- sub("Alpha_","",names(alldata_trialLevel))
  
  #only the delay variables at the subject level 
  alldf_sublevel <- DelayOnly_Sublevel()
  alldata_delayOnly <- alldf_sublevel$alldata_delayOnly
  
  
  # Correlation matrix 
  install.packages("corrplot")
  library(corrplot)
  
  gammaOnlyVars <- grep("Gamma|Subject|age|idvalues|Group|visitno", names(alldata_delayOnly), value = TRUE)
  gammaOnlyVars <- gammaOnlyVars[!grepl("Duration|a.Trial_Power", gammaOnlyVars)]
  gammaOnly <- na.omit(alldata_delayOnly[,gammaOnlyVars])
  
  youngestGroup <- gammaOnly[gammaOnly$Group == 1,]
  secondyoungestGroup <- gammaOnly[gammaOnly$Group == 2,]
  thirdyoungest <- gammaOnly[gammaOnly$Group == 3,]
  oldest <- gammaOnly[gammaOnly$Group == 4,]
  
  youngestGroup <- setnames(youngestGroup[c(2:5)], c("EN", "TP", "ENV", "TPV"))
  
  
  res <- cor(youngestGroup)
  round(res, 2)
  
  library(Hmisc)
  res2 <- rcorr(as.matrix(youngestGroup))
  pmat <- res2$P
  corrplot(res, type = "upper", method = "number", tl.cex = 0.8, sig.level = 0.05)
  flattenCorrMatrix(res2$r, res2$P)
  
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  
}

gammaBetaTrialCorrelations <- function () {
  
  alldf <- Only_Take_One_Delay_Bin_TrialLevel()
  alldata_TrialLevel <- alldf$alldata_TrialLevel_new
  
  # duration 
  truncatedMatrix <- subset(alldata_TrialLevel, select = c("Subject", "Trial", "Gamma.Event_Duration", "Beta.Event_Duration"))
  
  library(plyr)
  subjectCorr <- function(truncatedMatrix) {
    
    return(data.frame(durationCoeff = cor(truncatedMatrix$Gamma.Event_Duration, truncatedMatrix$Beta.Event_Duration, method = "pearson", use = "complete.obs")))
    
  }
  
  durationCoeff <- ddply(truncatedMatrix, .(Subject), subjectCorr)  
  
  
  # number
  truncatedMatrix <- subset(alldata_TrialLevel, select = c("Subject", "Trial", "Gamma.Event_Number", "Beta.Event_Number"))
  
  subjectCorr <- function(truncatedMatrix) {
    
    return(data.frame(numberCoeff = cor(truncatedMatrix$Gamma.Event_Number, truncatedMatrix$Beta.Event_Number, method = "pearson", use = "complete.obs")))
    
  }
  
  numberCoeff <- ddply(truncatedMatrix, .(Subject), subjectCorr)  
  
  
  
  # power
  truncatedMatrix <- subset(alldata_TrialLevel, select = c("Subject", "Trial", "Gamma.log1p_Trial_Power", "Beta.log1p_Trial_Power"))
  
  subjectCorr <- function(truncatedMatrix) {
    
    return(data.frame(powerCoeff = cor(truncatedMatrix$Gamma.log1p_Trial_Power, truncatedMatrix$Beta.log1p_Trial_Power, method = "pearson", use = "complete.obs")))
    
  }
  
  powerCoeff <- ddply(truncatedMatrix, .(Subject), subjectCorr) 
  
  # merge all coeffs
  subjectCoeffAges <- merge(durationCoeff, numberCoeff, by = "Subject") %>% merge(., powerCoeff, by = "Subject")%>% merge(., agefile, by = "Subject")
  subjectCoeffAges$inverseAge <- 1/subjectCoeffAges$age
  
  
  # plotting
  
  lunaize(ggplot(data = subjectCoeffAges[], aes(x = age, y = durationCoeff)) + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')) + ggtitle("Correlation Coefficient \nGamma Power and Beta Power per Subject vs. Age") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30))
  lm.model <- lm(data = subjectCoeffAges[], durationCoeff ~ inverseAge)
  print(anova(lm.model))
  

  
  
  # combine with behavior subject level 
  
  subjectBehavior_new <- Behavior_Sublevel_Maria()
  
  subjectCoeffAgesBehavior <- merge(subjectCoeffAges, subjectBehavior_new, by = "Subject")
  
  lunaize(ggplot(data = subjectCoeffAgesBehavior, aes(x = absPositionError, y = durationCoeff)) + geom_point() + stat_smooth(method = 'lm')) + ggtitle("Correlation Coefficient \nGamma Duration and Beta Duration per Subject vs. Accuracy") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30))
  
  lm.model <- lm(data = subjectCoeffAgesBehavior[], durationCoeff ~ absPositionError)
  print(anova(lm.model))
  
  
  # only look at adults 
  lunaize(ggplot(data = subjectCoeffAgesBehavior[subjectCoeffAgesBehavior$Group == 4,], aes(x = absPositionError, y = durationCoeff)) + geom_point() + stat_smooth(method = 'lm')) + ggtitle("Correlation Coefficient \nGamma Duration and Beta Duration per Subject (Adults) vs. Accuracy") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30))
  lm.model <- lm(data = subjectCoeffAgesBehavior[subjectCoeffAgesBehavior$Group == 4,], durationCoeff ~ absPositionError)
  print(anova(lm.model))
  
}

individualChannelCorrelations <- function () {
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)

  #Gamma
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
  
  
  #Beta
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

  
  allFrontal <- merge(Gamma_Frontal_avgChannel, Beta_Frontal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_Frontal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_Frontal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))

  allFrontal_sd <- merge(Gamma_Frontal_sdChannel, Beta_Frontal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_Frontal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_Frontal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
  

 res <- cor(allFrontal[c(6:17)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 flattenCorrMatrix <- function(cormat, pmat) {
   ut <- upper.tri(cormat)
   data.frame(
     row = rownames(cormat)[row(cormat)[ut]],
     column = rownames(cormat)[col(cormat)[ut]],
     cor  =(cormat)[ut],
     p = pmat[ut]
   )
 }
 
 
library(Hmisc)
  res <- rcorr(as.matrix(allFrontal[c(6:17)]))
 flattenCorrMatrix(res$r, res$P)
 
 
 
 res <- cor(allFrontal_sd[c(5:16)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 res <- rcorr(as.matrix(allFrontal_sd[c(5:16)]))
 flattenCorrMatrix(res$r, res$P)
 
 #Parietal 
 #Gamma
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
 
 
 #Beta
 Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
 Beta_Parietal <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "P"))
 
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
 
 
 
 # Alpha
 Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
 Alpha_Parietal <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "P"))
 
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
 
 
 
 # Theta
 Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
 Theta_Parietal <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "P"))
 
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
 
 
 

 allParietal <- merge(Gamma_Parietal_avgChannel, Beta_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))
 

 allParietal_sd <- merge(Gamma_Parietal_sdChannel, Beta_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
 
 
 res <- cor(allParietal[c(6:17)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 flattenCorrMatrix <- function(cormat, pmat) {
   ut <- upper.tri(cormat)
   data.frame(
     row = rownames(cormat)[row(cormat)[ut]],
     column = rownames(cormat)[col(cormat)[ut]],
     cor  =(cormat)[ut],
     p = pmat[ut]
   )
 }
 
 
 library(Hmisc)
 res <- rcorr(as.matrix(allParietal[c(6:17)]))
 flattenCorrMatrix(res$r, res$P)
 
 
 
 res <- cor(allParietal_sd[c(5:16)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 res <- rcorr(as.matrix(allParietal_sd[c(5:16)]))
 flattenCorrMatrix(res$r, res$P)
 
 
 
 #Occipital 
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
 
 
 #Beta
 Beta_Age_Channel <- individualChannelDF$Beta_Age_Channel_new
 Beta_Occipital <- filter(Beta_Age_Channel, str_detect(Beta_Age_Channel$Label, "O"))
 
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
 
 
 
 # Alpha
 Alpha_Age_Channel <- individualChannelDF$Alpha_Age_Channel_new
 Alpha_Occipital <- filter(Alpha_Age_Channel, str_detect(Alpha_Age_Channel$Label, "O"))
 
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
 
 
 
 # Theta
 Theta_Age_Channel <- individualChannelDF$Theta_Age_Channel_new
 Theta_Occipital <- filter(Theta_Age_Channel, str_detect(Theta_Age_Channel$Label, "O"))
 
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
 
 
 
 
 allOccipital <- merge(Gamma_Occipital_avgChannel, Beta_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Alpha_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge")) %>% merge(., Theta_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"))
 
 
 allOccipital_sd <- merge(Gamma_Occipital_sdChannel, Beta_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Alpha_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge")) %>% merge(., Theta_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge"))
 
 
 res <- cor(allOccipital[c(6:17)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 flattenCorrMatrix <- function(cormat, pmat) {
   ut <- upper.tri(cormat)
   data.frame(
     row = rownames(cormat)[row(cormat)[ut]],
     column = rownames(cormat)[col(cormat)[ut]],
     cor  =(cormat)[ut],
     p = pmat[ut]
   )
 }
 
 
 library(Hmisc)
 res <- rcorr(as.matrix(allOccipital[c(6:17)]))
 flattenCorrMatrix(res$r, res$P)
 
 
 
 res <- cor(allOccipital_sd[c(5:16)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 res <- rcorr(as.matrix(allOccipital_sd[c(5:16)]))
 flattenCorrMatrix(res$r, res$P)
 
 
 #Gamma
 
 allGamma <- merge(Gamma_Frontal_avgChannel, Gamma_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Gamma_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("","_Occipital")) 
 
 names(allGamma) <- c("Subject", "Trial", "age", "Group", "inverseAge", "GENF", "GPF", "GEDF", "GENP", "GPP", "GEDP","GENO", "GPO", "GEDO")
 
 allGamma_sd <- merge(Gamma_Frontal_sdChannel, Gamma_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Gamma_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("", "_Occipital")) 
 
 names(allGamma_sd) <- c("Subject", "age", "Group", "inverseAge", "GENF", "GPF", "GEDF", "GENP", "GPP", "GEDP","GENO", "GPO", "GEDO")
 
 
 res <- cor(allGamma[c(6:14)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 res <- cor(allGamma_sd[c(5:13)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 
 
 #Beta
 
 allBeta <- merge(Beta_Frontal_avgChannel, Beta_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Beta_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("","_Occipital")) 
 
 names(allBeta) <- c("Subject", "Trial", "age", "Group", "inverseAge", "BENF", "BPF", "BEDF", "BENP", "BPP", "BEDP","BENO", "BPO", "BEDO")
 
 allBeta_sd <- merge(Beta_Frontal_sdChannel, Beta_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Beta_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("", "_Occipital")) 
 
 names(allBeta_sd) <- c("Subject", "age", "Group", "inverseAge", "BENF", "BPF", "BEDF", "BENP", "BPP", "BEDP","BENO", "BPO", "BEDO")
 
 
 res <- cor(allBeta[c(6:14)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 res <- cor(allBeta_sd[c(5:13)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 
 
 #Alpha
 
 allAlpha <- merge(Alpha_Frontal_avgChannel, Alpha_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Alpha_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("","_Occipital")) 
 
 names(allAlpha) <- c("Subject", "Trial", "age", "Group", "inverseAge", "AENF", "APF", "AEDF", "AENP", "APP", "AEDP","AENO", "APO", "AEDO")
 
 allAlpha_sd <- merge(Alpha_Frontal_sdChannel, Alpha_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Alpha_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("", "_Occipital")) 
 
 names(allAlpha_sd) <- c("Subject", "age", "Group", "inverseAge", "AENF", "APF", "AEDF", "AENP", "APP", "AEDP","AENO", "APO", "AEDO")
 
 
 res <- cor(allAlpha[c(6:14)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 res <- cor(allAlpha_sd[c(5:13)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 
 #Theta
 
 allTheta <- merge(Theta_Frontal_avgChannel, Theta_Parietal_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Theta_Occipital_avgChannel, by = c("Subject", "Trial", "age", "Group", "inverseAge"), suffixes = c("","_Occipital")) 
 
 names(allTheta) <- c("Subject", "Trial", "age", "Group", "inverseAge",  "TENF", "TPF", "TEDF", "TENP", "TPP", "TEDP","TENO", "TPO", "TEDO")
 
 allTheta_sd <- merge(Theta_Frontal_sdChannel, Theta_Parietal_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("_Frontal", "_Parietal")) %>% merge(., Theta_Occipital_sdChannel, by = c("Subject", "age", "Group", "inverseAge"), suffixes = c("", "_Occipital")) 
 
 names(allTheta_sd) <- c("Subject", "age", "Group", "inverseAge", "TENF", "TPF", "TEDF", "TENP", "TPP", "TEDP","TENO", "TPO", "TEDO")
 
 
 res <- cor(allTheta[c(6:14)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
 res <- cor(allTheta_sd[c(5:13)], method = "pearson", use = "complete.obs")
 corrplot(res)
 
}