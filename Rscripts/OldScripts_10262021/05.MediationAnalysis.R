
# Mediation Analysis 
library(mediation)
library(lme4)
library(mlma)

#Trial Level 
MediationAnalysis <- function(dataFrame, other_vars) {

  TrialBehavior <- Behavior_TrialLevel_Maria()
  
  alldata_TrialLevel_behavior <- merge(dataFrame, TrialBehavior, by = c("Subject", "Trial"))
  
  alldata_TrialLevel_behavior$Subject <- as.factor(alldata_TrialLevel_behavior$Subject)
  
  # alldata_TrialLevel_behavior_impute <- alldata_TrialLevel_behavior %>% group_by(Subject) %>% VIM::kNN(alldata_TrialLevel_behavior, variable = colnames(alldata_TrialLevel_behavior), k = 10, impNA = T, imp_var = F, imp_suffix = F)
  # 
always_want <- c("Subject","age.x")

#other_vars is defined in rmarkdown 
mediationMatrix <- alldata_TrialLevel_behavior[,c(always_want, other_vars)]
mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]

# xvar <- sprintf("%s", other_vars[1])
# yvar <- sprintf("%s", other_vars[2])
# plotDF <- aggregate(. ~ Subject, data = mediationMatrix, mean)
# plotDF$logabsPE <- log(plotDF$absPositionError)
# 
# avgy <- mean(plotDF$logabsPE)
# sdy <- sd(plotDF$logabsPE)
# Cutoff <- 3* sdy + avgy
# 
# avgx <- mean(plotDF$Gamma_Event_Number)
# sdx <- sd(plotDF$Gamma_Event_Number)
# xCutoff <- 3* sdx + avgx
# xCutofflow <- avgx - 3* sdx
# 
# #plotting one measure against the other 
# lunaize(ggplot(data = mediationMatrix[], aes(x = Gamma_Event_Number, y = absPositionError)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
# 
# 
#  lunaize(ggplot(data = plotDF[plotDF$logabsPE < Cutoff & plotDF$Gamma_Event_Number <xCutoff & plotDF$Gamma_Event_Number > xCutofflow,], aes(x = Gamma_Event_Number, y = logabsPE)) + geom_point() + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
# 
# lm.model <- lm(data = mediationMatrix[], absPositionError ~ Gamma_Event_Number)
# print(anova(lm.model))

# random intercept model 
# modelFormula <- as.formula(sprintf("%s ~ 1 + (1|Subject)", other_vars[2]))
# model.00 <- lmer(modelFormula, REML = F, data = mediationMatrix)
# print(car::Anova(model.00))
# print(summary(model.00))
# confint(model.00)


#the effect of age on behavior (c)
model0Formula <- as.formula(sprintf("%s ~ 1 + age.x + (1|Subject)", other_vars[2]))
model.0 <- lmer(model0Formula, REML = F, data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
#anova(model.00, model.0) # significant improvement in fit when you include age 


#the effect of age on eeg activity  (a)
modelMformula <- as.formula(sprintf("%s ~ 1 + age.x + (1|Subject)", other_vars[1]))
model.M <- lmer(modelMformula, data = mediationMatrix)
print(car::Anova(model.M))
print(summary(model.M))

#the effect of eeg activity on behavior (b)
modelYformula <- as.formula(sprintf("%s ~ 1 + age.x +  %s + (1|Subject)", other_vars[2],  other_vars[1]))
model.Y <- lmer(modelYformula, REML = F, data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))


results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = other_vars[1], boot = FALSE, sims = 1000)
#pvalue <- results$d0.p

#plot(results, group.plots = T)


print(summary(results))

plotRelationshipsAge <- function(z) print(lunaize(ggplot(data = mediationMatrix, aes(x = age.x, y = z)) + stat_smooth(method = 'lm' ,formula='y~I(1/x)') + xlab("Age") + ylab(other_vars[1])))

plotRelationshipsAge(mediationMatrix[,other_vars[1]])

plotRelationships <- function(y,z) print(lunaize(ggplot(data = mediationMatrix, aes(x = y, y = z)) + stat_smooth(method = 'lm') + xlab(other_vars[1]) + ylab(other_vars[2])))

plotRelationships(mediationMatrix[,other_vars[1]], mediationMatrix[,other_vars[2]])


}

MediationAnalysis_SubjectLevel <- function(dataFrame, other_vars) {
  library(class)
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R")
  
  
  avgBehavior <- Behavior_Sublevel_Maria()
  
  alldata_SubLevel_behavior <- merge(dataFrame, avgBehavior, by = "Subject")
  
  alldata_SubLevel_behavior$Subject <- as.factor(alldata_SubLevel_behavior$Subject)

  always_want <- c("Subject","age.x")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- alldata_SubLevel_behavior[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # plotDF <- aggregate(. ~ Subject, data = mediationMatrix, mean)
  # plotDF$logabsPE <- log(plotDF$absPositionError)
  # 
  # avgy <- mean(plotDF$logabsPE)
  # sdy <- sd(plotDF$logabsPE)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$Gamma_Event_Number)
  # sdx <- sd(plotDF$Gamma_Event_Number)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other 
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Gamma_Event_Number, y = absPositionError)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
  # 
  # 
  #  lunaize(ggplot(data = plotDF[plotDF$logabsPE < Cutoff & plotDF$Gamma_Event_Number <xCutoff & plotDF$Gamma_Event_Number > xCutofflow,], aes(x = Gamma_Event_Number, y = logabsPE)) + geom_point() + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError ~ Gamma_Event_Number)
  # print(anova(lm.model))
  
  # random intercept model 
  # modelFormula <- as.formula(sprintf("%s ~ 1 + (1|Subject)", other_vars[2]))
  # model.00 <- lmer(modelFormula, REML = F, data = mediationMatrix)
  # print(car::Anova(model.00))
  # print(summary(model.00))
  # confint(model.00)
  
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ 1 + age.x", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ 1 + age.x", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ 1 + age.x +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = other_vars[1], boot = FALSE, sims = 1000)
  #pvalue <- results$d0.p
  
  #plot(results, group.plots = T)
  
  
  print(summary(results))
  
  
  
  plotRelationshipsAge <- function(z) print(lunaize(ggplot(data = mediationMatrix, aes(x = age.x, y = z)) + stat_smooth(method = 'lm' ,formula='y~I(1/x)') + xlab("Age") + ylab(other_vars[1])))
  
  plotRelationshipsAge(mediationMatrix[,other_vars[1]])
  
  plotRelationships <- function(y,z) print(lunaize(ggplot(data = mediationMatrix, aes(x = y, y = z)) + stat_smooth(method = 'lm') + geom_point() + xlab(other_vars[1]) + ylab(other_vars[2])))
  
  plotRelationships(mediationMatrix[,other_vars[1]], mediationMatrix[,other_vars[2]])
  
  
}


MediationAnalysis_RestingState <- function(dataFrame, other_vars) {
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  agefile$inverseAge <- 1/agefile$age
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  # Behavior outlier detection
  below_2sd <- function(x) abs(x - mean(x,na.rm=T)) < sd(x, na.rm=T)*2
  good <- below_2sd(Behavior[,other_vars[2] ])
  cleanBehavior <-  Behavior[ good,]
  
  AvgBehavior_SubLevel <- aggregate(. ~ Subject, cleanBehavior, mean)
  sdBehavior_Sublevel <- aggregate(. ~ Subject, cleanBehavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  Behavior_SubLevel_age <- merge(Behavior_SubLevel, agefile, by = "Subject")
  Behavior_SubLevel_age <- Behavior_SubLevel_age[Behavior_SubLevel_age$visitno < 2,]
  
  # EEG outlier detection
  below_2sd <- function(x) abs(x - mean(x,na.rm=T)) < sd(x, na.rm=T)*2
  good <- below_2sd(dataFrame[,other_vars[1] ])
  cleanEEG <-  dataFrame[ good,]
  
  avgRest_Sublevel <- aggregate(. ~ Subject, cleanEEG[1:18], mean)
  sdRest_Sublevel <- aggregate(. ~ Subject, cleanEEG[1:18], sd)
  EEG_SubLevel <- merge(avgRest_Sublevel, sdRest_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  cleanData <- merge(EEG_SubLevel, Behavior_SubLevel_age, by = "Subject")

  cleanData <- na.omit(cleanData) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  cleanData <- subset(cleanData, visitno <2)
  always_want <- c("Subject","age")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- cleanData[,c(always_want, other_vars)]
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # plotDF <- aggregate(. ~ Subject, data = mediationMatrix, mean)
  # plotDF$logabsPE <- log(plotDF$absPositionError)
  # 
  # avgy <- mean(plotDF$logabsPE)
  # sdy <- sd(plotDF$logabsPE)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$Gamma_Event_Number)
  # sdx <- sd(plotDF$Gamma_Event_Number)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other 
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Gamma_Event_Number, y = absPositionError)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
  # 
  # 
  #  lunaize(ggplot(data = plotDF[plotDF$logabsPE < Cutoff & plotDF$Gamma_Event_Number <xCutoff & plotDF$Gamma_Event_Number > xCutofflow,], aes(x = Gamma_Event_Number, y = logabsPE)) + geom_point() + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError ~ Gamma_Event_Number)
  # print(anova(lm.model))
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ age", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ age", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  age +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'age', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  print(summary(results))
  
  plotRelationshipsAge <- function(z) print(lunaize(ggplot(data = cleanData, aes(x = age, y = z)) + stat_smooth(method = 'lm' ,formula='y~I(1/x)') + xlab("Age") + ylab(other_vars[1])))
  
  plotRelationshipsAge(cleanData[,other_vars[1]])
  
  plotRelationships <- function(y,z) print(lunaize(ggplot(data = cleanData, aes(x = y, y = z)) + stat_smooth(method = 'lm') + xlab(other_vars[1]) + ylab(other_vars[2])))
  
  plotRelationships(cleanData[,other_vars[1]], cleanData[,other_vars[2]])
  
  
}

MediationAnalysis_bmlmPackage <- function(dataFrame, other_vars) {
  
  library("bmlm")
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  
  alldata_TrialLevel_behavior <- merge(dataFrame, Behavior, by = c("Subject", "Trial"))
  
  #alldata_TrialLevel_behavior$Subject <- as.factor(alldata_TrialLevel_behavior$Subject)
  
  #outlier detection
  below_2sd <- function(x) abs(x - mean(x,na.rm=T)) < sd(x, na.rm=T)*2
  good <- below_2sd(alldata_TrialLevel_behavior[,other_vars[1] ])
  cleanData <-  alldata_TrialLevel_behavior[ good,]
  
  # note! s.d. for second pass will have first pass outliers excluded!
  good <- below_2sd(cleanData[,other_vars[2] ])
  cleanData <-  cleanData[ good,]
  
  
  cleanData <- na.omit(cleanData) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  #cleanData <- subset(test_df, visitno <2)
  always_want <- c("Subject","Trial","age", "Delay")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- cleanData[,c(always_want, other_vars)]
  
  # creates a within- person component for abs(position error) and then gamma event number 
  mediationMatrix_isolated <- isolate(d = mediationMatrix, by = "Subject", value = c("absPositionError", "Gamma.Event_Number", "age"))
mediationMatrix_isolated$id <- as.integer((mediationMatrix_isolated$Subject))
  
  
  fit <- mlm(d = mediationMatrix_isolated, id = "id", x = "age", m = "Gamma.Event_Number_cw", y = "absPositionError_cw", iter = 2000, cores = 4)
  
  mlm_path_plot(fit, xlab = "Age", mlab = "Gamma Event Number", ylab = "abs(Position Error")
  
  mlm_pars_plot(fit, type = "coef", pars = "u_cp", level = .80) 
  
  
}


MediationAnalysis_inverseAge <- function(other_vars) {
  
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
  completeSubjects$inverseAge <- 1/completeSubjects$age
  
  test_df <- na.omit(completeSubjects) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno <2)
  always_want <- c("Subject","Trial","inverseAge", "Delay")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # plotDF <- aggregate(. ~ Subject, data = mediationMatrix, mean)
  # plotDF$logabsPE <- log(plotDF$absPositionError)
  # 
  # avgy <- mean(plotDF$logabsPE)
  # sdy <- sd(plotDF$logabsPE)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$Gamma_Event_Number)
  # sdx <- sd(plotDF$Gamma_Event_Number)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other 
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Gamma_Event_Number, y = absPositionError)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
  # 
  # 
  #  lunaize(ggplot(data = plotDF[plotDF$logabsPE < Cutoff & plotDF$Gamma_Event_Number <xCutoff & plotDF$Gamma_Event_Number > xCutofflow,], aes(x = Gamma_Event_Number, y = logabsPE)) + geom_point() + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Gamma Bursts") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError ~ Gamma_Event_Number)
  # print(anova(lm.model))
  
  #the effect of inverse age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge + (1|Subject)", other_vars[2]))
  model.0 <- lmer(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of inverse age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge + (1|Subject)", other_vars[1]))
  model.M <- lmer(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s + (1|Subject)", other_vars[2],  other_vars[1]))
  model.Y <- lmer(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
}

MediationAnalysisFixation <- function(other_vars) {
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200821.csv')
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  
  FixGammaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_Subs_FIX.csv')
  FixBetaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_Subs_dataFix.csv')
  FixAlphaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_Subs_datafix.csv')
  FixThetaTrials <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_Subs_dataFix.csv')
  
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues
  
  completeSubjects <-  merge(Behavior, FixGammaTrials, by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, FixBetaTrials, by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, FixAlphaTrials, by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, FixThetaTrials, by = c("Subject", "Trial"))
  completeSubjects <-  merge(completeSubjects, agefile, by = "Subject")
  completeSubjects$absPositionError <- abs(completeSubjects$PositionError)
  completeSubjects$Gamma_log1p_Trial_Power <- log1p(completeSubjects$Gamma_Trial_Power)
  completeSubjects$Beta_log1p_Trial_Power <- log1p(completeSubjects$Beta_Trial_Power)
  completeSubjects$Alpha_log1p_Trial_Power <- log1p(completeSubjects$Alpha_Trial_Power)
  completeSubjects$Theta_log1p_Trial_Power <- log1p(completeSubjects$Theta_Trial_Power)
  completeSubjects$inverseAge <- 1/completeSubjects$age
  
  
  
  completeSubjects$Subject <- as.factor(completeSubjects$Subject)
  
  test_df <- na.omit(completeSubjects) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno <2)

  always_want <- c("Subject","Trial","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # plotDF <- aggregate(. ~ Subject, data = mediationMatrix, mean)

  # avgy <- mean(mediationMatrix$mgsLatency)
  # sdy <- sd(mediationMatrix$mgsLatency)
  # Cutoff <- 2.5* sdy + avgy
  # 
  # avgx <- mean(mediationMatrix$Alpha_Event_Number)
  # sdx <- sd(mediationMatrix$Alpha_Event_Number)
  # xCutoff <- 2.5* sdx + avgx
  # xCutofflow <- avgx - 2.5* sdx
  # 
  # lunaize(ggplot(data = mediationMatrix[mediationMatrix$Alpha_Event_Number <= 3,], aes(x = Alpha_Event_Number, y = mgsLatency))+ stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Alpha Bursts") + ylab("mgsLatency (sec)"))
  # 
  # # 
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Alpha_Event_Number, y = mgsLatency)) + geom_point() + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Number of Alpha Bursts") + ylab("mgsLatency (sec)"))
  # # 
  # lm.model <- lm(data = mediationMatrix[mediationMatrix$Alpha_Event_Number], mgsLatency ~ Alpha_Event_Number)
  # print(anova(lm.model))
  # 
  model0Formula <- as.formula(sprintf("%s ~ inverseAge + (1|Subject)", other_vars[2]))
  model.0 <- lmer(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
 
  
  modelMformula <- as.formula(sprintf("%s ~ inverseAge + (1|Subject)", other_vars[1]))
  model.M <- lmer(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
 
  
  modelYformula <- as.formula(sprintf("%s ~ inverseAge + %s + (1|Subject)", other_vars[2],  other_vars[1]))
  model.Y <- lmer(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
 
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  summary(results)
  
 
 # pvalue <- med$d0.p
  
  #return(pvalue)  
  #return(mediationEffect)

}


#Subject level 

MediationAnalysis_SubjectLevel_Gamma <- function(dataFrame, other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  delayGamma_behav_SubLevel <- merge(dataFrame, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  delayGamma_behav_SubLevel <- merge(delayGamma_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  delayGamma_behav_SubLevel$Subject <- as.factor(delayGamma_behav_SubLevel$Subject)

  test_df <- na.omit(delayGamma_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge <- 1/test_df$age
  
  always_want <- c("Subject","age")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]

  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ age", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  summary(model.0)
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ age", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  summary(model.M)
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  age +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  summary(model.Y)
  
  results <- mediate(model.M, model.Y, treat = 'age', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Delay)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Delay)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Delay, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Delay)
  # print(anova(lm.model))
}

MediationAnalysis_SubjectLevel_Beta <- function(other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  allBeta_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/allBeta_3_4_20200821.csv')
 
  delayVars <- grep("Delay|Subject", names(allBeta_SubLevel), value = TRUE)
  delayVars <- delayVars[!grepl("Trial.x|Trial.y", delayVars)]
  delayBeta_SubLevel <- allBeta_SubLevel[,delayVars]
  
  delayBeta_behav_SubLevel <- merge(delayBeta_SubLevel, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  delayBeta_behav_SubLevel <- merge(delayBeta_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  delayBeta_behav_SubLevel$Subject <- as.factor(delayBeta_behav_SubLevel$Subject)
  
  test_df <- na.omit(delayBeta_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge <- 1/test_df$age
  
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Delay)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Delay)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Delay, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Delay)
  # print(anova(lm.model))
}

MediationAnalysis_SubjectLevel_Alpha <- function(other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  allAlpha_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/allAlpha_3_4_20200821.csv')

  delayVars <- grep("Delay|Subject", names(allAlpha_SubLevel), value = TRUE)
  delayVars <- delayVars[!grepl("Trial.x|Trial.y", delayVars)]
  delayAlpha_SubLevel <- allAlpha_SubLevel[,delayVars]
  
  delayAlpha_behav_SubLevel <- merge(delayAlpha_SubLevel, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  delayAlpha_behav_SubLevel <- merge(delayAlpha_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  delayAlpha_behav_SubLevel$Subject <- as.factor(delayAlpha_behav_SubLevel$Subject)
  
  test_df <- na.omit(delayAlpha_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge <- 1/test_df$age
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Delay)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Delay)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Delay, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Delay)
  # print(anova(lm.model))
}

MediationAnalysis_SubjectLevel_Theta <- function(other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  

  allTheta_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/allTheta_3_4_20200821.csv')
  
  delayVars <- grep("Delay|Subject", names(allTheta_SubLevel), value = TRUE)
  delayVars <- delayVars[!grepl("Trial.x|Trial.y", delayVars)]
  delayTheta_SubLevel <- allTheta_SubLevel[,delayVars]
  
  delayTheta_behav_SubLevel <- merge(delayTheta_SubLevel, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  delayTheta_behav_SubLevel <- merge(delayTheta_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  delayTheta_behav_SubLevel$Subject <- as.factor(delayTheta_behav_SubLevel$Subject)
  
  test_df <- na.omit(delayTheta_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge <- 1/test_df$age
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Delay)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Delay)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Delay, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Delay)
  # print(anova(lm.model))
}



MediationAnalysis_SubjectLevel_GammaFix <- function(other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  allGamma_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/allGamma_3_4_20200821.csv')
  allBeta_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/allBeta_3_4_20200821.csv')
  allAlpha_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/allAlpha_3_4_20200821.csv')
  allTheta_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/allTheta_3_4_20200821.csv')
  
  FixVars <- grep("Fix|Subject", names(allGamma_SubLevel), value = TRUE)
  FixVars <- FixVars[!grepl("Trial.x|Trial.y", FixVars)]
  FixGamma_SubLevel <- allGamma_SubLevel[,FixVars]
  
  FixGamma_behav_SubLevel <- merge(FixGamma_SubLevel, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  FixGamma_behav_SubLevel <- merge(FixGamma_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  FixGamma_behav_SubLevel$Subject <- as.factor(FixGamma_behav_SubLevel$Subject)
  
  test_df <- na.omit(FixGamma_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge <- 1/test_df$age
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Fix)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Fix)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Fix, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Fix)
  # print(anova(lm.model))
}

MediationAnalysis_SubjectLevel_BetaFix <- function(other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  allBeta_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/allBeta_3_4_20200821.csv')
  
  FixVars <- grep("Fix|Subject", names(allBeta_SubLevel), value = TRUE)
  FixVars <- FixVars[!grepl("Trial.x|Trial.y", FixVars)]
  FixBeta_SubLevel <- allBeta_SubLevel[,FixVars]
  
  FixBeta_behav_SubLevel <- merge(FixBeta_SubLevel, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  FixBeta_behav_SubLevel <- merge(FixBeta_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  FixBeta_behav_SubLevel$Subject <- as.factor(FixBeta_behav_SubLevel$Subject)
  
  test_df <- na.omit(FixBeta_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge <- 1/test_df$age
  
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Fix)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Fix)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Fix, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Fix)
  # print(anova(lm.model))
}

MediationAnalysis_SubjectLevel_AlphaFix <- function(other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  allAlpha_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/allAlpha_3_4_20200821.csv')
  
  FixVars <- grep("Fix|Subject", names(allAlpha_SubLevel), value = TRUE)
  FixVars <- FixVars[!grepl("Trial.x|Trial.y", FixVars)]
  FixAlpha_SubLevel <- allAlpha_SubLevel[,FixVars]
  
  FixAlpha_behav_SubLevel <- merge(FixAlpha_SubLevel, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  FixAlpha_behav_SubLevel <- merge(FixAlpha_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  FixAlpha_behav_SubLevel$Subject <- as.factor(FixAlpha_behav_SubLevel$Subject)
  
  test_df <- na.omit(FixAlpha_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge <- 1/test_df$age
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Fix)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Fix)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Fix, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Fix)
  # print(anova(lm.model))
}

MediationAnalysis_SubjectLevel_ThetaFix <- function(other_vars) {
  
  Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
  agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
  agefile$Subject <- agefile$idvalues 
  
  #behavioral data contains each trial 
  Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
  Behavior$absPositionError <- abs(Behavior$PositionError)
  AvgBehavior_SubLevel <- aggregate(.~ Subject, Behavior, mean)
  sdBehavior_Sublevel <- aggregate(.~Subject, Behavior, sd)
  
  Behavior_SubLevel <- merge(AvgBehavior_SubLevel, sdBehavior_Sublevel, by = "Subject", suffix = c("", "_Variability"))
  
  
  allTheta_SubLevel <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/allTheta_3_4_20200821.csv')
  
  FixVars <- grep("Fix|Subject", names(allTheta_SubLevel), value = TRUE)
  FixVars <- FixVars[!grepl("Trial.x|Trial.y", FixVars)]
  FixTheta_SubLevel <- allTheta_SubLevel[,FixVars]
  
  FixTheta_behav_SubLevel <- merge(FixTheta_SubLevel, Behavior_SubLevel, by = "Subject", all.x = TRUE, all.y = TRUE)
  FixTheta_behav_SubLevel <- merge(FixTheta_behav_SubLevel, agefile, by = "Subject", all.x = TRUE, all.y = TRUE)
  
  
  FixTheta_behav_SubLevel$Subject <- as.factor(FixTheta_behav_SubLevel$Subject)
  
  test_df <- na.omit(FixTheta_behav_SubLevel) #need to get rid of the rows with NAN values (some subjects have NAN for position error)
  test_df <- subset(test_df, visitno < 2)
  test_df$inverseAge<- 1/test_df$age
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- test_df[,c(always_want, other_vars)]
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~  inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  
  results <- mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  pvalue <- results$d0.p
  
  summary(results)
  
  
  #if you want to graph one variable against the other 
  
  # xvar <- sprintf("%s", other_vars[1])
  # yvar <- sprintf("%s", other_vars[2])
  # 
  # plotDF <- mediationMatrix
  # 
  # avgy <- mean(plotDF$absPositionError_Variability)
  # sdy <- sd(plotDF$absPositionError_Variability)
  # Cutoff <- 3* sdy + avgy
  # 
  # avgx <- mean(plotDF$log1p_Trial_Power_Variability_Fix)
  # sdx <- sd(plotDF$log1p_Trial_Power_Variability_Fix)
  # xCutoff <- 3* sdx + avgx
  # xCutofflow <- avgx - 3* sdx
  # 
  # #plotting one measure against the other
  # lunaize(ggplot(data = mediationMatrix[], aes(x = Event_Number_Variability_Fix, y = absPositionError_Variability)) + stat_smooth(method = "lm") + ggtitle(" ") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("abs(Position Error (degs))"))
  # 
  # lm.model <- lm(data = mediationMatrix[], absPositionError_Variability ~ Event_Number_Variability_Fix)
  # print(anova(lm.model))
}
