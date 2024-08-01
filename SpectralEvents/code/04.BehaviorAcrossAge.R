
TrialPowerBehav <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError, y = log1p_Trial_Power)) + stat_smooth(method='lm', se = T) + ggtitle("Trial Power vs Best Saccade") +theme_classic()))
  
  
  lm.model <- lm(data = delayBehavior[], log1p_Trial_Power ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power*age, data = delayBehavior[])
  print(summary(model))
  
  #age groups
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError, y = log1p_Trial_Power, group = Group, color = Group)) + stat_smooth(method='lm', se = T) + ggtitle("Trial Power vs Best Saccade") +theme_classic()))
  

}

TrialPowerVarBehav <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError, y = log1p_Trial_Power_Variability)) + stat_smooth(method='lm', se = T) + ggtitle("Trial Power Var vs Best Saccade")  +theme_classic()))
  
 
  
  lm.model <- lm(data = delayBehavior[], log1p_Trial_Power_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power_Variability*age, data = delayBehavior[])
  print(summary(model))
  
 
  
}

EventNumberBehav <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError, y = Event_Number))  + stat_smooth(method='loess', se = T) + ggtitle("Event Number vs Best Saccade") +theme_classic()))
  
 
  lm.model <- lm(data = delayBehavior[], Event_Number ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number*age, data = delayBehavior[])
  print(summary(model))
  
  
  
}

EventNumberVarBehav <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError , y = Event_Number_Variability)) + stat_smooth(method='lm', se = T) + ggtitle("Event Number Var vs Best Saccade") +theme_classic()))
  
 
  lm.model <- lm(data = delayBehavior[], Event_Number_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number_Variability*age, data = delayBehavior[])
  print(summary(model))
  
  
  
}

EventDurationBehav <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError , y = Event_Duration )) + stat_smooth(method='lm', se = T) + ggtitle("Event Duration vs Best Saccade") +theme_classic()))
  
  cat(paste("",
            "",
            "Gamma Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[], Event_Duration ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration*age, data = delayBehavior[])
  print(summary(model))
  

  
}

EventDurationVarBehav <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError , y = Event_Duration_Variability)) + stat_smooth(method='lm', se = T) + ggtitle("Event Duration Var vs Best Saccade") +theme_classic()))
  
 
  
  lm.model <- lm(data = delayBehavior[], Event_Duration_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration_Variability*age, data = delayBehavior[])
  print(summary(model))
  

  
}


TrialPowerBehavAllRegions <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError, y = log1p_Trial_Power, group = Region, color = Region)) + stat_smooth(method='lm', se = T) + ggtitle("Trial Power vs Best Saccade") +theme_classic()))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Whole Brain',], log1p_Trial_Power ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power*age, data = delayBehavior[delayBehavior$Region == 'Whole Brain',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',], log1p_Trial_Power ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power*age, data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Occipital',], log1p_Trial_Power ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power*age, data = delayBehavior[delayBehavior$Region == 'Occipital',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Parietal',], log1p_Trial_Power ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power*age, data = delayBehavior[delayBehavior$Region == 'Parietal',])
  print(summary(model))
  
}

TrialPowerVarBehavAllRegions <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError, y = log1p_Trial_Power_Variability, group = Region, color = Region)) + stat_smooth(method='lm', se = T) + ggtitle("Trial Power Var vs Best Saccade")  +theme_classic()))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Whole Brain',], log1p_Trial_Power_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power_Variability*age, data = delayBehavior[delayBehavior$Region == 'Whole Brain',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',], log1p_Trial_Power_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power_Variability*age, data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Occipital',], log1p_Trial_Power_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power_Variability*age, data = delayBehavior[delayBehavior$Region == 'Occipital',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Parietal',], log1p_Trial_Power_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ log1p_Trial_Power_Variability*age, data = delayBehavior[delayBehavior$Region == 'Parietal',])
  print(summary(model))
  
  
}

EventNumberBehavAllRegions <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError, y = Event_Number, group = Region, color = Region))  + stat_smooth(method='lm', se = T) + ggtitle("Event Number vs Best Saccade") +theme_classic()))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Whole Brain',], Event_Number ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number*age, data = delayBehavior[delayBehavior$Region == 'Whole Brain',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',], Event_Number ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number*age, data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Occipital',], Event_Number ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number*age, data = delayBehavior[delayBehavior$Region == 'Occipital',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Parietal',], Event_Number ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number*age, data = delayBehavior[delayBehavior$Region == 'Parietal',])
  print(summary(model))
  
  
}

EventNumberVarBehavAllRegions <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError , y = Event_Number_Variability, group = Region, color = Region)) + stat_smooth(method='lm', se = T) + ggtitle("Event Number Var vs Best Saccade") +theme_classic()))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Whole Brain',], Event_Number_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number_Variability*age, data = delayBehavior[delayBehavior$Region == 'Whole Brain',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',], Event_Number_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number_Variability*age, data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Occipital',], Event_Number_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number_Variability*age, data = delayBehavior[delayBehavior$Region == 'Occipital',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Parietal',], Event_Number_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Number_Variability*age, data = delayBehavior[delayBehavior$Region == 'Parietal',])
  print(summary(model))
  
  
}

EventDurationBehavAllRegions <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError , y = Event_Duration , group = Region, color = Region)) + stat_smooth(method='lm', se = T) + ggtitle("Event Duration vs Best Saccade") +theme_classic()))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Whole Brain',], Event_Duration ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration*age, data = delayBehavior[delayBehavior$Region == 'Whole Brain',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',], Event_Duration ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration*age, data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Occipital',], Event_Duration ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration*age, data = delayBehavior[delayBehavior$Region == 'Occipital',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Parietal',], Event_Duration ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration*age, data = delayBehavior[delayBehavior$Region == 'Parietal',])
  print(summary(model))
  
  
}

EventDurationVarBehavAllRegions <- function(Behavior, delay_alldata) {
  
  delayBehavior <- merge(Behavior, delay_alldata, by = c('Subject', 'age'))
  
  print(lunaize(ggplot(data = delayBehavior[], aes(x = absBestError , y = Event_Duration_Variability, group = Region, color = Region)) + stat_smooth(method='lm', se = T) + ggtitle("Event Duration Var vs Best Saccade") +theme_classic()))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Whole Brain',], Event_Duration_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration_Variability*age, data = delayBehavior[delayBehavior$Region == 'Whole Brain',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',], Event_Duration_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration_Variability*age, data = delayBehavior[delayBehavior$Region == 'SelectiveFrontal',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Occipital',], Event_Duration_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration_Variability*age, data = delayBehavior[delayBehavior$Region == 'Occipital',])
  print(summary(model))
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = delayBehavior[delayBehavior$Region == 'Parietal',], Event_Duration_Variability ~ absBestError)
  print(anova(lm.model))  
  
  # age interactions  
  model <- lm(absBestError ~ Event_Duration_Variability*age, data = delayBehavior[delayBehavior$Region == 'Parietal',])
  print(summary(model))
  
  
}
