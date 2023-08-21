
TrialPower <- function (alldata) {
  TP <- list()
  
  TPgraph <- lunaize(ggplot(data = alldata[], aes(x = age, y = log1p_Trial_Power))+ geom_point() + stat_smooth(method = "gam") + ggtitle("Trial Power") + xlab("Age") + ylab("Trial Power")) +theme(plot.title = element_text(hjust = 0.5))
  
  lm.model <- lm(data = alldata[], scale(log1p_Trial_Power) ~ scale(inverseAge))
  #print(anova(lm.model))  
 # print(summary(lm.model))
  
  TPanova <- anova(lm.model)
  TP$pvalue <- TPanova$`Pr(>F)`[1]
  TP$TPgraph <- TPgraph
  TP$model <- lm.model

  #show as bar graph 
  groupBehaviorAge <- aggregate(.~ Group, alldata[2:12], mean)
  sderror<-std.error(groupBehaviorAge$log1p_Trial_Power)
  
  YMIN <- groupBehaviorAge$log1p_Trial_Power-sderror
  
  YMAX <- groupBehaviorAge$log1p_Trial_Power+sderror
  
  TP$bargraph <- lunaize(ggplot(groupBehaviorAge, aes(x=Group, y=log1p_Trial_Power, ymin= YMIN, ymax=YMAX , group = Group, fill = Group)) + geom_col(position = 'dodge') + geom_errorbar(position='dodge'))
  
  lm(data = groupBehaviorAge, log1p_Trial_Power ~ Group)
  
  #significant periods of growth
  #gam.model <- gam(log1p_Trial_Power ~ s(age), data = alldata)
  #gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  #TP$sigGrowth <- gam_growthrate_plot(alldata, gam.model, gam.growthrate, agevar = 'age', yvar = 'log1p_Trial_Power', draw_points = F)
  
  
  return(TP)
}

TrialPowerVar <- function (alldata) {
  TPV <- list()
  
  
  TPVgraph <- (lunaize(ggplot(data = alldata, aes(x = age, y = log1p_Trial_Power_Variability)) + geom_point() + stat_smooth(method = "gam") + ggtitle("Trial Power Variability") + xlab("Age") + ylab("Trial Power Var.")) +theme(plot.title = element_text(hjust = 0.5)))
  
 
  lm.model <- lm(data = alldata[], scale(log1p_Trial_Power_Variability) ~ scale(inverseAge))
  #print(anova(lm.model))  
  #print(summary(lm.model))
  
  TPVanova <- anova(lm.model)
  TPV$pvalue <- TPVanova$`Pr(>F)`[1]
  
  TPV$TPVgraph <- TPVgraph
  TPV$model <- lm.model
  
  groupBehaviorAge <- aggregate(.~ Group, alldata[2:12], mean)
  sderror<-std.error(groupBehaviorAge$log1p_Trial_Power_Variability)
  
  YMIN <- groupBehaviorAge$log1p_Trial_Power_Variability-sderror
  
  YMAX <- groupBehaviorAge$log1p_Trial_Power_Variability+sderror
  
  
  TPV$bargraph <- lunaize(ggplot(groupBehaviorAge, aes(x=Group, y=log1p_Trial_Power_Variability, ymin= YMIN, ymax=YMAX , group = Group, fill = Group)) + geom_col(position = 'dodge') + geom_errorbar(position='dodge'))
  
  #significant periods of growth
  # gam.model <- gam(log1p_Trial_Power_Variability ~ s(age), data = alldata)
  # gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  # gam_growthrate_plot(alldata, gam.model, gam.growthrate, agevar = 'age', yvar = 'log1p_Trial_Power_Variability', draw_points = F)
  # 
  return(TPV)
}

EventNumber <- function (alldata) {
  EN <- list()
  
  ENgraph <- (lunaize(ggplot(data = alldata, aes(x = age, y = Event_Number)) + geom_point() + stat_smooth(method = "gam") + ggtitle("Event Number") + xlab("Age") + ylab("Event Number")) +theme(plot.title = element_text(hjust = 0.5)))
  
  
  lm.model <- lm(data = alldata[], scale(Event_Number) ~ scale(inverseAge))
  #print(anova(lm.model))  
  #print(summary(lm.model))
  
  
  ENanova <- anova(lm.model)
  EN$pvalue <- ENanova$`Pr(>F)`[1]
  
  EN$ENgraph <- ENgraph  
  EN$model <- lm.model
  
  
  groupBehaviorAge <- aggregate(.~ Group, alldata[2:12], mean)
  sderror<-std.error(groupBehaviorAge$Event_Number)
  
  YMIN <- groupBehaviorAge$Event_Number-sderror
  
  YMAX <- groupBehaviorAge$Event_Number+sderror
  
 
  EN$bargraph <- lunaize(ggplot(groupBehaviorAge, aes(x=Group, y=Event_Number, ymin= YMIN, ymax=YMAX , group = Group, fill = Group)) + geom_col(position = 'dodge') + geom_errorbar(position='dodge'))
  
  #significant periods of growth
  #gam.model <- gam(Event_Number ~ s(age), data = alldata)
  #gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  #gam_growthrate_plot(alldata, gam.model, gam.growthrate, agevar = 'age', yvar = 'Event_Number', draw_points = F)
 
  return(EN)
}

EventNumberVar <- function (alldata) {
  ENV <- list()

  ENVgraph <- (lunaize(ggplot(data = alldata, aes(x = age, y = Event_Number_Variability)) + geom_point() + stat_smooth(method = "gam") + ggtitle("Event Number Variability") + xlab("Age") + ylab("Event Number Var.")) +theme(plot.title = element_text(hjust = 0.5)))
  
  lm.model <- lm(data = alldata[], scale(Event_Number_Variability) ~ scale(inverseAge))
  #print(anova(lm.model))  
  #print(summary(lm.model))  
  
  ENVanova <- anova(lm.model)
  ENV$pvalue <- ENVanova$`Pr(>F)`[1]
  
  ENV$ENVgraph <- ENVgraph
  ENV$model <- lm.model
  
  
  groupBehaviorAge <- aggregate(.~ Group, alldata[2:12], mean)
  sderror<-std.error(groupBehaviorAge$Event_Number_Variability)
  
  YMIN <- groupBehaviorAge$Event_Number_Variability-sderror
  
  YMAX <- groupBehaviorAge$Event_Number_Variability+sderror
  
 
  
  ENV$bargraph <- lunaize(ggplot(groupBehaviorAge, aes(x=Group, y=Event_Number_Variability, ymin= YMIN, ymax=YMAX , group = Group, fill = Group)) + geom_col(position = 'dodge') + geom_errorbar(position='dodge'))
  
  #significant periods of growth
 # gam.model <- gam(Event_Number_Variability ~ s(age), data = alldata)
  #gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  #gam_growthrate_plot(alldata, gam.model, gam.growthrate, agevar = 'age', yvar = 'Event_Number_Variability', draw_points = F)
  
  
return(ENV)
  
}

EventDuration <- function (alldata) {
  ED <- list()
  
  EDgraph <- (lunaize(ggplot(data = alldata, aes(x = age, y = Event_Duration))  + geom_point() +  stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Event Duration") + xlab("Age") + ylab("Event Duration (s)")) +theme(plot.title = element_text(hjust = 0.5)))
  
  
  lm.model <- lm(data = alldata[], scale(Event_Duration) ~ scale(inverseAge))
 # print(anova(lm.model))  
 #print(summary(lm.model))  
  
  EDanova <- anova(lm.model)
  ED$pvalue <- EDanova$`Pr(>F)`[1]
  
  ED$EDgraph <- EDgraph   
  ED$model <- lm.model
  
  
  
  groupBehaviorAge <- aggregate(.~ Group, alldata[2:12], mean)
  sderror<-std.error(groupBehaviorAge$Event_Duration)
  
  YMIN <- groupBehaviorAge$Event_Duration-sderror
  
  YMAX <- groupBehaviorAge$Event_Duration+sderror
  

  
  ED$bargraph <- lunaize(ggplot(groupBehaviorAge, aes(x=Group, y=Event_Duration, ymin= YMIN, ymax=YMAX , group = Group, fill = Group)) + geom_col(position = 'dodge') + geom_errorbar(position='dodge'))
  
  #significant periods of growth
  #gam.model <- gam(Event_Duration ~ s(age), data = alldata)
  #gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  #gam_growthrate_plot(alldata, gam.model, gam.growthrate, agevar = 'age', yvar = 'Event_Duration', draw_points = F)
  
 return(ED)
  
}

EventDurationVar <- function (alldata) {
  EDV <- list()
  
 
  EDVgraph <- (lunaize(ggplot(data = alldata, aes(x = age, y = Event_Duration_Variability)) + geom_point() + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Event Duration Variability") + xlab("Age") + ylab("Event Duration Var (s)")) +theme(plot.title = element_text(hjust = 0.5)))
  

  lm.model <- lm(data = alldata[], scale(Event_Duration_Variability) ~ scale(inverseAge))
  #print(anova(lm.model))  
  #print(summary(lm.model)) 
  
  
  EDVanova <- anova(lm.model)
  EDV$pvalue <- EDVanova$`Pr(>F)`[1]
  
  EDV$EDVgraph <- EDVgraph    
  EDV$model <- lm.model
  
  
  
  groupBehaviorAge <- aggregate(.~ Group, alldata[2:12], mean)
  sderror<-std.error(groupBehaviorAge$Event_Duration_Variability)
  
  YMIN <- groupBehaviorAge$Event_Duration_Variability-sderror
  
  YMAX <- groupBehaviorAge$Event_Duration_Variability+sderror
  

  
  EDV$bargraph <- lunaize(ggplot(groupBehaviorAge, aes(x=Group, y=Event_Duration_Variability, ymin= YMIN, ymax=YMAX , group = Group, fill = Group)) + geom_col(position = 'dodge') + geom_errorbar(position='dodge'))

  #significant periods of growth
 # gam.model <- gam(Event_Duration_Variability ~ s(age), data = alldata)
  #gam.growthrate <- gam_growthrate(gam.model, agevar = 'age')
  #gam_growthrate_plot(alldata, gam.model, gam.growthrate, agevar = 'age', yvar = 'Event_Duration_Variability', draw_points = F)

  return(EDV)
 
}

TrialPowerAllRegions <- function (alldata) {
  TP <- list()
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = log1p_Trial_Power, group = Region, color = Region)) + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Trial Power") + xlab("Age") + ylab("Trial Power (log)")) +theme(plot.title = element_text(hjust = 0.5)))
  
  TPgraph <- lunaize(ggplot(data = alldata[], aes(x = age, y = log1p_Trial_Power, group = Region, color = Region)) + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Trial Power") + xlab("Age") + ylab("Trial Power (log)")) +theme(plot.title = element_text(hjust = 0.5))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Whole Brain',], log1p_Trial_Power ~inverseAge)
  print(anova(lm.model))  
  
  TPanova <- anova(lm.model)
  TP$pvalue <- TPanova$`Pr(>F)`[1]
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Frontal',], log1p_Trial_Power ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Selective Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'SelectiveFrontal',], log1p_Trial_Power ~inverseAge)
  print(anova(lm.model))  
  
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = alldata[alldata$Region == 'Occipital',], log1p_Trial_Power ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = alldata[alldata$Region == 'Parietal',], log1p_Trial_Power ~inverseAge)
  print(anova(lm.model))  
  
  return(TP)
}

TrialPowerVarAllRegions <- function (alldata) {
  TPV <- list()
  
  print(lunaize(ggplot(data = alldata, aes(x = age, y = log1p_Trial_Power_Variability, group = Region, color = Region)) + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Trial Power Variability") + xlab("Age") + ylab("Trial Power Variability (log)")) +theme(plot.title = element_text(hjust = 0.5)))
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Whole Brain',], log1p_Trial_Power_Variability ~inverseAge)
  print(anova(lm.model))  
  

  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Frontal',], log1p_Trial_Power_Variability ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Selective Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'SelectiveFrontal',], log1p_Trial_Power_Variability ~inverseAge)
  print(anova(lm.model)) 
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = alldata[alldata$Region == 'Occipital',], log1p_Trial_Power_Variability ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = alldata[alldata$Region == 'Parietal',], log1p_Trial_Power_Variability ~inverseAge)
  print(anova(lm.model)) 
  return(TPV)
  
}

EventNumberAllRegions <- function (alldata) {
  
  print(lunaize(ggplot(data = alldata, aes(x = age, y = Event_Number, group = Region, color = Region)) + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Event Number") + xlab("Age") + ylab("Event Number")) +theme(plot.title = element_text(hjust = 0.5)))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Whole Brain',], Event_Number ~inverseAge)
  print(anova(lm.model))  
  

  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Frontal',], Event_Number ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Selective Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'SelectiveFrontal',], Event_Number ~inverseAge)
  print(anova(lm.model))  
  
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = alldata[alldata$Region == 'Occipital',], Event_Number ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = alldata[alldata$Region == 'Parietal',], Event_Number ~inverseAge)
  print(anova(lm.model)) 
  
}

EventNumberVarAllRegions <- function (alldata) {
  
  
  print(lunaize(ggplot(data = alldata, aes(x = age, y = Event_Number_Variability, group = Region, color = Region)) + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Event Number Variability") + xlab("Age") + ylab("Event Number Var")) +theme(plot.title = element_text(hjust = 0.5)))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Whole Brain',], Event_Number_Variability ~inverseAge)
  print(anova(lm.model))  
  
  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Frontal',], Event_Number_Variability ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Selective Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'SelectiveFrontal',], Event_Number_Variability ~inverseAge)
  print(anova(lm.model))  
  
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = alldata[alldata$Region == 'Occipital',], Event_Number_Variability ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = alldata[alldata$Region == 'Parietal',], Event_Number_Variability ~inverseAge)
  print(anova(lm.model)) 
}

EventDurationAllRegions <- function (alldata) {
  
  
  print(lunaize(ggplot(data = alldata, aes(x = age, y = Event_Duration, group = Region, color = Region)) + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Event Duration") + xlab("Age") + ylab("Event Duration (s)")) +theme(plot.title = element_text(hjust = 0.5)))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Whole Brain',], Event_Duration ~inverseAge)
  print(anova(lm.model))  

  
  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Frontal',], Event_Duration ~inverseAge)
  print(anova(lm.model))  
  
  
  cat(paste("",
            "",
            "Selective Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'SelectiveFrontal',], Event_Duration ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = alldata[alldata$Region == 'Occipital',], Event_Duration ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = alldata[alldata$Region == 'Parietal',], Event_Duration ~inverseAge)
  print(anova(lm.model)) 
}

EventDurationVarAllRegions <- function (alldata) {
  
  
  print(lunaize(ggplot(data = alldata, aes(x = age, y = Event_Duration_Variability, group = Region, color = Region)) + stat_smooth(method = "lm",formula='y~I(1/x)') + ggtitle("Event Duration Variability") + xlab("Age") + ylab("Event Duration Var (s)")) +theme(plot.title = element_text(hjust = 0.5)))
  
  cat(paste("",
            "",
            "Whole Brain Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Whole Brain',], Event_Duration_Variability ~inverseAge)
  print(anova(lm.model))  

  cat(paste("",
            "",
            "Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'Frontal',], Event_Duration_Variability ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Selective Frontal Statistics",
            "", sep = "\n"))
  
  lm.model <- lm(data = alldata[alldata$Region == 'SelectiveFrontal',], Event_Duration_Variability ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Occipital Statistics",
            "", sep = "\n"))  
  
  lm.model <- lm(data = alldata[alldata$Region == 'Occipital',], Event_Duration_Variability ~inverseAge)
  print(anova(lm.model))  
  
  cat(paste("",
            "",
            "Parietal Statistics",
            "", sep = "\n")) 
  
  lm.model <- lm(data = alldata[alldata$Region == 'Parietal',], Event_Duration_Variability ~inverseAge)
  print(anova(lm.model)) 
  
}