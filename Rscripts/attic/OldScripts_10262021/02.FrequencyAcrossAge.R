


## Delay Period
### Average Trial Power- Trial Level 
AverageTrialPower <- function(delayOnly_SubLevel) {
  
  
  #gamma 
  gammaTP <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Gamma.log1p_Trial_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + xlab("Age") + ylab("log(Power)") + ggtitle("Average Trial Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "A"))
              
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Gamma.log1p_Trial_Power)) + stat_smooth(method = 'gam') +geom_point() + ggtitle("Delay Period: Gamma Power across Age")  + xlab("Age") + ylab("log(Power)")) + theme(plot.title = element_text(hjust = 0.5)))
  
  lm.model <- lmer(data = delayOnly_SubLevel[], Gamma.log1p_Trial_Power ~ inverseAge)
  print(car::Anova(lm.model))
  
  #beta
  betaTP <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Beta.log1p_Trial_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Beta.log1p_Trial_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Beta Power across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("log(Power)"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Beta.log1p_Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaTP <-  (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Alpha.log1p_Trial_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + xlab("Age") + ylab("log(Power)") + ggtitle("Average Trial Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "A"))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Alpha.log1p_Trial_Power))   + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Alpha Power across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("log(Power)"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Alpha.log1p_Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaTP <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Theta.log1p_Trial_Power))  + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Theta.log1p_Trial_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Theta Power across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("log(Power)"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Theta.log1p_Trial_Power ~ inverseAge)
  print(anova(lm.model))
}

### Trial Power Variability

TrialPowerVariability <- function (alldata) {
  
  
  #gamma 
  gammaTPV <-  (ggplot(data = alldata[], aes(x = age, y = Gamma.log1p_Trial_Power_Variability))  + stat_smooth(method = 'lm',  formula='y~I(1/x)')   + xlab("Age") + ylab("log(Power) (sd)") + ggtitle("Trial Variability: Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "B"))
  
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Gamma.log1p_Trial_Power_Variability))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Gamma Power Variabilitly across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("log(Power)"))
  
  lm.model <- lm(data = alldata[], Gamma.log1p_Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  betaTPV <- (ggplot(data = alldata[], aes(x = age, y = Beta.log1p_Trial_Power_Variability))+ stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Beta.log1p_Trial_Power_Variability))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Beta Power Variability across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("log(Power)"))
  
  lm.model <- lm(data = alldata[], Beta.log1p_Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaTPV <- (ggplot(data = alldata[], aes(x = age, y = Alpha.Trial_Power_Variability))  + stat_smooth(method = 'lm',  formula='y~I(1/x)')   + xlab("Age") + ylab("log(Power) (sd)") + ggtitle("Trial Variability: Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "B"))
  
  
 print(lunaize(ggplot(data = alldata[], aes(x = age, y = Alpha.Trial_Power_Variability))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Alpha Power Variability across Age") )+ theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("log(Power)"))
  
  
  lm.model <- lm(data = alldata[], Alpha.Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaTPV <- (ggplot(data = alldata[], aes(x = age, y = Theta.log1p_Trial_Power_Variability))  + stat_smooth(method = 'lm',  formula='y~I(1/x)')   + xlab("Age") + ylab("log(Power) (sd)") + ggtitle("Trial Variability: Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "B"))
  
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Theta.log1p_Trial_Power_Variability))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Theta Power Varbility across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("log(Power)"))
  
  lm.model <- lm(data = alldata[], Theta.log1p_Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
}

### Number of Events

NumberofEvents <- function(delayOnly_SubLevel) {
  
  #gamma 
 gammaNE <- (ggplot(data = delayOnly_SubLevel, aes(x = age, y = Gamma.Event_Number))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + xlab("Age") + ylab("Number of Bursts") + ggtitle("Average Number of Gamma Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "C"))
 
 
 print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Gamma.Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Gamma Events across Age")  + xlab("Age") + ylab("Number of Events")))+ theme(plot.title = element_text(hjust = 0.5))
  
  lm.model <- lmer(data = delayOnly_SubLevel[], Gamma.Event_Number ~ inverseAge)
  print(car::Anova(lm.model))
  
  #beta
  betaNE <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Beta.Event_Number))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Beta.Event_Number))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Beta Events across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Beta.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaNE <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Alpha.Event_Number))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + xlab("Age") + ylab("Number of Bursts") + ggtitle("Average Number of Alpha Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "C"))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Alpha.Event_Number))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Alpha Events across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Alpha.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaNE <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Theta.Event_Number))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Theta.Event_Number))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Theta Events across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Theta.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
}

### Number of Events Variabiltiy

NumberofEventsVariability <- function(alldata) {
  
  
  #gamma 
 gammaNEV <-  (ggplot(data = alldata[], aes(x = age, y = Gamma.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5))  + xlab("Age") + ylab("Number of Bursts (sd)") + ggtitle("Trial Variability: Number of Gamma Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "D"))
 
 print(lunaize(ggplot(data = alldata[], aes(x = age, y = Gamma.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Gamma Events (Variability) across Age") )+ theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata[], Gamma.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #beta
 betaNEV <- (ggplot(data = alldata[], aes(x = age, y = Beta.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
 
 print(lunaize(ggplot(data = alldata[], aes(x = age, y = Beta.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Beta Events (Variability) across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata[], Beta.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaNEV <- (ggplot(data = alldata[], aes(x = age, y = Alpha.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5))  + xlab("Age") + ylab("Number of Bursts (sd)") + ggtitle("Trial Variability: Number of Alpha Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "D"))
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Alpha.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Alpha Events (Variability) across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata[], Alpha.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
 thetaNEV <- (ggplot(data = alldata[], aes(x = age, y = Theta.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
 
 print(lunaize(ggplot(data = alldata[], aes(x = age, y = Theta.Event_Number_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Number of Theta Events (Variability) across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata[], Theta.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
}

### Duration of Events

DurationofEvents <- function(delayOnly_SubLevel) {
  
  #gamma 
  gammaDE <- (ggplot(data = delayOnly_SubLevel, aes(x = age, y = Gamma.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Gamma.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Duration of \nGamma Events across Age")  + xlab("Age") + ylab("Duration of Event"))+ theme(plot.title = element_text(hjust = 0.5)))
  
  lm.model <- lm(data = delayOnly_SubLevel, Gamma.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  #beta
 betaDE <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Beta.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
 
 print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Beta.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Average Duration of Beta Events across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Beta.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaDE <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Alpha.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')+ theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Alpha.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Average Duration of Alpha Events across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event"))
  
  lm.model <- lm(data = delayOnly_SubLevel[], Alpha.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
 thetaDE <- (ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Theta.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
 print(lunaize(ggplot(data = delayOnly_SubLevel[], aes(x = age, y = Theta.Event_Duration))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Average Duration of Theta Events across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration"))

 
  lm.model <- lm(data = delayOnly_SubLevel[], Theta.Event_Duration ~ inverseAge)
  print(anova(lm.model))
}

### Duration of Events Variability
DurationofEventsVariability <- function(alldata) {
  
  
  #gamma 
  gammaDEV <- (ggplot(data = alldata[], aes(x = age, y = Gamma.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Gamma.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Duration of  Gamma Events (Variability) across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event (sd)"))
  
  lm.model <- lm(data = alldata[], Gamma.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  betaDEV <- (ggplot(data = alldata[], aes(x = age, y = Beta.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Beta.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Duration of Beta Events (Variability) across Age") )+ theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event (sd)"))
  
  lm.model <- lm(data = alldata[], Beta.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaDEV <- (ggplot(data = alldata[], aes(x = age, y = Alpha.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Alpha.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Duration of Alpha Events (Variability) across Age") )+ theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event (sd)"))
  
  lm.model <- lm(data = alldata[], Alpha.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaDEV <- (ggplot(data = alldata[], aes(x = age, y = Theta.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(lunaize(ggplot(data = alldata[], aes(x = age, y = Theta.Event_Duration_Variability))  +  stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Duration of Theta Events (Variability) across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event (sd)"))
  
  lm.model <- lm(data = alldata[], Theta.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
}

CreatePaperFiguresDelay <- function() {
  
  cowplot::plot_grid(lunaize(gammaTP),lunaize(gammaNE),lunaize(gammaDE), lunaize(betaTP), lunaize(betaNE), lunaize(betaDE), lunaize(alphaTP), lunaize(alphaNE), lunaize(alphaDE), lunaize(thetaTP),lunaize(thetaNE), lunaize(thetaDE), nrow = 4, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G','H','I','J','K','L'), label_size = 20) 
  
#variability
  cowplot::plot_grid(lunaize(gammaTPV),lunaize(gammaNEV),lunaize(gammaDEV), lunaize(betaTPV), lunaize(betaNEV), lunaize(betaDEV), lunaize(alphaTPV), lunaize(alphaNEV), lunaize(alphaDEV), lunaize(thetaTPV),lunaize(thetaNEV), lunaize(thetaDEV), nrow = 4, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G','H','I','J','K','L'), label_size = 20) 
  
  #beta only 
  cowplot::plot_grid(lunaize(betaTP), lunaize(betaNE), lunaize(betaDE), lunaize(betaTPV), lunaize(betaNEV), lunaize(betaDEV), nrow = 2, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G'), label_size = 20) 
  
  #gamma only
  cowplot::plot_grid(lunaize(gammaTP), lunaize(gammaNE), lunaize(gammaDE), lunaize(gammaTPV), lunaize(gammaNEV), lunaize(gammaDEV), nrow = 2, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G'), label_size = 20) 
  
  library(gridExtra)
  grid.arrange(gammaTP, gammaTPV, gammaNE, gammaNEV, nrow = 2, ncol = 2, widths = c(1,1), heights = c(2,2), layout_matrix = rbind(c(1,2),c(3,4)))
  
  
  #alpha only
  cowplot::plot_grid(lunaize(alphaTP),lunaize(alphaTPV), lunaize(alphaNE), lunaize(alphaNEV),  nrow = 2, ncol=2, labels= c('A', 'B', 'C', 'D'), label_size = 20) 
  
  grid.arrange(alphaTP, alphaTPV, alphaNE, alphaNEV, nrow = 2, ncol = 2, widths = c(1,1), heights = c(2,2), layout_matrix = rbind(c(1,2),c(3,4)))
  
 
  
 
  
  
}  
  

## Fixation Period
### Average Trial Power 
AverageTrialPowerFix <- function() {
  
  avgTheta <- aggregate(Theta.log1p_Trial_Power_Fix ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.log1p_Trial_Power_Fix ~ visitno , alldata, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.log1p_Trial_Power_Fix ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.log1p_Trial_Power_Fix ~ visitno , alldata, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.log1p_Trial_Power_Fix ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.log1p_Trial_Power_Fix ~ visitno , alldata, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.log1p_Trial_Power_Fix ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.log1p_Trial_Power_Fix ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.log1p_Trial_Power_Fix < gammaCutoff,], aes(x = age, y = Gamma.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Gamma Power across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 gammaTP <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.log1p_Trial_Power_Fix < gammaCutoff,], aes(x = age, y = Gamma.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.log1p_Trial_Power_Fix < gammaCutoff,], Gamma.log1p_Trial_Power_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.log1p_Trial_Power_Fix < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Beta Power across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))

  betaTP <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.log1p_Trial_Power_Fix < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
    
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.log1p_Trial_Power_Fix < betaCutoff,], Beta.log1p_Trial_Power_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.log1p_Trial_Power_Fix < alphaCutoff,], aes(x = age, y = Alpha.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Alpha Power across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  alphaTP <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.log1p_Trial_Power_Fix < alphaCutoff,], aes(x = age, y = Alpha.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.log1p_Trial_Power_Fix < alphaCutoff,], Alpha.log1p_Trial_Power_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.log1p_Trial_Power_Fix < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Theta Power across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 thetaTP <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.log1p_Trial_Power_Fix < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.log1p_Trial_Power_Fix < thetaCutoff,], Theta.log1p_Trial_Power_Fix ~ inverseAge)
  print(anova(lm.model))
}

### Trial Power Variability

TrialPowerVariabilityFix <- function() {
  
  avgTheta <- aggregate(Theta.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.log1p_Trial_Power_Variability_Fix ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.log1p_Trial_Power_Variability_Fix < gammaCutoff,], aes(x = age, y = Gamma.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Gamma Power Variabilitly across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  gammaTPV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.log1p_Trial_Power_Variability_Fix < gammaCutoff,], aes(x = age, y = Gamma.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.log1p_Trial_Power_Variability_Fix < gammaCutoff,], Gamma.log1p_Trial_Power_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.log1p_Trial_Power_Variability_Fix < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Beta Power Variability across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  betaTPV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.log1p_Trial_Power_Variability_Fix < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.log1p_Trial_Power_Variability_Fix < betaCutoff,], Beta.log1p_Trial_Power_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.log1p_Trial_Power_Variability_Fix < alphaCutoff,], aes(x = age, y = Alpha.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Alpha Power Variability across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  alphaTPV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.log1p_Trial_Power_Variability_Fix < alphaCutoff,], aes(x = age, y = Alpha.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.log1p_Trial_Power_Variability_Fix < alphaCutoff,], Alpha.log1p_Trial_Power_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.log1p_Trial_Power_Variability_Fix < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Trial Theta Power Varbility across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 thetaTPV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.log1p_Trial_Power_Variability_Fix < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.log1p_Trial_Power_Variability_Fix < thetaCutoff,], Theta.log1p_Trial_Power_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
}

### Number of Events

NumberofEventsFix <- function() {
  
  avgTheta <- aggregate(Theta.Event_Number_Fix ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Number_Fix ~ visitno , alldata, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.Event_Number_Fix ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Number_Fix ~ visitno , alldata, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.Event_Number_Fix ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Number_Fix ~ visitno , alldata, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.Event_Number_Fix ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Number_Fix ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  gammaCutofflower <-  avgGamma[1,2] - 2.5* sdGamma[1,2]
  
  
  #gamma 
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Fix < gammaCutoff  & alldata$Gamma.Event_Number_Fix > gammaCutofflower,], aes(x = age, y = Gamma.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Number of  Gamma Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 gammaNE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Fix < gammaCutoff & alldata$Gamma.Event_Number_Fix > gammaCutofflower,], aes(x = age, y = Gamma.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Fix < gammaCutoff & alldata$Gamma.Event_Number_Fix > gammaCutofflower,], Gamma.Event_Number_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Number of Beta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  betaNE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Fix < betaCutoff,], Beta.Event_Number_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Number of Alpha Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  alphaNE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Fix < alphaCutoff,], Alpha.Event_Number_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Number of Theta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  thetaNE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Number_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Fix < thetaCutoff,], Theta.Event_Number_Fix ~ inverseAge)
  print(anova(lm.model))
  
}

### Number of Events Variabiltiy

NumberofEventsVariabilityFix <- function() {
  
  avgTheta <- aggregate(Theta.Event_Number_Variability_Fix ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Number_Variability_Fix ~ visitno , alldata, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.Event_Number_Variability_Fix ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Number_Variability_Fix ~ visitno , alldata, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.Event_Number_Variability_Fix ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Number_Variability_Fix ~ visitno , alldata, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.Event_Number_Variability_Fix ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Number_Variability_Fix ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Variability_Fix < gammaCutoff,], aes(x = age, y = Gamma.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Number of Gamma Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 gammaNEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Variability_Fix < gammaCutoff,], aes(x = age, y = Gamma.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Variability_Fix < gammaCutoff,], Gamma.Event_Number_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Variability_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Number of Beta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  betaNEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Variability_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Variability_Fix < betaCutoff,], Beta.Event_Number_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Variability_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Number of Alpha Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 alphaNEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Variability_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Variability_Fix < alphaCutoff,], Alpha.Event_Number_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Variability_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Number of Theta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  thetaNEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Variability_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Number_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Variability_Fix < thetaCutoff,], Theta.Event_Number_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
}

### Duration of Events

DurationofEventsFix <- function() {
  
  avgTheta <- aggregate(Theta.Event_Duration_Fix ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Duration_Fix ~ visitno , alldata, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.Event_Duration_Fix ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration_Fix ~ visitno , alldata, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.Event_Duration_Fix ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Duration_Fix ~ visitno , alldata, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.Event_Duration_Fix ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Duration_Fix ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Fix < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Duration of  Gamma Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  gammaDE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Fix < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Fix < gammaCutoff,], Gamma.Event_Duration_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Duration of Beta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  betaDE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Fix < betaCutoff,], Beta.Event_Duration_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Duration of Alpha Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 alphaDE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Fix < alphaCutoff,], Alpha.Event_Duration_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Average Duration of Theta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  thetaDE <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Duration_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Fix < thetaCutoff,], Theta.Event_Duration_Fix ~ inverseAge)
  print(anova(lm.model))
}

### Duration of Events Variability

DurationofEventsVariabilityFix <- function() {
  
  avgTheta <- aggregate(Theta.Event_Duration_Variability_Fix ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Duration_Variability_Fix ~ visitno , alldata, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.Event_Duration_Variability_Fix ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration_Variability_Fix ~ visitno , alldata, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.Event_Duration_Variability_Fix ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Duration_Variability_Fix ~ visitno , alldata, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.Event_Duration_Variability_Fix ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Duration_Variability_Fix ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Variability_Fix < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Duration of  Gamma Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
 gammaDEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Variability_Fix < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Variability_Fix < gammaCutoff,], Gamma.Event_Duration_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Variability_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Duration of Beta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  betaDEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Variability_Fix < betaCutoff,], aes(x = age, y = Beta.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Variability_Fix < betaCutoff,], Beta.Event_Duration_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Variability_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Duration of Alpha Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  alphaDEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Variability_Fix < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Variability_Fix < alphaCutoff,], Alpha.Event_Duration_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Variability_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Fix Period: Duration of Theta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  thetaDEV <- (ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Variability_Fix < thetaCutoff,], aes(x = age, y = Theta.Event_Duration_Variability_Fix))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Variability_Fix < thetaCutoff,], Theta.Event_Duration_Variability_Fix ~ inverseAge)
  print(anova(lm.model))
  
}

CreatePaperFiguresFixation <- function() {

cowplot::plot_grid(lunaize(gammaTP),lunaize(gammaNE),lunaize(gammaDE), lunaize(betaTP), lunaize(betaNE), lunaize(betaDE), lunaize(alphaTP), lunaize(alphaNE), lunaize(alphaDE), lunaize(thetaTP),lunaize(thetaNE), lunaize(thetaDE), nrow = 4, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G','H','I','J','K','L'), label_size = 20) 

#variability

cowplot::plot_grid(lunaize(gammaTPV),lunaize(gammaNEV),lunaize(gammaDEV), lunaize(betaTPV), lunaize(betaNEV), lunaize(betaDEV), lunaize(alphaTPV), lunaize(alphaNEV), lunaize(alphaDEV), lunaize(thetaTPV),lunaize(thetaNEV), lunaize(thetaDEV), nrow = 4, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G','H','I','J','K','L'), label_size = 20) 

#beta only 
cowplot::plot_grid(lunaize(betaTP), lunaize(betaNE), lunaize(betaDE), lunaize(betaTPV), lunaize(betaNEV), lunaize(betaDEV), nrow = 2, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G'), label_size = 20) 

#gamma only
cowplot::plot_grid(lunaize(gammaTP), lunaize(gammaNE), lunaize(gammaDE), lunaize(gammaTPV), lunaize(gammaNEV), lunaize(gammaDEV), nrow = 2, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G'), label_size = 20) 

#alpha only
cowplot::plot_grid(lunaize(alphaTP), lunaize(alphaNE), lunaize(alphaDE), lunaize(alphaTPV), lunaize(alphaNEV), lunaize(alphaDEV), nrow = 2, ncol=3, labels= c('A', 'B', 'C', 'D','E','F','G'), label_size = 20) 
}


ZscoreStats <- function () {
  
  test_df <- subset(AllGamma[], select = c(Subject,Age, Gamma_Event_Duration_sd)) 
  test_df$inverseAge <- 1/test_df$Age
  zscore <- as.data.frame(scale(test_df[,3-4]))
  
  print(ggplot(data = zscore[], aes(x = age, y = Gamma.Trial_Power_Variability_Delay))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay Period: Average Trial Gamma Power Variabilitly across Age (z-scored)") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = zscore[], Gamma_Event_Duration_sd ~ inverseAge)
  print(anova(lm.model))
  summary(lm.model)
  
}

DelayMinusFix <- function() {
  
  alldata$Alpha.Trial_Power <- (alldata$Alpha.Trial_Power_Delay - alldata$Alpha.Trial_Power_Fix)
  alldata$Theta.Trial_Power <- (alldata$Theta.Trial_Power_Delay - alldata$Theta.Trial_Power_Fix)
  alldata$Beta.Trial_Power <- (alldata$Beta.Trial_Power_Delay - alldata$Beta.Trial_Power_Fix)
  alldata$Gamma.Trial_Power <- (alldata$Gamma.Trial_Power_Delay - alldata$Gamma.Trial_Power_Fix)
  
  
  alldata$Alpha.Trial_Power_Variability <- (alldata$Alpha.Trial_Power_Variability_Delay - alldata$Alpha.Trial_Power_Variability_Fix)
  alldata$Theta.Trial_Power_Variability <- (alldata$Theta.Trial_Power_Variability_Delay - alldata$Theta.Trial_Power_Variability_Fix)
  alldata$Beta.Trial_Power_Variability <- (alldata$Beta.Trial_Power_Variability_Delay - alldata$Beta.Trial_Power_Variability_Fix)
  alldata$Gamma.Trial_Power_Variability <- (alldata$Gamma.Trial_Power_Variability_Delay - alldata$Gamma.Trial_Power_Variability_Fix)
  
  alldata$Alpha.Event_Number <- (alldata$Alpha.Event_Number_Delay - alldata$Alpha.Event_Number_Fix)
  alldata$Theta.Event_Number <- (alldata$Theta.Event_Number_Delay - alldata$Theta.Event_Number_Fix)
  alldata$Beta.Event_Number <- (alldata$Beta.Event_Number_Delay - alldata$Beta.Event_Number_Fix)
  alldata$Gamma.Event_Number <- (alldata$Gamma.Event_Number_Delay - alldata$Gamma.Event_Number_Fix)
  
  alldata$Alpha.Event_Number_Variability <- (alldata$Alpha.Event_Number_Variability_Delay - alldata$Alpha.Event_Number_Variability_Fix)
  alldata$Theta.Event_Number_Variability <- (alldata$Theta.Event_Number_Variability_Delay - alldata$Theta.Event_Number_Variability_Fix)
  alldata$Beta.Event_Number_Variability <- (alldata$Beta.Event_Number_Variability_Delay - alldata$Beta.Event_Number_Variability_Fix)
  alldata$Gamma.Event_Number_Variability <- (alldata$Gamma.Event_Number_Variability_Delay - alldata$Gamma.Event_Number_Variability_Fix)
  
  alldata$Alpha.Event_Duration <- (alldata$Alpha.Event_Duration_Delay - alldata$Alpha.Event_Duration_Fix)
  alldata$Theta.Event_Duration <- (alldata$Theta.Event_Duration_Delay - alldata$Theta.Event_Duration_Fix)
  alldata$Beta.Event_Duration <- (alldata$Beta.Event_Duration_Delay - alldata$Beta.Event_Duration_Fix)
  alldata$Gamma.Event_Duration <- (alldata$Gamma.Event_Duration_Delay - alldata$Gamma.Event_Duration_Fix)
  
  alldata$Alpha.Event_Duration_Variability <- (alldata$Alpha.Event_Duration_Variability_Delay - alldata$Alpha.Event_Duration_Variability_Fix)
  alldata$Theta.Event_Duration_Variability <- (alldata$Theta.Event_Duration_Variability_Delay - alldata$Theta.Event_Duration_Variability_Fix)
  alldata$Beta.Event_Duration_Variability <- (alldata$Beta.Event_Duration_Variability_Delay - alldata$Beta.Event_Duration_Variability_Fix)
  alldata$Gamma.Event_Duration_Variability <- (alldata$Gamma.Event_Duration_Variability_Delay - alldata$Gamma.Event_Duration_Variability_Fix)
  
  return(alldata)
}

DelayMinusFix_TrialPower <- function () {
  
  avgGamma <- aggregate(Gamma.Trial_Power ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Trial_Power ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Trial_Power ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Trial_Power ~ visitno , alldata, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Trial_Power ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Trial_Power ~ visitno , alldata, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Trial_Power ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Trial_Power ~ visitno , alldata, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Trial_Power < gammaCutoff & alldata$Gamma.Trial_Power > gammaCutofflower,], aes(x = age, y = Gamma.Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Trial_Power < gammaCutoff & alldata$Gamma.Trial_Power > gammaCutofflower,], Gamma.Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Trial_Power < BetaCutoff & alldata$Beta.Trial_Power > BetaCutofflower,], aes(x = age, y = Beta.Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Trial_Power < BetaCutoff & alldata$Beta.Trial_Power > BetaCutofflower,], Beta.Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Trial_Power < AlphaCutoff & alldata$Alpha.Trial_Power > AlphaCutofflower,], aes(x = age, y = Alpha.Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Trial_Power < AlphaCutoff & alldata$Alpha.Trial_Power > AlphaCutofflower,], Alpha.Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Trial_Power < ThetaCutoff & alldata$Theta.Trial_Power > ThetaCutofflower,], aes(x = age, y = Theta.Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Trial_Power < ThetaCutoff & alldata$Theta.Trial_Power > ThetaCutofflower,], Theta.Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  
}


DelayMinusFix_TrialPowerVar <- function () {
  
 
  avgGamma <- aggregate(Gamma.Trial_Power_Variability ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Trial_Power_Variability ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Trial_Power_Variability ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Trial_Power_Variability ~ visitno , alldata, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Trial_Power_Variability ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Trial_Power_Variability ~ visitno , alldata, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Trial_Power_Variability ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Trial_Power_Variability ~ visitno , alldata, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Trial_Power_Variability < gammaCutoff & alldata$Gamma.Trial_Power_Variability > gammaCutofflower,], aes(x = age, y = Gamma.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Trial_Power_Variability < gammaCutoff & alldata$Gamma.Trial_Power > gammaCutofflower,], Gamma.Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Trial_Power_Variability < BetaCutoff & alldata$Beta.Trial_Power_Variability > BetaCutofflower,], aes(x = age, y = Beta.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Trial_Power_Variability < BetaCutoff & alldata$Beta.Trial_Power_Variability > BetaCutofflower,], Beta.Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Trial_Power_Variability < AlphaCutoff & alldata$Alpha.Trial_Power_Variability > AlphaCutofflower,], aes(x = age, y = Alpha.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Trial_Power_Variability < AlphaCutoff & alldata$Alpha.Trial_Power_Variability > AlphaCutofflower,], Alpha.Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Trial_Power_Variability < ThetaCutoff & alldata$Theta.Trial_Power_Variability > ThetaCutofflower,], aes(x = age, y = Theta.Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Trial_Power_Variability < ThetaCutoff & alldata$Theta.Trial_Power_Variability > ThetaCutofflower,], Theta.Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  
}


DelayMinusFix_EventNumber <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Number ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Number ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Number ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Number ~ visitno , alldata, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Number ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Number ~ visitno , alldata, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Number ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Number ~ visitno , alldata, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number < gammaCutoff & alldata$Gamma.Event_Number > gammaCutofflower,], aes(x = age, y = Gamma.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number < gammaCutoff & alldata$Gamma.Event_Number > gammaCutofflower,], Gamma.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number < BetaCutoff & alldata$Beta.Event_Number > BetaCutofflower,], aes(x = age, y = Beta.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number < BetaCutoff & alldata$Beta.Event_Number > BetaCutofflower,], Beta.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number < AlphaCutoff & alldata$Alpha.Event_Number > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number < AlphaCutoff & alldata$Alpha.Event_Number > AlphaCutofflower,], Alpha.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number < ThetaCutoff & alldata$Theta.Event_Number > ThetaCutofflower,], aes(x = age, y = Theta.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number < ThetaCutoff & alldata$Theta.Event_Number > ThetaCutofflower,], Theta.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  
  
  
}


DelayMinusFix_EventNumberVar <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Number_Variability ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Number_Variability ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Number_Variability ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Number_Variability ~ visitno , alldata, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Number_Variability ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Number_Variability ~ visitno , alldata, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Number_Variability ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Number_Variability ~ visitno , alldata, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Variability < gammaCutoff & alldata$Gamma.Event_Number_Variability > gammaCutofflower,], aes(x = age, y = Gamma.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Number_Variability < gammaCutoff & alldata$Gamma.Event_Number_Variability > gammaCutofflower,], Gamma.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Variability < BetaCutoff & alldata$Beta.Event_Number_Variability > BetaCutofflower,], aes(x = age, y = Beta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Number_Variability < BetaCutoff & alldata$Beta.Event_Number_Variability > BetaCutofflower,], Beta.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Variability < AlphaCutoff & alldata$Alpha.Event_Number_Variability > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Number_Variability < AlphaCutoff & alldata$Alpha.Event_Number_Variability > AlphaCutofflower,], Alpha.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Variability < ThetaCutoff & alldata$Theta.Event_Number_Variability > ThetaCutofflower,], aes(x = age, y = Theta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Number_Variability < ThetaCutoff & alldata$Theta.Event_Number_Variability > ThetaCutofflower,], Theta.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  
  
}


DelayMinusFix_EventDuration <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Duration ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Duration ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Duration ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Duration ~ visitno , alldata, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Duration ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration ~ visitno , alldata, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Duration ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Duration ~ visitno , alldata, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration < gammaCutoff & alldata$Gamma.Event_Duration > gammaCutofflower,], aes(x = age, y = Gamma.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration < gammaCutoff & alldata$Gamma.Event_Duration > gammaCutofflower,], Gamma.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration < BetaCutoff & alldata$Beta.Event_Duration > BetaCutofflower,], aes(x = age, y = Beta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration < BetaCutoff & alldata$Beta.Event_Duration > BetaCutofflower,], Beta.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration < AlphaCutoff & alldata$Alpha.Event_Duration > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration < AlphaCutoff & alldata$Alpha.Event_Duration > AlphaCutofflower,], Alpha.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration < ThetaCutoff & alldata$Theta.Event_Duration > ThetaCutofflower,], aes(x = age, y = Theta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration < ThetaCutoff & alldata$Theta.Event_Duration > ThetaCutofflower,], Theta.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  
  
}


DelayMinusFix_EventDurationVar <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Duration_Variability ~ visitno , alldata, mean)
  sdGamma <- aggregate(Gamma.Event_Duration_Variability ~ visitno , alldata, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Duration_Variability ~ visitno , alldata, mean)
  sdBeta <- aggregate(Beta.Event_Duration_Variability ~ visitno , alldata, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Duration_Variability ~ visitno , alldata, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration_Variability ~ visitno , alldata, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Duration_Variability ~ visitno , alldata, mean)
  sdTheta <- aggregate(Theta.Event_Duration_Variability ~ visitno , alldata, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Variability < gammaCutoff & alldata$Gamma.Event_Duration_Variability > gammaCutofflower,], aes(x = age, y = Gamma.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Gamma.Event_Duration_Variability < gammaCutoff & alldata$Gamma.Event_Duration_Variability > gammaCutofflower,], Gamma.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Variability < BetaCutoff & alldata$Beta.Event_Duration_Variability > BetaCutofflower,], aes(x = age, y = Beta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Beta.Event_Duration_Variability < BetaCutoff & alldata$Beta.Event_Duration_Variability > BetaCutofflower,], Beta.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Variability < AlphaCutoff & alldata$Alpha.Event_Duration_Variability > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Alpha.Event_Duration_Variability < AlphaCutoff & alldata$Alpha.Event_Duration_Variability > AlphaCutofflower,], Alpha.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Variability < ThetaCutoff & alldata$Theta.Event_Duration_Variability > ThetaCutofflower,], aes(x = age, y = Theta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = alldata[alldata$visitno < 2 & alldata$Theta.Event_Duration_Variability < ThetaCutoff & alldata$Theta.Event_Duration_Variability > ThetaCutofflower,], Theta.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  
  
}


## Resting State 
### Average Trial Power - trial level 
AverageTrialPowerRS <- function() {
  
  alldf_RS <- Resting_State_Data_TrialLevel()
  alldata_RS_trialLevel <- alldf_RS$alldata_RS_TrialLevel
  
  Gamma_cleandata <-  alldata_RS_trialLevel %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2))
  
  
  #gamma 
  gammaTP <- (ggplot(data = Gamma_cleandata[], aes(x = age, y = Gamma.log1p_Trial_Power))  + stat_smooth(method = 'lm',  formula='y~I(1/x)') + xlab("Age") + ylab("log(Power)") + ggtitle("Average Trial Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "A"))
  
  
  lunaize(print(ggplot(data = Gamma_cleandata[], aes(x = age, y = Gamma.log1p_Trial_Power)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Trial Gamma Power across Age")  + xlab("Age") + ylab("log(Power)")))+ theme(plot.title = element_text(hjust = 0.5))
  
  lm.model <- lmer(data = Gamma_cleandata[], Gamma.log1p_Trial_Power ~ inverseAge)
  print(car::Anova(lm.model))
  
  #beta
  betaTP <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.log1p_Trial_Power < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.log1p_Trial_Power < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Trial Beta Power across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.log1p_Trial_Power < betaCutoff,], Beta.log1p_Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaTP <-  (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.log1p_Trial_Power < alphaCutoff,], aes(x = age, y = Alpha.log1p_Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + xlab("Age") + ylab("log(Power)") + ggtitle("Average Trial Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "A"))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.log1p_Trial_Power < alphaCutoff,], aes(x = age, y = Alpha.log1p_Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Trial Alpha Power across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.log1p_Trial_Power < alphaCutoff,], Alpha.log1p_Trial_Power ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaTP <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.log1p_Trial_Power < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.log1p_Trial_Power < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Trial Theta Power across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.log1p_Trial_Power < thetaCutoff,], Theta.log1p_Trial_Power ~ inverseAge)
  print(anova(lm.model));
}

### Trial Power Variability

TrialPowerVariabilityRS <- function () {
  
  avgTheta <- aggregate(Theta.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdTheta <- aggregate(Theta.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdAlpha <- aggregate(Alpha.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdBeta <- aggregate(Beta.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdGamma <- aggregate(Gamma.log1p_Trial_Power_Variability ~ visitno , alldata_RS_SubLevel, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  gammaTPV <-  (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.log1p_Trial_Power_Variability < gammaCutoff,], aes(x = age, y = Gamma.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')   + xlab("Age") + ylab("log(Power) (sd)") + ggtitle("Trial Variability: Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "B"))
  
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.log1p_Trial_Power_Variability < gammaCutoff,], aes(x = age, y = Gamma.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Trial Gamma Power Variabilitly across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.log1p_Trial_Power_Variability < gammaCutoff,], Gamma.log1p_Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  betaTPV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.log1p_Trial_Power_Variability < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.log1p_Trial_Power_Variability < betaCutoff,], aes(x = age, y = Beta.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Trial Beta Power Variability across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.log1p_Trial_Power_Variability < betaCutoff,], Beta.log1p_Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaTPV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.log1p_Trial_Power_Variability < 30,], aes(x = age, y = Alpha.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')   + xlab("Age") + ylab("log(Power) (sd)") + ggtitle("Trial Variability: Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "B"))
  
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.log1p_Trial_Power_Variability < 30,], aes(x = age, y = Alpha.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Trial Alpha Power Variability across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.log1p_Trial_Power_Variability < 30,], Alpha.log1p_Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaTPV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.log1p_Trial_Power_Variability < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')   + xlab("Age") + ylab("log(Power) (sd)") + ggtitle("Trial Variability: Power") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "B"))
  
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.log1p_Trial_Power_Variability < thetaCutoff,], aes(x = age, y = Theta.log1p_Trial_Power_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Trial Theta Power Varbility across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.log1p_Trial_Power_Variability < thetaCutoff,], Theta.log1p_Trial_Power_Variability ~ inverseAge)
  print(anova(lm.model))
}

### Number of Events

NumberofEventsRS <- function() {
  
  alldf_RS <- Resting_State_Data_TrialLevel()
  alldata_RS_trialLevel <- alldf_RS$alldata_RS_TrialLevel
  
  Gamma_cleandata <-  alldata_RS_trialLevel %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2))
  
  
  #gamma 
  gammaNE <- (ggplot(data = Gamma_cleandata[], aes(x = age, y = Gamma.Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + xlab("Age") + ylab("Number of Bursts") + ggtitle("Average Number of Gamma Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "C"))
  
  
  lunaize(print(ggplot(data = Gamma_cleandata[], aes(x = age, y = Gamma.Event_Number)) + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Number of  Gamma Events across Age")  + xlab("Age") + ylab("Number of Events")))+ theme(plot.title = element_text(hjust = 0.5))
  
  lm.model <- lmer(data = Gamma_cleandata[], Gamma.Event_Number ~ inverseAge)
  print(car::Anova(lm.model))
  
  #beta
  betaNE <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Number < betaCutoff,], aes(x = age, y = Beta.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Number < betaCutoff,], aes(x = age, y = Beta.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Number of Beta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Number < betaCutoff,], Beta.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaNE <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Number < alphaCutoff,], aes(x = age, y = Alpha.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + xlab("Age") + ylab("Number of Bursts") + ggtitle("Average Number of Alpha Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "C"))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Number < alphaCutoff,], aes(x = age, y = Alpha.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Number of Alpha Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Number < alphaCutoff,], Alpha.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaNE <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Number < thetaCutoff,], aes(x = age, y = Theta.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Number < thetaCutoff,], aes(x = age, y = Theta.Event_Number))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Number of Theta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Number < thetaCutoff,], Theta.Event_Number ~ inverseAge)
  print(anova(lm.model))
  
}

### Number of Events Variabiltiy

NumberofEventsVariabilityRS <- function() {
  
  avgTheta <- aggregate(Theta.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdTheta <- aggregate(Theta.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdAlpha <- aggregate(Alpha.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdBeta <- aggregate(Beta.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdGamma <- aggregate(Gamma.Event_Number_Variability ~ visitno , alldata_RS_SubLevel, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  gammaNEV <-  (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Number_Variability < gammaCutoff,], aes(x = age, y = Gamma.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5))  + xlab("Age") + ylab("Number of Bursts (sd)") + ggtitle("Trial Variability: Number of Gamma Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "D"))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Number_Variability < gammaCutoff,], aes(x = age, y = Gamma.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Number of Gamma Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Number_Variability < gammaCutoff,], Gamma.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  betaNEV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Number_Variability < betaCutoff,], aes(x = age, y = Beta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Number_Variability < betaCutoff,], aes(x = age, y = Beta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Number of Beta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Number_Variability < betaCutoff,], Beta.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaNEV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Number_Variability < alphaCutoff,], aes(x = age, y = Alpha.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5))  + xlab("Age") + ylab("Number of Bursts (sd)") + ggtitle("Trial Variability: Number of Alpha Bursts") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs(tag= "D"))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Number_Variability < alphaCutoff,], aes(x = age, y = Alpha.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Number of Alpha Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Number_Variability < alphaCutoff,], Alpha.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaNEV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Number_Variability < thetaCutoff,], aes(x = age, y = Theta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Number_Variability < thetaCutoff,], aes(x = age, y = Theta.Event_Number_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Number of Theta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Number of Events (sd)"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Number_Variability < thetaCutoff,], Theta.Event_Number_Variability ~ inverseAge)
  print(anova(lm.model))
}

### Duration of Events

DurationofEventsRS <- function() {
  
  avgTheta <- aggregate(Theta.Event_Duration ~ visitno , alldata_RS_SubLevel, mean)
  sdTheta <- aggregate(Theta.Event_Duration ~ visitno , alldata_RS_SubLevel, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.Event_Duration ~ visitno , alldata_RS_SubLevel, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration ~ visitno , alldata_RS_SubLevel, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.Event_Duration ~ visitno , alldata_RS_SubLevel, mean)
  sdBeta <- aggregate(Beta.Event_Duration ~ visitno , alldata_RS_SubLevel, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.Event_Duration ~ visitno , alldata_RS_SubLevel, mean)
  sdGamma <- aggregate(Gamma.Event_Duration ~ visitno , alldata_RS_SubLevel, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  gammaDE <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Duration < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Duration < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Duration of  Gamma Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Duration < gammaCutoff,], Gamma.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  betaDE <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Duration < betaCutoff,], aes(x = age, y = Beta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Duration < betaCutoff,], aes(x = age, y = Beta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Duration of Beta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Duration < betaCutoff,], Beta.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaDE <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Duration < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')+ theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Duration < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Duration of Alpha Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Duration < alphaCutoff,], Alpha.Event_Duration ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaDE <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Duration < thetaCutoff,], aes(x = age, y = Theta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Duration < thetaCutoff,], aes(x = age, y = Theta.Event_Duration))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Average Duration of Theta Events across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Duration < thetaCutoff,], Theta.Event_Duration ~ inverseAge)
  print(anova(lm.model))
}

### Duration of Events Variability
DurationofEventsVariabilityRS <- function() {
  
  avgTheta <- aggregate(Theta.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdTheta <- aggregate(Theta.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, sd)
  thetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  
  avgAlpha <- aggregate(Alpha.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, sd)
  alphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  
  avgBeta <- aggregate(Beta.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdBeta <- aggregate(Beta.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, sd)
  betaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  
  avgGamma <- aggregate(Gamma.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, mean)
  sdGamma <- aggregate(Gamma.Event_Duration_Variability ~ visitno , alldata_RS_SubLevel, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  
  #gamma 
  gammaDEV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Duration_Variability < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)')  + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Duration_Variability < gammaCutoff,], aes(x = age, y = Gamma.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Duration of  Gamma Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event (sd)"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Gamma.Event_Duration_Variability < gammaCutoff,], Gamma.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #beta
  betaDEV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Duration_Variability < betaCutoff,], aes(x = age, y = Beta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Duration_Variability < betaCutoff,], aes(x = age, y = Beta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Duration of Beta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event (sd)"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Beta.Event_Duration_Variability < betaCutoff,], Beta.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Alpha
  alphaDEV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Duration_Variability < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Duration_Variability < alphaCutoff,], aes(x = age, y = Alpha.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Duration of Alpha Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Duration of Event (sd)"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Alpha.Event_Duration_Variability < alphaCutoff,], Alpha.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
  
  #Theta
  thetaDEV <- (ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Duration_Variability < thetaCutoff,], aes(x = age, y = Theta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab(""))
  
  print(ggplot(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Duration_Variability < thetaCutoff,], aes(x = age, y = Theta.Event_Duration_Variability))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Resting State: Duration of Theta Events (Variability) across Age") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = alldata_RS_SubLevel[alldata_RS_SubLevel$visitno < 2 & alldata_RS_SubLevel$Theta.Event_Duration_Variability < thetaCutoff,], Theta.Event_Duration_Variability ~ inverseAge)
  print(anova(lm.model))
}


DelayMinusRest <- function() {
  
  delayRest$Alpha.log1p_Trial_Power_delayRest <- (delayRest$Alpha.log1p_Trial_Power_Delay - delayRest$Alpha.log1p_Trial_Power)
  delayRest$Theta.log1p_Trial_Power_delayRest <- (delayRest$Theta.log1p_Trial_Power_Delay - delayRest$Theta.log1p_Trial_Power)
  delayRest$Beta.log1p_Trial_Power_delayRest <- (delayRest$Beta.log1p_Trial_Power_Delay - delayRest$Beta.log1p_Trial_Power)
  delayRest$Gamma.log1p_Trial_Power_delayRest <- (delayRest$Gamma.log1p_Trial_Power_Delay - delayRest$Gamma.log1p_Trial_Power)
  
  
  delayRest$Alpha.log1p_Trial_Power_Variability_delayRest <- (delayRest$Alpha.log1p_Trial_Power_Variability_Delay - delayRest$Alpha.log1p_Trial_Power_Variability)
  delayRest$Theta.log1p_Trial_Power_Variability_delayRest <- (delayRest$Theta.log1p_Trial_Power_Variability_Delay - delayRest$Theta.log1p_Trial_Power_Variability)
  delayRest$Beta.log1p_Trial_Power_Variability_delayRest <- (delayRest$Beta.log1p_Trial_Power_Variability_Delay - delayRest$Beta.log1p_Trial_Power_Variability)
  delayRest$Gamma.log1p_Trial_Power_Variability_delayRest <- (delayRest$Gamma.log1p_Trial_Power_Variability_Delay - delayRest$Gamma.log1p_Trial_Power_Variability)
  
  delayRest$Alpha.Event_Number_delayRest <- (delayRest$Alpha.Event_Number_Delay - delayRest$Alpha.Event_Number)
  delayRest$Theta.Event_Number_delayRest <- (delayRest$Theta.Event_Number_Delay - delayRest$Theta.Event_Number)
  delayRest$Beta.Event_Number_delayRest <- (delayRest$Beta.Event_Number_Delay - delayRest$Beta.Event_Number)
  delayRest$Gamma.Event_Number_delayRest <- (delayRest$Gamma.Event_Number_Delay - delayRest$Gamma.Event_Number)
  
  delayRest$Alpha.Event_Number_Variability_delayRest <- (delayRest$Alpha.Event_Number_Variability_Delay - delayRest$Alpha.Event_Number_Variability)
  delayRest$Theta.Event_Number_Variability_delayRest <- (delayRest$Theta.Event_Number_Variability_Delay - delayRest$Theta.Event_Number_Variability)
  delayRest$Beta.Event_Number_Variability_delayRest <- (delayRest$Beta.Event_Number_Variability_Delay - delayRest$Beta.Event_Number_Variability)
  delayRest$Gamma.Event_Number_Variability_delayRest <- (delayRest$Gamma.Event_Number_Variability_Delay - delayRest$Gamma.Event_Number_Variability)
  
  delayRest$Alpha.Event_Duration_delayRest <- (delayRest$Alpha.Event_Duration_Delay - delayRest$Alpha.Event_Duration)
  delayRest$Theta.Event_Duration_delayRest <- (delayRest$Theta.Event_Duration_Delay - delayRest$Theta.Event_Duration)
  delayRest$Beta.Event_Duration_delayRest <- (delayRest$Beta.Event_Duration_Delay - delayRest$Beta.Event_Duration)
  delayRest$Gamma.Event_Duration_delayRest <- (delayRest$Gamma.Event_Duration_Delay - delayRest$Gamma.Event_Duration)
  
  delayRest$Alpha.Event_Duration_Variability_delayRest <- (delayRest$Alpha.Event_Duration_Variability_Delay - delayRest$Alpha.Event_Duration_Variability)
  delayRest$Theta.Event_Duration_Variability_delayRest <- (delayRest$Theta.Event_Duration_Variability_Delay - delayRest$Theta.Event_Duration_Variability)
  delayRest$Beta.Event_Duration_Variability_delayRest <- (delayRest$Beta.Event_Duration_Variability_Delay - delayRest$Beta.Event_Duration_Variability)
  delayRest$Gamma.Event_Duration_Variability_delayRest <- (delayRest$Gamma.Event_Duration_Variability_Delay - delayRest$Gamma.Event_Duration_Variability)
  
  return(delayRest)
}

DelayMinusRest_TrialPower <- function () {
  
  avgGamma <- aggregate(Gamma.log1p_Trial_Power_delayRest ~ visitno , delayRest, mean)
  sdGamma <- aggregate(Gamma.log1p_Trial_Power_delayRest ~ visitno , delayRest, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.log1p_Trial_Power_delayRest ~ visitno , delayRest, mean)
  sdBeta <- aggregate(Beta.log1p_Trial_Power_delayRest ~ visitno , delayRest, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.log1p_Trial_Power_delayRest ~ visitno , delayRest, mean)
  sdAlpha <- aggregate(Alpha.log1p_Trial_Power_delayRest ~ visitno , delayRest, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.log1p_Trial_Power_delayRest ~ visitno , delayRest, mean)
  sdTheta <- aggregate(Theta.log1p_Trial_Power_delayRest ~ visitno , delayRest, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.log1p_Trial_Power_delayRest < gammaCutoff & delayRest$Gamma.log1p_Trial_Power_delayRest > gammaCutofflower,], aes(x = age, y = Gamma.log1p_Trial_Power_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.log1p_Trial_Power_delayRest < gammaCutoff & delayRest$Gamma.log1p_Trial_Power_delayRest > gammaCutofflower,], Gamma.log1p_Trial_Power_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.log1p_Trial_Power_delayRest < BetaCutoff & delayRest$Beta.log1p_Trial_Power_delayRest > BetaCutofflower,], aes(x = age, y = Beta.log1p_Trial_Power_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.log1p_Trial_Power_delayRest < BetaCutoff & delayRest$Beta.log1p_Trial_Power_delayRest > BetaCutofflower,], Beta.log1p_Trial_Power_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.log1p_Trial_Power_delayRest < AlphaCutoff & delayRest$Alpha.log1p_Trial_Power_delayRest > AlphaCutofflower,], aes(x = age, y = Alpha.log1p_Trial_Power_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.log1p_Trial_Power_delayRest < AlphaCutoff & delayRest$Alpha.log1p_Trial_Power_delayRest > AlphaCutofflower,], Alpha.log1p_Trial_Power_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.log1p_Trial_Power_delayRest < ThetaCutoff & delayRest$Theta.log1p_Trial_Power_delayRest > ThetaCutofflower,], aes(x = age, y = Theta.log1p_Trial_Power_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Trial Power") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power"))
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.log1p_Trial_Power_delayRest < ThetaCutoff & delayRest$Theta.log1p_Trial_Power_delayRest > ThetaCutofflower,], Theta.log1p_Trial_Power_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
}


DelayMinusRest_TrialPowerVar <- function () {
  
  
  avgGamma <- aggregate(Gamma.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, mean)
  sdGamma <- aggregate(Gamma.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, mean)
  sdBeta <- aggregate(Beta.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, mean)
  sdAlpha <- aggregate(Alpha.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, mean)
  sdTheta <- aggregate(Theta.log1p_Trial_Power_Variability_delayRest ~ visitno , delayRest, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.log1p_Trial_Power_Variability_delayRest < gammaCutoff & delayRest$Gamma.log1p_Trial_Power_Variability_delayRest > gammaCutofflower,], aes(x = age, y = Gamma.log1p_Trial_Power_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.log1p_Trial_Power_Variability_delayRest < gammaCutoff & delayRest$Gamma.log1p_Trial_Power_delayRest > gammaCutofflower,], Gamma.log1p_Trial_Power_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.log1p_Trial_Power_Variability_delayRest < BetaCutoff & delayRest$Beta.log1p_Trial_Power_Variability_delayRest > BetaCutofflower,], aes(x = age, y = Beta.log1p_Trial_Power_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.log1p_Trial_Power_Variability_delayRest < BetaCutoff & delayRest$Beta.log1p_Trial_Power_Variability_delayRest > BetaCutofflower,], Beta.log1p_Trial_Power_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.log1p_Trial_Power_Variability_delayRest < AlphaCutoff & delayRest$Alpha.log1p_Trial_Power_Variability_delayRest > AlphaCutofflower,], aes(x = age, y = Alpha.log1p_Trial_Power_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.log1p_Trial_Power_Variability_delayRest < AlphaCutoff & delayRest$Alpha.log1p_Trial_Power_Variability_delayRest > AlphaCutofflower,], Alpha.log1p_Trial_Power_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.log1p_Trial_Power_Variability_delayRest < ThetaCutoff & delayRest$Theta.log1p_Trial_Power_Variability_delayRest > ThetaCutofflower,], aes(x = age, y = Theta.log1p_Trial_Power_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Trial Power Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Power (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.log1p_Trial_Power_Variability_delayRest < ThetaCutoff & delayRest$Theta.log1p_Trial_Power_Variability_delayRest > ThetaCutofflower,], Theta.log1p_Trial_Power_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
}


DelayMinusRest_EventNumber <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Number_delayRest ~ visitno , delayRest, mean)
  sdGamma <- aggregate(Gamma.Event_Number_delayRest ~ visitno , delayRest, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Number_delayRest ~ visitno , delayRest, mean)
  sdBeta <- aggregate(Beta.Event_Number_delayRest ~ visitno , delayRest, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Number_delayRest ~ visitno , delayRest, mean)
  sdAlpha <- aggregate(Alpha.Event_Number_delayRest ~ visitno , delayRest, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Number_delayRest ~ visitno , delayRest, mean)
  sdTheta <- aggregate(Theta.Event_Number_delayRest ~ visitno , delayRest, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Number_delayRest < gammaCutoff & delayRest$Gamma.Event_Number_delayRest > gammaCutofflower,], aes(x = age, y = Gamma.Event_Number_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Number_delayRest < gammaCutoff & delayRest$Gamma.Event_Number_delayRest > gammaCutofflower,], Gamma.Event_Number_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Number_delayRest < BetaCutoff & delayRest$Beta.Event_Number_delayRest > BetaCutofflower,], aes(x = age, y = Beta.Event_Number_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Number_delayRest < BetaCutoff & delayRest$Beta.Event_Number_delayRest > BetaCutofflower,], Beta.Event_Number_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Number_delayRest < AlphaCutoff & delayRest$Alpha.Event_Number_delayRest > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Number_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Number_delayRest < AlphaCutoff & delayRest$Alpha.Event_Number_delayRest > AlphaCutofflower,], Alpha.Event_Number_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Number_delayRest < ThetaCutoff & delayRest$Theta.Event_Number_delayRest > ThetaCutofflower,], aes(x = age, y = Theta.Event_Number_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Number") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Number_delayRest < ThetaCutoff & delayRest$Theta.Event_Number_delayRest > ThetaCutofflower,], Theta.Event_Number_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
  
  
}


DelayMinusRest_EventNumberVar <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Number_Variability_delayRest ~ visitno , delayRest, mean)
  sdGamma <- aggregate(Gamma.Event_Number_Variability_delayRest ~ visitno , delayRest, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Number_Variability_delayRest ~ visitno , delayRest, mean)
  sdBeta <- aggregate(Beta.Event_Number_Variability_delayRest ~ visitno , delayRest, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Number_Variability_delayRest ~ visitno , delayRest, mean)
  sdAlpha <- aggregate(Alpha.Event_Number_Variability_delayRest ~ visitno , delayRest, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Number_Variability_delayRest ~ visitno , delayRest, mean)
  sdTheta <- aggregate(Theta.Event_Number_Variability_delayRest ~ visitno , delayRest, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Number_Variability_delayRest < gammaCutoff & delayRest$Gamma.Event_Number_Variability_delayRest > gammaCutofflower,], aes(x = age, y = Gamma.Event_Number_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Number_Variability_delayRest < gammaCutoff & delayRest$Gamma.Event_Number_Variability_delayRest > gammaCutofflower,], Gamma.Event_Number_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Number_Variability_delayRest < BetaCutoff & delayRest$Beta.Event_Number_Variability_delayRest > BetaCutofflower,], aes(x = age, y = Beta.Event_Number_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Number_Variability_delayRest < BetaCutoff & delayRest$Beta.Event_Number_Variability_delayRest > BetaCutofflower,], Beta.Event_Number_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Number_Variability_delayRest < AlphaCutoff & delayRest$Alpha.Event_Number_Variability_delayRest > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Number_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Number_Variability_delayRest < AlphaCutoff & delayRest$Alpha.Event_Number_Variability_delayRest > AlphaCutofflower,], Alpha.Event_Number_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Number_Variability_delayRest < ThetaCutoff & delayRest$Theta.Event_Number_Variability_delayRest > ThetaCutofflower,], aes(x = age, y = Theta.Event_Number_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Number Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Number (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Number_Variability_delayRest < ThetaCutoff & delayRest$Theta.Event_Number_Variability_delayRest > ThetaCutofflower,], Theta.Event_Number_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
  
}


DelayMinusRest_EventDuration <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Duration_delayRest ~ visitno , delayRest, mean)
  sdGamma <- aggregate(Gamma.Event_Duration_delayRest ~ visitno , delayRest, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Duration_delayRest ~ visitno , delayRest, mean)
  sdBeta <- aggregate(Beta.Event_Duration_delayRest ~ visitno , delayRest, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Duration_delayRest ~ visitno , delayRest, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration_delayRest ~ visitno , delayRest, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Duration_delayRest ~ visitno , delayRest, mean)
  sdTheta <- aggregate(Theta.Event_Duration_delayRest ~ visitno , delayRest, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Duration_delayRest < gammaCutoff & delayRest$Gamma.Event_Duration_delayRest > gammaCutofflower,], aes(x = age, y = Gamma.Event_Duration_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Duration_delayRest < gammaCutoff & delayRest$Gamma.Event_Duration_delayRest > gammaCutofflower,], Gamma.Event_Duration_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Duration_delayRest < BetaCutoff & delayRest$Beta.Event_Duration_delayRest > BetaCutofflower,], aes(x = age, y = Beta.Event_Duration_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Duration_delayRest < BetaCutoff & delayRest$Beta.Event_Duration_delayRest > BetaCutofflower,], Beta.Event_Duration_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Duration_delayRest < AlphaCutoff & delayRest$Alpha.Event_Duration_delayRest > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Duration_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Duration_delayRest < AlphaCutoff & delayRest$Alpha.Event_Duration_delayRest > AlphaCutofflower,], Alpha.Event_Duration_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Duration_delayRest < ThetaCutoff & delayRest$Theta.Event_Duration_delayRest > ThetaCutofflower,], aes(x = age, y = Theta.Event_Duration_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Duration") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Duration_delayRest < ThetaCutoff & delayRest$Theta.Event_Duration_delayRest > ThetaCutofflower,], Theta.Event_Duration_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
  
}


DelayMinusRest_EventDurationVar <- function () {
  
  avgGamma <- aggregate(Gamma.Event_Duration_Variability_delayRest ~ visitno , delayRest, mean)
  sdGamma <- aggregate(Gamma.Event_Duration_Variability_delayRest ~ visitno , delayRest, sd)
  gammaCutoff <- 2.5* sdGamma[1,2] + avgGamma[1,2]
  gammaCutofflower <- avgGamma[1,2] - 2.5* sdGamma[1,2]  
  
  avgBeta <- aggregate(Beta.Event_Duration_Variability_delayRest ~ visitno , delayRest, mean)
  sdBeta <- aggregate(Beta.Event_Duration_Variability_delayRest ~ visitno , delayRest, sd)
  BetaCutoff <- 2.5* sdBeta[1,2] + avgBeta[1,2]
  BetaCutofflower <- avgBeta[1,2] - 2.5* sdBeta[1,2] 
  
  avgAlpha <- aggregate(Alpha.Event_Duration_Variability_delayRest ~ visitno , delayRest, mean)
  sdAlpha <- aggregate(Alpha.Event_Duration_Variability_delayRest ~ visitno , delayRest, sd)
  AlphaCutoff <- 2.5* sdAlpha[1,2] + avgAlpha[1,2]
  AlphaCutofflower <- avgAlpha[1,2] - 2.5* sdAlpha[1,2]  
  
  avgTheta <- aggregate(Theta.Event_Duration_Variability_delayRest ~ visitno , delayRest, mean)
  sdTheta <- aggregate(Theta.Event_Duration_Variability_delayRest ~ visitno , delayRest, sd)
  ThetaCutoff <- 2.5* sdTheta[1,2] + avgTheta[1,2]
  ThetaCutofflower <- avgTheta[1,2] - 2.5* sdTheta[1,2]  
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Duration_Variability_delayRest < gammaCutoff & delayRest$Gamma.Event_Duration_Variability_delayRest > gammaCutofflower,], aes(x = age, y = Gamma.Event_Duration_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Gamma Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Gamma.Event_Duration_Variability_delayRest < gammaCutoff & delayRest$Gamma.Event_Duration_Variability_delayRest > gammaCutofflower,], Gamma.Event_Duration_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Duration_Variability_delayRest < BetaCutoff & delayRest$Beta.Event_Duration_Variability_delayRest > BetaCutofflower,], aes(x = age, y = Beta.Event_Duration_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Beta Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Beta.Event_Duration_Variability_delayRest < BetaCutoff & delayRest$Beta.Event_Duration_Variability_delayRest > BetaCutofflower,], Beta.Event_Duration_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Duration_Variability_delayRest < AlphaCutoff & delayRest$Alpha.Event_Duration_Variability_delayRest > AlphaCutofflower,], aes(x = age, y = Alpha.Event_Duration_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Alpha Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Alpha.Event_Duration_Variability_delayRest < AlphaCutoff & delayRest$Alpha.Event_Duration_Variability_delayRest > AlphaCutofflower,], Alpha.Event_Duration_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  print(ggplot(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Duration_Variability_delayRest < ThetaCutoff & delayRest$Theta.Event_Duration_Variability_delayRest > ThetaCutofflower,], aes(x = age, y = Theta.Event_Duration_Variability_delayRest))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Delay-Fix: Theta Event Duration Variability") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Event Duration (sd)"))
  
  lm.model <- lm(data = delayRest[delayRest$visitno < 2 & delayRest$Theta.Event_Duration_Variability_delayRest < ThetaCutoff & delayRest$Theta.Event_Duration_Variability_delayRest > ThetaCutofflower,], Theta.Event_Duration_Variability_delayRest ~ inverseAge)
  print(anova(lm.model))
  
  
  
}




