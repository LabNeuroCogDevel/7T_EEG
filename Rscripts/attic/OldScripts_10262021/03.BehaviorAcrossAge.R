
sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
subjectBehavior_new <- Behavior_Sublevel_Maria()


# All run on the subject level 
absPositionError <- function (subjectBehavior_new) {
  

 absPE <- (ggplot(data = subjectBehavior_new, aes(x = age, y = absPositionError))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Accuracy")) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),  axis.title = element_text(size = 25), axis.text = element_text(size = 20),  axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)), plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Age") + ylab("Accuracy (degs)") + labs(tag= "A")
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absPositionError))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Accuracy across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Accuracy (degs)"))
  
  lm.model <- lm(data = subjectBehavior_new[], absPositionError ~ inverseAge)
  print(anova(lm.model))
}

absPositionErrorVariability <- function (subjectBehavior_new) {
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absPositionError_sd))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Accuracy Variability across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("SD of Accuracy (degs)"))
  
 absPEV <- (ggplot(data = subjectBehavior_new[], aes(x = age, y = absPositionError_sd))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Variability: Accuracy") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),  axis.title = element_text(size = 25), axis.text = element_text(size = 20),  axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)), plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Age") + ylab("SD of Accuracy (degs)"))+ labs(tag= "B")
  
  lm.model <- lm(data = subjectBehavior_new[], absPositionError_sd ~ inverseAge)
  print(anova(lm.model))
}



mgsLatency <- function (subjectBehavior_new) {
  

  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("MGS Latency across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Latency"))
  
 mgsLat <- (ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Response Latency") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)), plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Age") + ylab("Latency (sec)"))+ labs(tag= "C")
  
  lm.model <- lm(data = subjectBehavior_new[], mgsLatency ~ inverseAge)
  print(anova(lm.model))
  
}

mgsLatencyVariability <- function (subjectBehavior_new) {
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency_sd))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("MGS Latency Variability across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Latency (SD)"))
  
  mgsLatV <- (ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency_sd))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Variability: Response Latency") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), axis.title = element_text(size = 25), axis.text = element_text(size = 20), axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)),  plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Age") + ylab("SD of Latency (sec)")) + labs(tag= "D") 
  
  lm.model <- lm(data = subjectBehavior_new[], mgsLatency_sd ~ inverseAge)
  print(anova(lm.model))
}


# library(gridExtra)
# grid.arrange(absPE, absPEV, mgsLat, mgsLatV, nrow = 2, widths = c(1,1), heights = c(2,2), layout_matrix = rbind(c(1,2),c(3,4)))
# 
# 
