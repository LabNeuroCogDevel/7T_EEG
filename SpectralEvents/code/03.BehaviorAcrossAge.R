
sys.source("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
subjectBehavior_new <- Behavior_Sublevel_Maria()


# All run on the subject level 
absPositionError <- function (subjectBehavior_new) {
  

 absPE <- (ggplot(data = subjectBehavior_new, aes(x = age, y = absPositionError))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Average Accuracy")) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),  axis.title = element_text(size = 25), axis.text = element_text(size = 20),  axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)), plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Age") + ylab("Accuracy (degs)") + labs(tag= "A")
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absPositionError))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Accuracy across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Accuracy (degs)"))
  
  absPE.model <- lm(data = subjectBehavior_new[], absPositionError ~ inverseAge)
  print(anova(absPE.model))
  
  anovadf <- anova(absPE.model)
  PE <- list()
  PE$Pvalue <- anovadf$`Pr(>F)`
  PE$PEgraph <- absPE
  
  return(PE)
}

absPositionErrorVariability <- function (subjectBehavior_new) {
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absPositionError_sd))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Accuracy Variability across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("SD of Accuracy (degs)"))
  
 absPEV <- (ggplot(data = subjectBehavior_new[], aes(x = age, y = absPositionError_sd))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("Trial Variability: Accuracy") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),  axis.title = element_text(size = 25), axis.text = element_text(size = 20),  axis.title.y = element_text(margin = margin(t=0, r= 10, b= 0, l= 0)), plot.tag = element_text(size = 30, face = "bold"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Age") + ylab("SD of Accuracy (degs)"))+ labs(tag= "B")
  
 absPEV.model <- lm(data = subjectBehavior_new[], absPositionError_sd ~ inverseAge)
  print(anova(absPEV.model))
  
  anovadf <- anova(absPEV.model)
  PEV <- list()
  PEV$Pvalue <- anovadf$`Pr(>F)`
  PEV$PEVgraph <- absPEV
  
  return(PEV)
}



mgsLatency <- function (subjectBehavior_new) {
  

  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("MGS Latency across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Latency"))
  scaleFUN <- function(x) sprintf("%.2f", x)
  
 mgsLat <- lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency))  + geom_point() + stat_smooth(method = 'gam') + xlab("") + ylab("")+ scale_y_continuous(labels = scaleFUN))
  
 lat.model <- lm(data = subjectBehavior_new[], mgsLatency ~ inverseAge)
  print(anova(lat.model))
  
  anovadf <- anova(lat.model)
  lat <- list()
  lat$Pvalue <- anovadf$`Pr(>F)`
  lat$latgraph <- mgsLat
  
  return(lat)
  
}

mgsLatencyVariability <- function (subjectBehavior_new) {
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency_sd))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + ggtitle("MGS Latency Variability across Age")) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Age") + ylab("Latency (SD)"))
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  mgsLatV <- (lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = mgsLatency_sd))  + geom_point() + stat_smooth(method = 'gam') + xlab("") + ylab("")) + scale_y_continuous(labels = scaleFUN))
  
  latV.model <- lm(data = subjectBehavior_new[], mgsLatency_sd ~ inverseAge)
  print(anova(latV.model))
  
  anovadf <- anova(latV.model)
  latv <- list()
  latv$Pvalue <- anovadf$`Pr(>F)`
  latv$latvgraph <- mgsLatV
  
  return(latv)
}

bestSaccade <- function(subjectBehavior_new) {
  
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absBestError))  + geom_point() + stat_smooth(method = 'lm',  formula='y~I(1/x)') + theme_classic()))
  
  BSgraph <- (lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absBestError))  + geom_point() + stat_smooth(method = 'gam') + xlab("") + ylab(""))+ theme(plot.margin = margin(1,1,1,1)))
  
  
  BS.model <- lm(data = subjectBehavior_new[], absBestError ~ inverseAge)
  print(anova(BS.model))
  
  anovadf <- anova(BS.model)
  BS <- list()
  BS$Pvalue <- anovadf$`Pr(>F)`
  BS$BSgraph <- BSgraph
  
  return(BS)
  
}

bestSaccadeVariability <- function(subjectBehavior_new) {
  
  
  print(lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absBestError_sd))  + geom_point() + stat_smooth(method = 'gam') + theme_classic()))
  
  BSVgraph <- (lunaize(ggplot(data = subjectBehavior_new[], aes(x = age, y = absBestError_sd))  + geom_point() + stat_smooth(method = 'gam') + xlab("") + ylab("")))
  
  BSV.model <- lm(data = subjectBehavior_new[], absBestError_sd ~ inverseAge)
  print(anova(BSV.model))
  
  anovadf <- anova(BSV.model)
  BSV <- list()
  BSV$Pvalue <- anovadf$`Pr(>F)`
  BSV$BSVgraph <- BSVgraph
  

  return(BSV)
}


# library(gridExtra)
# grid.arrange(absPE, absPEV, mgsLat, mgsLatV, nrow = 2, widths = c(1,1), heights = c(2,2), layout_matrix = rbind(c(1,2),c(3,4)))
# 
# 
