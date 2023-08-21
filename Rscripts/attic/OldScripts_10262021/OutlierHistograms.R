library(RPMG)
library(ggplot2)
library(car)
library(cowplot)

Only_Take_One_Delay_Bin_TrialLevel()
alldata_TrialLevel <- alldf$alldata_TrialLevel

Behavior <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Behavior_20200921.csv')
agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')
agefile$Subject <- agefile$idvalues 

#behavioral data contains each trial 
Behavior$Subject <- paste0(Behavior$LunaID, '_', Behavior$ScanDate)
Behavior$absPositionError <- abs(Behavior$PositionError)
Behavior <- Behavior %>% filter(absPositionError < 40)  # ANYTHING MORE THAN 40 DEGS WOULD BE OFF THE SCREEN 


alldata_TrialLevel_behavior <- merge(alldata_TrialLevel, Behavior, by = c("Subject", "Trial"))


hist2d <- function(mydata, xvar, yvar, groupvar='Subject', xlabel=xvar, ylabel=yvar) {
  
  PEcriteria <- mean(mydata$PositionError, na.rm= T) + sd(mydata$PositionError, na.rm = T) * 2
  gammaCriteria <- mean(mydata$Gamma.Gamma_Event_Number, na.rm= T) + sd(mydata$Gamma.Gamma_Event_Number, na.rm = T) * 2
  
  mydata_outlier <- mydata %>% mutate(isgood = absPositionError < PEcriteria & Gamma.Gamma_Event_Number < gammaCriteria)
  
  sp <- ggscatter(mydata_outlier, x = xvar, y = yvar, color = 'isgood',
                  palette = "jco",
                  size = 1.5, alpha = 0.6)+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(x=xlabel, y=ylabel)+
    border()                                         
  # Marginal density plot of x (top panel) and y (right panel)
  xplot <- ggdensity(sp$data, xvar, fill = 'isgood',
                     palette = "jco")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  
  yplot <- ggdensity(sp$data, yvar, fill = 'isgood',
                     palette = "jco")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    rotate()
  # Cleaning the plots
  sp <- sp + rremove("legend")
  yplot <- yplot + clean_theme() + rremove("legend")
  xplot <- xplot + clean_theme() + rremove("legend")
  # Arranging the plot using cowplot
  p <- cowplot::plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
                     rel_widths = c(3, 1), rel_heights = c(1, 4))
  
  name <- paste(d,".png")
  png(name)
  print(p)
  dev.off()

  
}



mydata <- alldata_TrialLevel_behavior
xvar <- 'absPositionError'
yvar <- 'Gamma.Gamma_Event_Number'

for (d in unique(mydata$Subject)) {
  print(d)
  print(hist2d(mydata[mydata$Subject==d,],  xvar, yvar, groupvar='Subject', xlabel=xvar, ylabel=yvar))
}

