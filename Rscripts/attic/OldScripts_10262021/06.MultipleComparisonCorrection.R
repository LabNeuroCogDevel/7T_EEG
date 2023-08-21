
#Multiple Comparison Correction 

Pvalues <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/FrequencyAcrossAge_Pvalues.csv')
Pvalues <- subset(Pvalues[1])
Pvalues <- data.matrix(Pvalues)

newP <- p.adjust(Pvalues, method = "fdr", n = length(Pvalues))
