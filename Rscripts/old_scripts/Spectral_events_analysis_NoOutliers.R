
gammapeakPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_all_allevents_NoOutliers_peakPower.csv')
gammaavgPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_all_allevents_NoOutliers_AvgPower.csv')
gammaEventNumber <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_all_allevents_NoOutliers_NumberofEvents.csv')
gammaDuration <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_all_allevents_NoOutliers_MeanEventDuration.csv')


Gamma_Spectral_Analysis_Table <- merge(gammapeakPower, gammaavgPower,  by='idvalues')
Gamma_Spectral_Analysis_Table2 <- merge(gammaEventNumber, gammaDuration,  by='idvalues')
Final_Gamma = merge(Gamma_Spectral_Analysis_Table, Gamma_Spectral_Analysis_Table2, by= 'idvalues')

agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')

Gamma_Spectral_Analysis_Table_all = merge(agefile, Final_Gamma, by= 'idvalues')
write.csv(Gamma_Spectral_Analysis_Table_all, "H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_all_allevents.csv")



thetapeakPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_all_allevents_NoOutliers_peakPower.csv')
thetaavgPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_all_allevents_NoOutliers_AvgPower.csv')
thetaEventNumber <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_all_allevents_NoOutliers_NumberofEvents.csv')
ThetaDuration <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_all_allevents_NoOutliers_MeanEventDuration.csv')


Theta_Spectral_Analysis_Table <- merge(thetapeakPower, thetaavgPower,  by='idvalues')
Theta_Spectral_Analysis_Table2 <- merge(thetaEventNumber, ThetaDuration,  by='idvalues')
Final_Theta = merge(Theta_Spectral_Analysis_Table, Theta_Spectral_Analysis_Table2, by= 'idvalues')

agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')

Theta_Spectral_Analysis_Table_all = merge(agefile, Final_Theta, by= 'idvalues')
write.csv(Theta_Spectral_Analysis_Table_all, "H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_all_allevents.csv")




betapeakPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_all_allevents_NoOutliers_peakPower.csv')
betaavgPower <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_all_allevents_NoOutliers_AvgPower.csv')
betaEventNumber <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_all_allevents_NoOutliers_NumberofEvents.csv')
betaDuration <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_all_allevents_NoOutliers_MeanEventDuration.csv')


Beta_Spectral_Analysis_Table <- merge(betapeakPower, betaavgPower,  by='idvalues')
Beta_Spectral_Analysis_Table2 <- merge(betaEventNumber, betaDuration,  by='idvalues')
Final_Beta = merge(Beta_Spectral_Analysis_Table, Beta_Spectral_Analysis_Table2, by= 'idvalues')



agefile <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv')

Beta_Spectral_Analysis_Table_all = merge(agefile, Final_Beta, by= 'idvalues')
write.csv(Beta_Spectral_Analysis_Table_all, "H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_all_allevents.csv")
