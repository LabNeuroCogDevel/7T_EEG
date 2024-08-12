
# Looking into relationships with MRS data 

library(readxl)
library(tidyr)
library(dplyr)

prepMRSdata <- function() {

MRS <- read.csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/13MP20200207_LCMv2fixidx.csv')

# Step 1 Outlier Detection - visual inspection of LCModel fits/spectra
# create a list of who to remove and remove them
lcm <- read_excel('H:/Projects/7TBrainMech/scripts/eeg/Shane/MRS/lcm.xlsx', col_names = FALSE)

lcm <- separate(lcm, "...1", c("ld8", "junk","y","x"),extra="merge", sep = "[-.]")
lcm <- dplyr::select(lcm, -junk)
lcm$bad <- TRUE
MRS <- MRS %>% mutate(x=216+1-x,y=216+1-y)
MRS <- merge(MRS, lcm, by=c("ld8", "x", "y"), all=T) 
MRS <- filter(MRS, is.na(bad))
MRS <- dplyr::select(MRS, -bad)

#keep only visit 1 people
MRS <- MRS %>% filter(visitnum==1)
#keep people's correct coordinates
MRS <- MRS %>% filter(!is.na(roi))
#get rid of people who are actually visit 2 but for some reason aren't filtered out
MRS <- MRS %>% filter(ld8!="10195_20191205")

# Step 2 Outlier Detection - get rid of peole who have bad data for 3 major metabolite peaks - GPC+Cho, NAA+NAAG, Cr
MRS<- filter(MRS, GPC.Cho.SD <= 10 | is.na(GPC.Cho.SD))
MRS <- filter(MRS, NAA.NAAG.SD <= 10 | is.na(NAA.NAAG.SD))
MRS <- filter(MRS, Cr.SD <= 10 | is.na(Cr.SD))

# Step 3 Outlier Detection - get rid of people who have lots of macromolecule in their spectra, as that can create distortions
MRS <- filter(MRS, MM20.Cr <= 3 | is.na(MM20.Cr))
#make inverse age column
MRS$invage <- 1/MRS$age
#make age^2 column
MRS$age2 <- (MRS$age - mean(MRS$age))^2
z_thres = 2

# Create dataframe with good quality Glutamate data 
MRS_glu <- MRS %>% filter(Glu.Cr != 0)
MRS_glu <- MRS_glu %>% filter(Glu.SD <=20)

# Create dataframe with good quality GABA data 
MRS_GABA <- MRS %>% filter(GABA.SD <=20)
}

# Load in eeg data, subject level 
sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)

all_df <- DelayOnly_Sublevel()
alldata_delayOnly <- all_df$alldata_delayOnly
alldata_delayOnly$Subject <- as.character(alldata_delayOnly$Subject)
alldata_delayOnly <- alldata_delayOnly %>% separate(Subject, c("ID", "Date"), sep = "_" )

# Looking at max freq vs gaba 
alldf_peaks <- DelayOnly_peakFreq()
gammaPeaks <- alldf_peaks$delayGammaPeaks

# combine all eeg data together 
alldata_delayOnly <- merge(alldata_delayOnly, gammaPeaks, by = c("idvalues", "age"))


# mediation age, power, GABA
source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/05.MediationAnalysis.R")

other_vars <-c("GABA", "Gamma.log1p_Trial_Power") #define what you want in the mediation analysis 

EEG_MRS_data_ACC_cleanPowerGABA_Subset <- subset(EEG_MRS_data_ACC_cleanPowerGABA, select = c(ID,Date.x, GABA, age.x, Gamma.log1p_Trial_Power) )
EEG_MRS_data_ACC_cleanPowerGABA_Subset$Subject <- EEG_MRS_data_ACC_cleanPowerGABA$idvalues

MediationAnalysis_SubjectLevel_Gamma(EEG_MRS_data_ACC_cleanPowerGABA_Subset, other_vars)

other_vars <-c("GABA", "Gamma.log1p_Trial_Power") #define what you want in the mediation analysis 
MediationAnalysis_SubjectLevel_Gamma(EEG_MRS_data_ACC_cleanPowerGABA_Subset, other_vars)

WholeBrainEEG_vs_DLPFC_MRS <- function () {

# DLPFC ROI
DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
DLPFC_GABA <- DLPFC_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
DLPFC_GABA <- DLPFC_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)
DLPFC_GABA$Subject <- DLPFC_GABA$ld8
DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)

EEG_MRS_data_DLPFC <- merge(alldata_delayOnly, Avg_DLPFC_GABA, by = "ID")
EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)

# Gamma number of events vs GABA
EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = Gamma.Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
summary(lm.model)

# Gamma number of events variability vs GABA
EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA, aes(x = Gamma.Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA)
summary(lm.model)

# Gamma power vs GABA
EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = Gamma.log1p_Trial_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
summary(lm.model)

# Gamma power variability vs GABA
EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA, aes(x = Gamma.log1p_Trial_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA)
summary(lm.model)

# gamma max frequency vs GABA
EEG_MRS_data_DLPFC_cleanFreqGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Peak_Frequency - mean(Gamma.Peak_Frequency, na.rm= T)) < (sd(Gamma.Peak_Frequency, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanFreqGABA, aes(x = Gamma.Peak_Frequency, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Peak_Frequency + age, data = EEG_MRS_data_DLPFC_cleanFreqGABA)
summary(lm.model)

# gamma max frequency variability vs GABA
EEG_MRS_data_DLPFC_cleanFreqVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Peak_Frequency_Variability - mean(Gamma.Peak_Frequency_Variability, na.rm= T)) < (sd(Gamma.Peak_Frequency_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanFreqVarGABA, aes(x = Gamma.Peak_Frequency_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Peak_Frequency_Variability + age, data = EEG_MRS_data_DLPFC_cleanFreqVarGABA)
summary(lm.model)

# DLPFC Glutamate 
DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
DLPFC_Glu <- DLPFC_Glu%>%
  mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
DLPFC_Glu <- DLPFC_Glu %>% 
  filter(abs(zscore_Glu) <= z_thres)
DLPFC_Glu$Subject <- DLPFC_Glu$ld8
DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)

EEG_MRS_data_DLPFC <- merge(alldata_delayOnly, Avg_DLPFC_Glu, by = "ID")
EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group.y)


# Gamma number of events vs Glu
EEG_MRS_data_DLPFC_cleanEventGlu <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventGlu, aes(x = Gamma.Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGlu)
summary(lm.model)

# Gamma number of events variability vs Glu
EEG_MRS_data_DLPFC_cleanEventVarGlu <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGlu, aes(x = Gamma.Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGlu)
summary(lm.model)


# Gamma power vs Glu
EEG_MRS_data_DLPFC_cleanPowerGlu <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGlu, aes(x = Gamma.log1p_Trial_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGlu)
summary(lm.model)


# Gamma power variability vs Glu
EEG_MRS_data_DLPFC_cleanPowerVarGlu <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGlu, aes(x = Gamma.log1p_Trial_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGlu)
summary(lm.model)

# gamma max frequency vs glu
EEG_MRS_data_DLPFC_cleanFreqGlu <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Peak_Frequency - mean(Gamma.Peak_Frequency, na.rm= T)) < (sd(Gamma.Peak_Frequency, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanFreqGlu, aes(x = Gamma.Peak_Frequency, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Peak_Frequency + age, data = EEG_MRS_data_DLPFC_cleanFreqGlu)
summary(lm.model)


# gamma max frequency variability vs glu
EEG_MRS_data_DLPFC_cleanFreqVarGlu <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Peak_Frequency_Variability - mean(Gamma.Peak_Frequency_Variability, na.rm= T)) < (sd(Gamma.Peak_Frequency_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanFreqVarGlu, aes(x = Gamma.Peak_Frequency_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Peak_Frequency_Variability + age, data = EEG_MRS_data_DLPFC_cleanFreqVarGlu)
summary(lm.model)

# looking at glu/gaba
# DLPFC ROI
DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
DLPFC_GABA <- DLPFC_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
DLPFC_GABA <- DLPFC_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)
DLPFC_GABA$Subject <- DLPFC_GABA$ld8
DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)

EEG_MRS_data_DLPFC <- merge(alldata_delayOnly, Avg_DLPFC_GABA, by = "ID")
EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)


# DLPFC Glutamate 
DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
DLPFC_Glu <- DLPFC_Glu%>%
  mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
DLPFC_Glu <- DLPFC_Glu %>% 
  filter(abs(zscore_Glu) <= z_thres)
DLPFC_Glu$Subject <- DLPFC_Glu$ld8
DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)

EEG_MRS_data_DLPFC <- merge(EEG_MRS_data_DLPFC, Avg_DLPFC_Glu, by = "ID")

EEG_MRS_data_DLPFC$gluGABAratio <- EEG_MRS_data_DLPFC$Glu.Cr / EEG_MRS_data_DLPFC$GABA.Cr

# Gamma number of events vs glu/GABA
EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = Gamma.Event_Number, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
summary(lm.model)

# Gamma number of events variability vs glu/GABA
EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA, aes(x = Gamma.Event_Number_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA)
summary(lm.model)

# Gamma power vs glu/GABA
EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = Gamma.log1p_Trial_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
summary(lm.model)

# Gamma power variability vs glu/GABA
EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA, aes(x = Gamma.log1p_Trial_Power_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA)
summary(lm.model)

}

WholeBrainEEG_vs_MPFC_MRS <- function () {
# MPFC ROI
MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
MPFC_GABA <- MPFC_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
MPFC_GABA <- MPFC_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)
MPFC_GABA$Subject <- MPFC_GABA$ld8
MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)

EEG_MRS_data_MPFC <- merge(alldata_delayOnly, Avg_MPFC_GABA, by = "ID")
EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group.y)

# Gamma number of events vs GABA
EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma.Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
summary(lm.model)

# Gamma number of events variability vs GABA
EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA, aes(x = Gamma.Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA)
summary(lm.model)

# Gamma power vs GABA
EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma.log1p_Trial_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
summary(lm.model)

# Gamma power variability vs GABA
EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA, aes(x = Gamma.log1p_Trial_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA)
summary(lm.model)


# gamma max frequency vs GABA
EEG_MRS_data_MPFC_cleanFreqGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Peak_Frequency - mean(Gamma.Peak_Frequency, na.rm= T)) < (sd(Gamma.Peak_Frequency, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanFreqGABA, aes(x = Gamma.Peak_Frequency, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Peak_Frequency + age, data = EEG_MRS_data_MPFC_cleanFreqGABA)
summary(lm.model)

# gamma max frequency variability vs GABA
EEG_MRS_data_MPFC_cleanFreqVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Peak_Frequency_Variability - mean(Gamma.Peak_Frequency_Variability, na.rm= T)) < (sd(Gamma.Peak_Frequency_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanFreqVarGABA, aes(x = Gamma.Peak_Frequency_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Peak_Frequency_Variability + age, data = EEG_MRS_data_MPFC_cleanFreqVarGABA)
summary(lm.model)


# MPFC Glutamate 

MPFC_Glu <- MRS_glu %>% filter(roi == 8)
MPFC_Glu <- MPFC_Glu%>%
  mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
MPFC_Glu <- MPFC_Glu %>% 
  filter(abs(zscore_Glu) <= z_thres)
MPFC_Glu$Subject <- MPFC_Glu$ld8
MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)

EEG_MRS_data_MPFC <- merge(alldata_delayOnly, Avg_MPFC_Glu, by = "ID")
EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group.y)


# Gamma number of events vs Glu
EEG_MRS_data_MPFC_cleanEventGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventGlu, aes(x = Gamma.Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGlu)
summary(lm.model)

# Gamma number of events variability vs Glu
EEG_MRS_data_MPFC_cleanEventVarGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGlu, aes(x = Gamma.Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGlu)
summary(lm.model)


# Gamma power vs Glu
EEG_MRS_data_MPFC_cleanPowerGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventGlu, aes(x = Gamma.log1p_Trial_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_MPFC_cleanEventGlu)
summary(lm.model)


# Gamma power variability vs Glu
EEG_MRS_data_MPFC_cleanPowerVarGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGlu, aes(x = Gamma.log1p_Trial_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGlu)
summary(lm.model)


# gamma max frequency vs glu
EEG_MRS_data_MPFC_cleanFreqGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Peak_Frequency - mean(Gamma.Peak_Frequency, na.rm= T)) < (sd(Gamma.Peak_Frequency, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanFreqGlu, aes(x = Gamma.Peak_Frequency, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Peak_Frequency + age, data = EEG_MRS_data_MPFC_cleanFreqGlu)
summary(lm.model)


# gamma max frequency variability vs glu
EEG_MRS_data_MPFC_cleanFreqVarGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Peak_Frequency_Variability - mean(Gamma.Peak_Frequency_Variability, na.rm= T)) < (sd(Gamma.Peak_Frequency_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanFreqVarGlu, aes(x = Gamma.Peak_Frequency_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Peak_Frequency_Variability + age, data = EEG_MRS_data_MPFC_cleanFreqVarGlu)
summary(lm.model)

# looking at glu/gaba
# MPFC ROI
MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
MPFC_GABA <- MPFC_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
MPFC_GABA <- MPFC_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)
MPFC_GABA$Subject <- MPFC_GABA$ld8
MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)

EEG_MRS_data_MPFC <- merge(alldata_delayOnly, Avg_MPFC_GABA, by = "ID")
EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)


# MPFC Glutamate 
MPFC_Glu <- MRS_glu %>% filter(roi == 8)
MPFC_Glu <- MPFC_Glu%>%
  mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
MPFC_Glu <- MPFC_Glu %>% 
  filter(abs(zscore_Glu) <= z_thres)
MPFC_Glu$Subject <- MPFC_Glu$ld8
MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)

EEG_MRS_data_MPFC <- merge(EEG_MRS_data_MPFC, Avg_MPFC_Glu, by = "ID")

EEG_MRS_data_MPFC$gluGABAratio <- EEG_MRS_data_MPFC$Glu.Cr / EEG_MRS_data_MPFC$GABA.Cr

# Gamma number of events vs glu/GABA
EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma.Event_Number, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
summary(lm.model)

# Gamma number of events variability vs glu/GABA
EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA, aes(x = Gamma.Event_Number_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA)
summary(lm.model)

# Gamma power vs glu/GABA
EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma.log1p_Trial_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
summary(lm.model)

# Gamma power variability vs glu/GABA
EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA, aes(x = Gamma.log1p_Trial_Power_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA)
summary(lm.model)
}

WholeBrainEEG_vs_ACC_MRS <- function () {
# ACC ROI
ACC_GABA <- MRS_GABA %>% filter(roi == 7)
ACC_GABA <- ACC_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ACC_GABA <- ACC_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)
ACC_GABA$Subject <- ACC_GABA$ld8
ACC_GABA <- ACC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_ACC_GABA <- aggregate(GABA.Cr ~ ID, data = ACC_GABA, mean)

EEG_MRS_data_ACC <- merge(alldata_delayOnly, Avg_ACC_GABA, by = "ID")
EEG_MRS_data_ACC$Group <- factor(EEG_MRS_data_ACC$Group.y)

# Gamma number of events vs GABA
EEG_MRS_data_ACC_cleanEventGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventGABA, aes(x = Gamma.Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Event_Number + age, data = EEG_MRS_data_ACC_cleanEventGABA)
summary(lm.model)

# Gamma number of events variability vs GABA
EEG_MRS_data_ACC_cleanEventVarGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventVarGABA, aes(x = Gamma.Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_ACC_cleanEventVarGABA)
summary(lm.model)

# Gamma power vs GABA
EEG_MRS_data_ACC_cleanPowerGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventGABA, aes(x = Gamma.log1p_Trial_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_ACC_cleanEventGABA)
summary(lm.model)

# Gamma power variability vs GABA
EEG_MRS_data_ACC_cleanPowerVarGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanPowerVarGABA, aes(x = Gamma.log1p_Trial_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_ACC_cleanPowerVarGABA)
summary(lm.model)


# gamma max frequency vs GABA
EEG_MRS_data_ACC_cleanFreqGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Peak_Frequency - mean(Gamma.Peak_Frequency, na.rm= T)) < (sd(Gamma.Peak_Frequency, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanFreqGABA, aes(x = Gamma.Peak_Frequency, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Peak_Frequency + age, data = EEG_MRS_data_ACC_cleanFreqGABA)
summary(lm.model)

# gamma max frequency variability vs GABA
EEG_MRS_data_ACC_cleanFreqVarGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Peak_Frequency_Variability - mean(Gamma.Peak_Frequency_Variability, na.rm= T)) < (sd(Gamma.Peak_Frequency_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanFreqVarGABA, aes(x = Gamma.Peak_Frequency_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(GABA.Cr ~ Gamma.Peak_Frequency_Variability + age, data = EEG_MRS_data_ACC_cleanFreqVarGABA)
summary(lm.model)

# ACC Glutamate 

ACC_Glu <- MRS_glu %>% filter(roi == 7)
ACC_Glu <- ACC_Glu%>%
  mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ACC_Glu <- ACC_Glu %>% 
  filter(abs(zscore_Glu) <= z_thres)
ACC_Glu$Subject <- ACC_Glu$ld8
ACC_Glu <- ACC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_ACC_Glu <- aggregate(Glu.Cr ~ ID, data = ACC_Glu, mean)

EEG_MRS_data_ACC <- merge(alldata_delayOnly, Avg_ACC_Glu, by = "ID")
EEG_MRS_data_ACC$Group <- factor(EEG_MRS_data_ACC$Group.y)


# Gamma number of events vs Glu
EEG_MRS_data_ACC_cleanEventGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventGlu, aes(x = Gamma.Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Event_Number + age, data = EEG_MRS_data_ACC_cleanEventGlu)
summary(lm.model)

# Gamma number of events variability vs Glu
EEG_MRS_data_ACC_cleanEventVarGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventVarGlu, aes(x = Gamma.Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_ACC_cleanEventVarGlu)
summary(lm.model)


# Gamma power vs Glu
EEG_MRS_data_ACC_cleanPowerGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventGlu, aes(x = Gamma.log1p_Trial_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_ACC_cleanEventGlu)
summary(lm.model)


# Gamma power variability vs Glu
EEG_MRS_data_ACC_cleanPowerVarGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanPowerVarGlu, aes(x = Gamma.log1p_Trial_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_ACC_cleanPowerVarGlu)
summary(lm.model)


# gamma max frequency vs glu
EEG_MRS_data_ACC_cleanFreqGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Peak_Frequency - mean(Gamma.Peak_Frequency, na.rm= T)) < (sd(Gamma.Peak_Frequency, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanFreqGlu, aes(x = Gamma.Peak_Frequency, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Peak_Frequency + age, data = EEG_MRS_data_ACC_cleanFreqGlu)
summary(lm.model)


# gamma max frequency variability vs glu
EEG_MRS_data_ACC_cleanFreqVarGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Peak_Frequency_Variability - mean(Gamma.Peak_Frequency_Variability, na.rm= T)) < (sd(Gamma.Peak_Frequency_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanFreqVarGlu, aes(x = Gamma.Peak_Frequency_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(Glu.Cr ~ Gamma.Peak_Frequency_Variability + age, data = EEG_MRS_data_ACC_cleanFreqVarGlu)
summary(lm.model)

# looking at glu/gaba
# ACC ROI
ACC_GABA <- MRS_GABA %>% filter(roi == 7)
ACC_GABA <- ACC_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ACC_GABA <- ACC_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)
ACC_GABA$Subject <- ACC_GABA$ld8
ACC_GABA <- ACC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_ACC_GABA <- aggregate(GABA.Cr ~ ID, data = ACC_GABA, mean)

EEG_MRS_data_ACC <- merge(alldata_delayOnly, Avg_ACC_GABA, by = "ID")
EEG_MRS_data_ACC$Group <- factor(EEG_MRS_data_ACC$Group)


# ACC Glutamate 
ACC_Glu <- MRS_glu %>% filter(roi == 7)
ACC_Glu <- ACC_Glu%>%
  mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ACC_Glu <- ACC_Glu %>% 
  filter(abs(zscore_Glu) <= z_thres)
ACC_Glu$Subject <- ACC_Glu$ld8
ACC_Glu <- ACC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )

Avg_ACC_Glu <- aggregate(Glu.Cr ~ ID, data = ACC_Glu, mean)

EEG_MRS_data_ACC <- merge(EEG_MRS_data_ACC, Avg_ACC_Glu, by = "ID")

EEG_MRS_data_ACC$gluGABAratio <- EEG_MRS_data_ACC$Glu.Cr / EEG_MRS_data_ACC$GABA.Cr

# Gamma number of events vs glu/GABA
EEG_MRS_data_ACC_cleanEventGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventGABA, aes(x = Gamma.Event_Number, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.Event_Number + age, data = EEG_MRS_data_ACC_cleanEventGABA)
summary(lm.model)

# Gamma number of events variability vs glu/GABA
EEG_MRS_data_ACC_cleanEventVarGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventVarGABA, aes(x = Gamma.Event_Number_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_ACC_cleanEventVarGABA)
summary(lm.model)

# Gamma power vs glu/GABA
EEG_MRS_data_ACC_cleanPowerGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanEventGABA, aes(x = Gamma.log1p_Trial_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_ACC_cleanEventGABA)
summary(lm.model)

# Gamma power variability vs glu/GABA
EEG_MRS_data_ACC_cleanPowerVarGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))

ggplot(data = EEG_MRS_data_ACC_cleanPowerVarGABA, aes(x = Gamma.log1p_Trial_Power_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")

lm.model <- lm(gluGABAratio ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_ACC_cleanPowerVarGABA)
summary(lm.model)
}

WholeBrainEEG_vs_PCA_MRS <- function () {
  
  all_df <- DelayOnly_Sublevel()
  alldata_delayOnly <- all_df$alldata_delayOnly
  alldata_delayOnly$Subject <- as.character(alldata_delayOnly$Subject)
  alldata_delayOnly <- alldata_delayOnly %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  
  library(missMDA) 
  library(FactoMineR)
  #### PCA ####
  # Glutamate 
  Glu_wide_PCA <- pivot_wider(MRS_glu, id_cols=ld8, names_from=label, values_from=Glu.Cr)
  Glu_wide_PCA <- select(Glu_wide_PCA, -`R STS`, -`L STS`, -`R Thalamus`, -`R Caudate`, -`L Caudate`, -`R Posterior Insula`, -`L Posterior Insula`)
  nb_glu_kfold <- estim_ncpPCA(Glu_wide_PCA %>% select(-ld8),method.cv = "Kfold", verbose = FALSE) 
  nb_glu_kfold$ncp
  nb_glu_kfold$ncp <- 2
  
  plot(0:5, nb_glu_kfold$criterion, xlab = "nb dim", ylab = "MSEP")
  
  res_comp_glu <- imputePCA(Glu_wide_PCA %>% select(-ld8) %>% as.data.frame, ncp = nb_glu_kfold$ncp) # iterativePCA algorithm
  
  glu_pca <- PCA(res_comp_glu, ncp = nb_glu_kfold$ncp)
  plot(glu_pca)
  plot(glu_pca, choix="var")
  Glu_wide_PCA$dim1 <- glu_pca$ind$coord[,1]
  dim1_lai <- lm(data=Glu_wide_PCA, `MPFC` ~ dim1)
  summary(dim1_lai)
  ggplot(Glu_wide_PCA) + aes(y=dim1, x=`L Anterior Insula`) +
    geom_point() + geom_smooth(method="lm") + 
    theme_classic(base_size = 13) +xlab("L Anterior Insula") + ylab("dim1") + 
    theme(legend.key = element_rect(fill = "white", colour = "black"))
  
  
  # GABA
  GABA_wide_PCA <- pivot_wider(MRS_GABA, id_cols=ld8, names_from=label, values_from=GABA.Cr)
  GABA_wide_PCA <- dplyr::select(GABA_wide_PCA, -`R STS`, -`L STS`, -`R Thalamus`, -`R Caudate`, -`L Caudate`, -`R Posterior Insula`, -`L Posterior Insula`)
  
  nb_GABA_kfold <- estim_ncpPCA(GABA_wide_PCA %>% dplyr::select(-ld8), method.cv = "Kfold", verbose = FALSE) 
  nb_GABA_kfold$ncp
  nb_GABA_kfold$ncp <- 2
  
  plot(0:5, nb_GABA_kfold$criterion, xlab = "nb dim", ylab = "MSEP")
  
  res_comp_GABA <- imputePCA(GABA_wide_PCA %>% dplyr::select(-ld8) %>% as.data.frame, ncp = nb_GABA_kfold$ncp) # iterativePCA algorithm
  
  GABA_pca <- PCA(res_comp_GABA, ncp = nb_GABA_kfold$ncp)
  plot(GABA_pca)
  plot(GABA_pca, choix="var")
  GABA_wide_PCA$dim1 <- GABA_pca$ind$coord[,1]
  dim1_lai <- lm(data=GABA_wide_PCA, `MPFC` ~ dim1)
  summary(dim1_lai)
  ggplot(GABA_wide_PCA) + aes(y=dim1, x=`L Anterior Insula`) +
    geom_point() + geom_smooth(method="lm") + 
    theme_classic(base_size = 13) +xlab("L Anterior Insula") + ylab("dim1") + 
    theme(legend.key = element_rect(fill = "white", colour = "black"))
  

 PCA_measures <- merge(Glu_wide_PCA, GABA_wide_PCA, by = "ld8", suffixes = c("_Glu", "_GABA"))
 PCA_measures <- PCA_measures %>% separate(ld8, c("ID", "Date"), sep = "_" )

 EEG_PCA <- merge(alldata_delayOnly, PCA_measures, by = "ID")  
 
 
 # Gamma number of events vs GABA PCA
 EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_PCA %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = Gamma.Event_Number, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_GABA ~ Gamma.Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
 summary(lm.model)
 
 # Gamma number of events variability vs GABA PCA
 EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_PCA %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA, aes(x = Gamma.Event_Number_Variability, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_GABA ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA)
 summary(lm.model)
 
 # Gamma power vs GABA PCA
 EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = Gamma.log1p_Trial_Power, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_GABA ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
 summary(lm.model)
 
 # Gamma power variability vs GABA PCA
 EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA, aes(x = Gamma.log1p_Trial_Power_Variability, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_GABA ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA)
 summary(lm.model)
 
 
 #glutamate
 # Gamma number of events vs Glu PCA
 EEG_MRS_data_DLPFC_cleanEventGlu <- EEG_PCA %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanEventGlu, aes(x = Gamma.Event_Number, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_Glu ~ Gamma.Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGlu)
 summary(lm.model)
 
 # Gamma number of events variability vs Glu PCA
 EEG_MRS_data_DLPFC_cleanEventVarGlu <- EEG_PCA %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGlu, aes(x = Gamma.Event_Number_Variability, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_Glu ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGlu)
 summary(lm.model)
 
 # Gamma power vs Glu PCA
 EEG_MRS_data_DLPFC_cleanPowerGlu <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanEventGlu, aes(x = Gamma.log1p_Trial_Power, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_Glu ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanEventGlu)
 summary(lm.model)
 
 # Gamma power variability vs Glu PCA
 EEG_MRS_data_DLPFC_cleanPowerVarGlu <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
 
 ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGlu, aes(x = Gamma.log1p_Trial_Power_Variability, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
 
 lm.model <- lm(dim1_Glu ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGlu)
 summary(lm.model)
 
}



DLPFC_EEG_vs_DLPFC_MRS <- function () {

  # Prep EEG Data 
  
  GammaDelay_Sublevel_Channel <- DelayOnly_IndividualChannels_SubLevel()
  
  Gamma_DLPFC <- GammaDelay_Sublevel_Channel %>% filter(Label == "'F3'")
  Gamma_DRPFC <- GammaDelay_Sublevel_Channel %>% filter(Label == "'F4'")
  
  Gamma_DLPFC_all <- merge(Gamma_DLPFC, Gamma_DRPFC, by = c("Subject", "Trial", "age"), suffixes = c("_Left","_Right"))
  
  #create avg gamma number between the left and right DLPFC
  Gamma_DLPFC_all$avg_Gamma_Event_Number <- (Gamma_DLPFC_all$Gamma_Event_Number_Left + Gamma_DLPFC_all$Gamma_Event_Number_Right)/2
  
  Gamma_DLPFC_all$avg_Gamma_Event_Number_Variability <- (Gamma_DLPFC_all$Gamma_Event_Number_Variability_Left + Gamma_DLPFC_all$Gamma_Event_Number_Variability_Right)/2
  
  #create avg gamma power between the left and right DLPFC
  Gamma_DLPFC_all$avg_log_Gamma_Trial_Power <- (Gamma_DLPFC_all$log_Gamma_Power_Left + Gamma_DLPFC_all$log_Gamma_Power_Right)/2
 
  Gamma_DLPFC_all$avg_log_Gamma_Trial_Power_Variability <- (Gamma_DLPFC_all$log_Gamma_Power_Variability_Left + Gamma_DLPFC_all$log_Gamma_Power_Variability_Right)/2
   
  Gamma_DLPFC_all_avgEventSub <- aggregate(avg_Gamma_Event_Number ~ Subject, Gamma_DLPFC_all, mean)
  Gamma_DLPFC_all_avgEventVarSub <- aggregate(avg_Gamma_Event_Number_Variability ~ Subject, Gamma_DLPFC_all, mean)
  
  Gamma_DLPFC_all_avgPowerSub <- aggregate(avg_log_Gamma_Trial_Power ~ Subject, Gamma_DLPFC_all, mean)
  Gamma_DLPFC_all_avgPowerVarSub <- aggregate(avg_log_Gamma_Trial_Power_Variability ~ Subject, Gamma_DLPFC_all, mean)
  
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_all_avgEventSub, Gamma_DLPFC_all_avgPowerSub, by = "Subject")
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_Gamma_SubLevel, Gamma_DLPFC_all_avgEventVarSub, by = "Subject")
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_Gamma_SubLevel, Gamma_DLPFC_all_avgPowerVarSub, by = "Subject")
  
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_Gamma_SubLevel, agefile, by = "Subject")
  
  Gamma_DLPFC_Gamma_SubLevel <- Gamma_DLPFC_Gamma_SubLevel %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  # Prep MRS GABA data 
  DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
  DLPFC_GABA <- DLPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_GABA <- DLPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  DLPFC_GABA$Subject <- DLPFC_GABA$ld8
  DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_GABA, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  # Gamma power of events vs GABA
  EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGABA, aes(x = avg_log_Gamma_Trial_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA
  EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA, aes(x = avg_log_Gamma_Trial_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs GABA
  EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = avg_Gamma_Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs GABA
  EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA, aes(x = avg_Gamma_Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA)
  summary(lm.model)
  
  
  # DLPFC Glutamate 
  
  DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
  DLPFC_Glu <- DLPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_Glu <- DLPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  DLPFC_Glu$Subject <- DLPFC_Glu$ld8
  DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_Glu, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  
  # Gamma number of events vs Glu
  EEG_MRS_data_DLPFC_cleanEventGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGlu, aes(x = avg_Gamma_Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGlu)
  summary(lm.model)
  
  # Gamma number of events variability vs Glu
  EEG_MRS_data_DLPFC_cleanEventVarGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGlu, aes(x = avg_Gamma_Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGlu)
  summary(lm.model)
  
  
  # Gamma power vs Glu
  EEG_MRS_data_DLPFC_cleanPowerGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGlu, aes(x = avg_log_Gamma_Trial_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGlu)
  summary(lm.model)
  
  
  # Gamma power variability vs Glu
  EEG_MRS_data_DLPFC_cleanPowerVarGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGlu, aes(x = avg_log_Gamma_Trial_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGlu)
  summary(lm.model)
  
  
  # Looking at glu/gaba ratio 
  
  # Prep MRS GABA data 
  DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
  DLPFC_GABA <- DLPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_GABA <- DLPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  DLPFC_GABA$Subject <- DLPFC_GABA$ld8
  DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_GABA, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  # DLPFC Glutamate 
  
  DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
  DLPFC_Glu <- DLPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_Glu <- DLPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  DLPFC_Glu$Subject <- DLPFC_Glu$ld8
  DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)
  
  EEG_MRS_data_DLPFC <- merge(EEG_MRS_data_DLPFC, Avg_DLPFC_Glu, by = "ID")
  
  EEG_MRS_data_DLPFC$gluGABAratio <- EEG_MRS_data_DLPFC$Glu.Cr / EEG_MRS_data_DLPFC$GABA.Cr

  # Gamma power of events vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGABA, aes(x = avg_log_Gamma_Trial_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA, aes(x = avg_log_Gamma_Trial_Power_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = avg_Gamma_Event_Number, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA, aes(x = avg_Gamma_Event_Number_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA)
  summary(lm.model)
  
  
  # DLPFC EEG vs. Gaba Glu Residuals 
  # Prep MRS GABA data 
  DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
  DLPFC_GABA <- DLPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_GABA <- DLPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  DLPFC_GABA$Subject <- DLPFC_GABA$ld8
  DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_GABA, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  # DLPFC Glutamate 
  
  DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
  DLPFC_Glu <- DLPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_Glu <- DLPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  DLPFC_Glu$Subject <- DLPFC_Glu$ld8
  DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)
  
  EEG_MRS_data_DLPFC <- merge(EEG_MRS_data_DLPFC, Avg_DLPFC_Glu, by = "ID")
  
 # extract residuals 
  model <- lm(data = EEG_MRS_data_DLPFC, GABA.Cr ~ Glu.Cr * age) 
  EEG_MRS_data_DLPFC$residuals <- residuals(model) 
  
  # Gamma power of events vs GABA Glu residuals 
  EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGABA, aes(x = avg_log_Gamma_Trial_Power, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA Glu residuals
  EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA, aes(x = avg_log_Gamma_Trial_Power_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs GABA Glu residuals
  EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = avg_Gamma_Event_Number, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs GABA Glu residuals
  EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA, aes(x = avg_Gamma_Event_Number_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA)
  summary(lm.model)
  
}

MPFC_EEG_vs_MPFC_MRS <- function () {
  
  # Prep EEG Data 
  
  GammaDelay_Sublevel_Channel <- DelayOnly_IndividualChannels_SubLevel()
  
  Gamma_MPFC_all <- filter(GammaDelay_Sublevel_Channel, str_detect(GammaDelay_Sublevel_Channel$Label, "F"))
  
  #create avg gamma number across MPFC
  avg_Gamma_Event_Number <- aggregate(Gamma_Event_Number ~ Subject, Gamma_MPFC_all, mean)
  avg_Gamma_Event_Number_Variability <- aggregate(Gamma_Event_Number_Variability ~ Subject, Gamma_MPFC_all, mean)
  
  #create avg gamma power between the left and right MPFC
  avg_log_Gamma_Trial_Power <- aggregate(log_Gamma_Power ~ Subject, Gamma_MPFC_all, mean)
  avg_log_Gamma_Trial_Power_Variability <- aggregate(log_Gamma_Power_Variability ~ Subject, Gamma_MPFC_all, mean)
  
  Gamma_MPFC_Gamma_SubLevel <- merge(avg_Gamma_Event_Number, avg_Gamma_Event_Number_Variability, by = "Subject")
  Gamma_MPFC_Gamma_SubLevel <- merge(Gamma_MPFC_Gamma_SubLevel, avg_log_Gamma_Trial_Power, by = "Subject")
  Gamma_MPFC_Gamma_SubLevel <- merge(Gamma_MPFC_Gamma_SubLevel, avg_log_Gamma_Trial_Power_Variability, by = "Subject")
  
  Gamma_MPFC_Gamma_SubLevel <- merge(Gamma_MPFC_Gamma_SubLevel, agefile, by = "Subject")
  
  Gamma_MPFC_Gamma_SubLevel <- Gamma_MPFC_Gamma_SubLevel %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  # Prep MRS GABA data 
  MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
  MPFC_GABA <- MPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_GABA <- MPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  MPFC_GABA$Subject <- MPFC_GABA$ld8
  MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_GABA, by = "ID")
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  # Gamma power of events vs GABA
  EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGABA, aes(x = log_Gamma_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA
  EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA, aes(x = log_Gamma_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs GABA
  EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma_Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs GABA
  EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA, aes(x = Gamma_Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA)
  summary(lm.model)
  
  
  # MPFC Glutamate 
  
  MPFC_Glu <- MRS_glu %>% filter(roi == 8)
  MPFC_Glu <- MPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_Glu <- MPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  MPFC_Glu$Subject <- MPFC_Glu$ld8
  MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_Glu, by = "ID")
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  
  # Gamma number of events vs Glu
  EEG_MRS_data_MPFC_cleanEventGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGlu, aes(x = Gamma_Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGlu)
  summary(lm.model)
  
  # Gamma number of events variability vs Glu
  EEG_MRS_data_MPFC_cleanEventVarGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGlu, aes(x = Gamma_Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGlu)
  summary(lm.model)
  
  
  # Gamma power vs Glu
  EEG_MRS_data_MPFC_cleanPowerGlu <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGlu, aes(x = log_Gamma_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGlu)
  summary(lm.model)
  
  
  # Gamma power variability vs Glu
  EEG_MRS_data_MPFC_cleanPowerVarGlu <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGlu, aes(x = log_Gamma_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGlu)
  summary(lm.model)
  
  
  
  
  # looking at Glu/GABA
  
  # MPFC Glutamate 
  MPFC_Glu <- MRS_glu %>% filter(roi == 8)
  MPFC_Glu <- MPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_Glu <- MPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  MPFC_Glu$Subject <- MPFC_Glu$ld8
  MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)
  
  
  # Prep MRS GABA data 
  MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
  MPFC_GABA <- MPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_GABA <- MPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  MPFC_GABA$Subject <- MPFC_GABA$ld8
  MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_GABA, by = "ID")
  EEG_MRS_data_MPFC <- merge(EEG_MRS_data_MPFC, Avg_MPFC_Glu, by = "ID")
  
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  EEG_MRS_data_MPFC$gluGABAratio <- EEG_MRS_data_MPFC$Glu.Cr / EEG_MRS_data_MPFC$GABA.Cr

  
  # Gamma power of events vs Glu/GABA
  EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGABA, aes(x = log_Gamma_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGABA)
  summary(lm.model)
  
  # Gamma power of events variability vs Glu/GABA
  EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA, aes(x = log_Gamma_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs Glu/GABA
  EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma_Event_Number, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs Glu/GABA
  EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA, aes(x = Gamma_Event_Number_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA)
  summary(lm.model)
  
  
  # MPFC EEG vs. Gaba Glu Residuals 
  # Prep MRS GABA data 
  
  
  # MPFC Glutamate 
  MPFC_Glu <- MRS_glu %>% filter(roi == 8)
  MPFC_Glu <- MPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_Glu <- MPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  MPFC_Glu$Subject <- MPFC_Glu$ld8
  MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)
  
  
  # Prep MRS GABA data 
  MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
  MPFC_GABA <- MPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_GABA <- MPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  MPFC_GABA$Subject <- MPFC_GABA$ld8
  MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_GABA, by = "ID")
  EEG_MRS_data_MPFC <- merge(EEG_MRS_data_MPFC, Avg_MPFC_Glu, by = "ID")
  
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  
  
  # extract residuals 
  model <- lm(data = EEG_MRS_data_MPFC, GABA.Cr ~ Glu.Cr * age) 
  EEG_MRS_data_MPFC$residuals <- residuals(model) 
  
  # Gamma power of events vs GABA Glu residuals 
  EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGABA, aes(x = log_Gamma_Power, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA Glu residuals
  EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA, aes(x = log_Gamma_Power_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs GABA Glu residuals
  EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma_Event_Number, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs GABA Glu residuals
  EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA, aes(x = Gamma_Event_Number_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA)
  summary(lm.model)
  
  
}

ACC_EEG_vs_ACC_MRS <- function () {
  
  # Prep EEG Data 
  
  GammaDelay_Sublevel_Channel <- DelayOnly_IndividualChannels_SubLevel()
  
  Gamma_ACC_all <- filter(GammaDelay_Sublevel_Channel, str_detect(GammaDelay_Sublevel_Channel$Label, "F"))
  
  #create avg gamma number across ACC
  avg_Gamma_Event_Number <- aggregate(Gamma_Event_Number ~ Subject, Gamma_ACC_all, mean)
  avg_Gamma_Event_Number_Variability <- aggregate(Gamma_Event_Number_Variability ~ Subject, Gamma_ACC_all, mean)
  
  #create avg gamma power between the left and right ACC
  avg_log_Gamma_Trial_Power <- aggregate(log_Gamma_Power ~ Subject, Gamma_ACC_all, mean)
  avg_log_Gamma_Trial_Power_Variability <- aggregate(log_Gamma_Power_Variability ~ Subject, Gamma_ACC_all, mean)
  
  Gamma_ACC_Gamma_SubLevel <- merge(avg_Gamma_Event_Number, avg_Gamma_Event_Number_Variability, by = "Subject")
  Gamma_ACC_Gamma_SubLevel <- merge(Gamma_ACC_Gamma_SubLevel, avg_log_Gamma_Trial_Power, by = "Subject")
  Gamma_ACC_Gamma_SubLevel <- merge(Gamma_ACC_Gamma_SubLevel, avg_log_Gamma_Trial_Power_Variability, by = "Subject")
  
  Gamma_ACC_Gamma_SubLevel <- merge(Gamma_ACC_Gamma_SubLevel, agefile, by = "Subject")
  
  Gamma_ACC_Gamma_SubLevel <- Gamma_ACC_Gamma_SubLevel %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  # Prep MRS GABA data 
  ACC_GABA <- MRS_GABA %>% filter(roi == 8)
  ACC_GABA <- ACC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  ACC_GABA <- ACC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  ACC_GABA$Subject <- ACC_GABA$ld8
  ACC_GABA <- ACC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_ACC_GABA <- aggregate(GABA.Cr ~ ID, data = ACC_GABA, mean)
  
  EEG_MRS_data_ACC <- merge(Gamma_ACC_Gamma_SubLevel, Avg_ACC_GABA, by = "ID")
  EEG_MRS_data_ACC$Group <- factor(EEG_MRS_data_ACC$Group)
  
  # Gamma power of events vs GABA
  EEG_MRS_data_ACC_cleanPowerGABA <- EEG_MRS_data_ACC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanPowerGABA, aes(x = log_Gamma_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ log_Gamma_Power + age, data = EEG_MRS_data_ACC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA
  EEG_MRS_data_ACC_cleanPowerVarGABA <- EEG_MRS_data_ACC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanPowerVarGABA, aes(x = log_Gamma_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_ACC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs GABA
  EEG_MRS_data_ACC_cleanEventGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanEventGABA, aes(x = Gamma_Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ Gamma_Event_Number + age, data = EEG_MRS_data_ACC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs GABA
  EEG_MRS_data_ACC_cleanEventVarGABA <- EEG_MRS_data_ACC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanEventVarGABA, aes(x = Gamma_Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_ACC_cleanEventVarGABA)
  summary(lm.model)
  
  
  # ACC Glutamate 
  
  ACC_Glu <- MRS_glu %>% filter(roi == 8)
  ACC_Glu <- ACC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  ACC_Glu <- ACC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  ACC_Glu$Subject <- ACC_Glu$ld8
  ACC_Glu <- ACC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_ACC_Glu <- aggregate(Glu.Cr ~ ID, data = ACC_Glu, mean)
  
  EEG_MRS_data_ACC <- merge(Gamma_ACC_Gamma_SubLevel, Avg_ACC_Glu, by = "ID")
  EEG_MRS_data_ACC$Group <- factor(EEG_MRS_data_ACC$Group)
  
  
  # Gamma number of events vs Glu
  EEG_MRS_data_ACC_cleanEventGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanEventGlu, aes(x = Gamma_Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ Gamma_Event_Number + age, data = EEG_MRS_data_ACC_cleanEventGlu)
  summary(lm.model)
  
  # Gamma number of events variability vs Glu
  EEG_MRS_data_ACC_cleanEventVarGlu <- EEG_MRS_data_ACC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanEventVarGlu, aes(x = Gamma_Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_ACC_cleanEventVarGlu)
  summary(lm.model)
  
  
  # Gamma power vs Glu
  EEG_MRS_data_ACC_cleanPowerGlu <- EEG_MRS_data_ACC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanPowerGlu, aes(x = log_Gamma_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ log_Gamma_Power + age, data = EEG_MRS_data_ACC_cleanPowerGlu)
  summary(lm.model)
  
  
  # Gamma power variability vs Glu
  EEG_MRS_data_ACC_cleanPowerVarGlu <- EEG_MRS_data_ACC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_ACC_cleanPowerVarGlu, aes(x = log_Gamma_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_ACC_cleanPowerVarGlu)
  summary(lm.model)
  
  
  
}



DLPFC_EEG_vs_DLPFC_MRS_AgeGroups <- function () {
  
  # Prep EEG Data 
  
  GammaDelay_Sublevel_Channel <- DelayOnly_IndividualChannels_SubLevel()
  
  Gamma_DLPFC <- GammaDelay_Sublevel_Channel %>% filter(Label == "'F3'")
  Gamma_DRPFC <- GammaDelay_Sublevel_Channel %>% filter(Label == "'F4'")
  
  Gamma_DLPFC_all <- merge(Gamma_DLPFC, Gamma_DRPFC, by = c("Subject", "Trial", "age"), suffixes = c("_Left","_Right"))
  
  #create avg gamma number between the left and right DLPFC
  Gamma_DLPFC_all$avg_Gamma_Event_Number <- (Gamma_DLPFC_all$Gamma_Event_Number_Left + Gamma_DLPFC_all$Gamma_Event_Number_Right)/2
  
  Gamma_DLPFC_all$avg_Gamma_Event_Number_Variability <- (Gamma_DLPFC_all$Gamma_Event_Number_Variability_Left + Gamma_DLPFC_all$Gamma_Event_Number_Variability_Right)/2
  
  #create avg gamma power between the left and right DLPFC
  Gamma_DLPFC_all$avg_log_Gamma_Trial_Power <- (Gamma_DLPFC_all$log_Gamma_Power_Left + Gamma_DLPFC_all$log_Gamma_Power_Right)/2
  
  Gamma_DLPFC_all$avg_log_Gamma_Trial_Power_Variability <- (Gamma_DLPFC_all$log_Gamma_Power_Variability_Left + Gamma_DLPFC_all$log_Gamma_Power_Variability_Right)/2
  
  Gamma_DLPFC_all_avgEventSub <- aggregate(avg_Gamma_Event_Number ~ Subject, Gamma_DLPFC_all, mean)
  Gamma_DLPFC_all_avgEventVarSub <- aggregate(avg_Gamma_Event_Number_Variability ~ Subject, Gamma_DLPFC_all, mean)
  
  Gamma_DLPFC_all_avgPowerSub <- aggregate(avg_log_Gamma_Trial_Power ~ Subject, Gamma_DLPFC_all, mean)
  Gamma_DLPFC_all_avgPowerVarSub <- aggregate(avg_log_Gamma_Trial_Power_Variability ~ Subject, Gamma_DLPFC_all, mean)
  
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_all_avgEventSub, Gamma_DLPFC_all_avgPowerSub, by = "Subject")
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_Gamma_SubLevel, Gamma_DLPFC_all_avgEventVarSub, by = "Subject")
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_Gamma_SubLevel, Gamma_DLPFC_all_avgPowerVarSub, by = "Subject")
  
  Gamma_DLPFC_Gamma_SubLevel <- merge(Gamma_DLPFC_Gamma_SubLevel, agefile, by = "Subject")
  
  Gamma_DLPFC_Gamma_SubLevel <- Gamma_DLPFC_Gamma_SubLevel %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  # Prep MRS GABA data 
  DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
  DLPFC_GABA <- DLPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_GABA <- DLPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  DLPFC_GABA$Subject <- DLPFC_GABA$ld8
  DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_GABA, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  # Gamma power of events vs GABA
  EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGABA[EEG_MRS_data_DLPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerGABA$age <= 16,], aes(x = avg_log_Gamma_Trial_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGABA[EEG_MRS_data_DLPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA
  EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA[EEG_MRS_data_DLPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerVarGABA$age <= 16,], aes(x = avg_log_Gamma_Trial_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA[EEG_MRS_data_DLPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerVarGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma number of events vs GABA
  EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA[EEG_MRS_data_DLPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventGABA$age <= 16,], aes(x = avg_Gamma_Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA[EEG_MRS_data_DLPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventGABA$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs GABA
  EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA[EEG_MRS_data_DLPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventVarGABA$age <= 16,], aes(x = avg_Gamma_Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA[EEG_MRS_data_DLPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventVarGABA$age <= 16,])
  summary(lm.model)
  
  
  # DLPFC Glutamate 
  
  DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
  DLPFC_Glu <- DLPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_Glu <- DLPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  DLPFC_Glu$Subject <- DLPFC_Glu$ld8
  DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_Glu, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  
  # Gamma number of events vs Glu
  EEG_MRS_data_DLPFC_cleanEventGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGlu[EEG_MRS_data_DLPFC_cleanEventGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanEventGlu$age <= 16,], aes(x = avg_Gamma_Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGlu[EEG_MRS_data_DLPFC_cleanEventGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanEventGlu$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs Glu
  EEG_MRS_data_DLPFC_cleanEventVarGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGlu[EEG_MRS_data_DLPFC_cleanEventVarGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanEventVarGlu$age <= 16,], aes(x = avg_Gamma_Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGlu[EEG_MRS_data_DLPFC_cleanEventVarGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanEventVarGlu$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power vs Glu
  EEG_MRS_data_DLPFC_cleanPowerGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGlu[EEG_MRS_data_DLPFC_cleanPowerGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerGlu$age <= 16,], aes(x = avg_log_Gamma_Trial_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGlu[EEG_MRS_data_DLPFC_cleanPowerGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerGlu$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power variability vs Glu
  EEG_MRS_data_DLPFC_cleanPowerVarGlu <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGlu[EEG_MRS_data_DLPFC_cleanPowerVarGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerVarGlu$age <= 16,], aes(x = avg_log_Gamma_Trial_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGlu[EEG_MRS_data_DLPFC_cleanPowerVarGlu$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerVarGlu$age <= 16,])
  summary(lm.model)
  
  
  # Looking at glu/gaba ratio 
  
  # Prep MRS GABA data 
  DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
  DLPFC_GABA <- DLPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_GABA <- DLPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  DLPFC_GABA$Subject <- DLPFC_GABA$ld8
  DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_GABA, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  # DLPFC Glutamate 
  
  DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
  DLPFC_Glu <- DLPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_Glu <- DLPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  DLPFC_Glu$Subject <- DLPFC_Glu$ld8
  DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)
  
  EEG_MRS_data_DLPFC <- merge(EEG_MRS_data_DLPFC, Avg_DLPFC_Glu, by = "ID")
  
  EEG_MRS_data_DLPFC$gluGABAratio <- EEG_MRS_data_DLPFC$Glu.Cr / EEG_MRS_data_DLPFC$GABA.Cr
  
  # Gamma power of events vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGABA, aes(x = avg_log_Gamma_Trial_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA, aes(x = avg_log_Gamma_Trial_Power_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA, aes(x = avg_Gamma_Event_Number, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs Glu/GABA
  EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA, aes(x = avg_Gamma_Event_Number_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA)
  summary(lm.model)
  
  
  # DLPFC EEG vs. Gaba Glu Residuals 
  # Prep MRS GABA data 
  DLPFC_GABA <- MRS_GABA %>% filter(roi == c(9,10))
  DLPFC_GABA <- DLPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_GABA <- DLPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  DLPFC_GABA$Subject <- DLPFC_GABA$ld8
  DLPFC_GABA <- DLPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_GABA <- aggregate(GABA.Cr ~ ID, data = DLPFC_GABA, mean)
  
  EEG_MRS_data_DLPFC <- merge(Gamma_DLPFC_Gamma_SubLevel, Avg_DLPFC_GABA, by = "ID")
  EEG_MRS_data_DLPFC$Group <- factor(EEG_MRS_data_DLPFC$Group)
  
  # DLPFC Glutamate 
  
  DLPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
  DLPFC_Glu <- DLPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  DLPFC_Glu <- DLPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  DLPFC_Glu$Subject <- DLPFC_Glu$ld8
  DLPFC_Glu <- DLPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_DLPFC_Glu <- aggregate(Glu.Cr ~ ID, data = DLPFC_Glu, mean)
  
  EEG_MRS_data_DLPFC <- merge(EEG_MRS_data_DLPFC, Avg_DLPFC_Glu, by = "ID")
  
  # extract residuals 
  model <- lm(data = EEG_MRS_data_DLPFC, GABA.Cr ~ Glu.Cr * age) 
  EEG_MRS_data_DLPFC$residuals <- residuals(model) 
  
  # Gamma power of events vs GABA Glu residuals 
  EEG_MRS_data_DLPFC_cleanPowerGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power - mean(avg_log_Gamma_Trial_Power, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerGABA[EEG_MRS_data_DLPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerGABA$age <= 16,], aes(x = avg_log_Gamma_Trial_Power, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_log_Gamma_Trial_Power + age, data = EEG_MRS_data_DLPFC_cleanPowerGABA[EEG_MRS_data_DLPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA Glu residuals
  EEG_MRS_data_DLPFC_cleanPowerVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_log_Gamma_Trial_Power_Variability - mean(avg_log_Gamma_Trial_Power_Variability, na.rm= T)) < (sd(avg_log_Gamma_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanPowerVarGABA[EEG_MRS_data_DLPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerVarGABA$age <= 16,], aes(x = avg_log_Gamma_Trial_Power_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_log_Gamma_Trial_Power_Variability + age, data = EEG_MRS_data_DLPFC_cleanPowerVarGABA[EEG_MRS_data_DLPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanPowerVarGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma number of events vs GABA Glu residuals
  EEG_MRS_data_DLPFC_cleanEventGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number - mean(avg_Gamma_Event_Number, na.rm= T)) < (sd(avg_Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventGABA[EEG_MRS_data_DLPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventGABA$age <= 16,], aes(x = avg_Gamma_Event_Number, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_Gamma_Event_Number + age, data = EEG_MRS_data_DLPFC_cleanEventGABA[EEG_MRS_data_DLPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventGABA$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs GABA Glu residuals
  EEG_MRS_data_DLPFC_cleanEventVarGABA <- EEG_MRS_data_DLPFC %>% filter(abs(avg_Gamma_Event_Number_Variability - mean(avg_Gamma_Event_Number_Variability, na.rm= T)) < (sd(avg_Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_DLPFC_cleanEventVarGABA[EEG_MRS_data_DLPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventVarGABA$age <= 16,], aes(x = avg_Gamma_Event_Number_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ avg_Gamma_Event_Number_Variability + age, data = EEG_MRS_data_DLPFC_cleanEventVarGABA[EEG_MRS_data_DLPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_DLPFC_cleanEventVarGABA$age <= 16,])
  summary(lm.model)
  
}

MPFC_EEG_vs_MPFC_MRS_AgeGroups <- function () {
  
  # Prep EEG Data 
  
  GammaDelay_Sublevel_Channel <- DelayOnly_IndividualChannels_SubLevel()
  
  Gamma_MPFC_all <- filter(GammaDelay_Sublevel_Channel, str_detect(GammaDelay_Sublevel_Channel$Label, "F"))
  
  #create avg gamma number across MPFC
  avg_Gamma_Event_Number <- aggregate(Gamma_Event_Number ~ Subject, Gamma_MPFC_all, mean)
  avg_Gamma_Event_Number_Variability <- aggregate(Gamma_Event_Number_Variability ~ Subject, Gamma_MPFC_all, mean)
  
  #create avg gamma power between the left and right MPFC
  avg_log_Gamma_Trial_Power <- aggregate(log_Gamma_Power ~ Subject, Gamma_MPFC_all, mean)
  avg_log_Gamma_Trial_Power_Variability <- aggregate(log_Gamma_Power_Variability ~ Subject, Gamma_MPFC_all, mean)
  
  Gamma_MPFC_Gamma_SubLevel <- merge(avg_Gamma_Event_Number, avg_Gamma_Event_Number_Variability, by = "Subject")
  Gamma_MPFC_Gamma_SubLevel <- merge(Gamma_MPFC_Gamma_SubLevel, avg_log_Gamma_Trial_Power, by = "Subject")
  Gamma_MPFC_Gamma_SubLevel <- merge(Gamma_MPFC_Gamma_SubLevel, avg_log_Gamma_Trial_Power_Variability, by = "Subject")
  
  Gamma_MPFC_Gamma_SubLevel <- merge(Gamma_MPFC_Gamma_SubLevel, agefile, by = "Subject")
  
  Gamma_MPFC_Gamma_SubLevel <- Gamma_MPFC_Gamma_SubLevel %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  # Prep MRS GABA data 
  MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
  MPFC_GABA <- MPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_GABA <- MPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  MPFC_GABA$Subject <- MPFC_GABA$ld8
  MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_GABA, by = "ID")
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  # Gamma power of events vs GABA
  EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGABA[EEG_MRS_data_MPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerGABA$age <= 16,], aes(x = log_Gamma_Power, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGABA[EEG_MRS_data_MPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA
  EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA[EEG_MRS_data_MPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerVarGABA$age <= 16,], aes(x = log_Gamma_Power_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA[EEG_MRS_data_MPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerVarGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma number of events vs GABA
  EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA[EEG_MRS_data_MPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventGABA$age <= 16,], aes(x = Gamma_Event_Number, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA[EEG_MRS_data_MPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventGABA$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs GABA
  EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(GABA.Cr - mean(GABA.Cr, na.rm= T)) < (sd(GABA.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA[EEG_MRS_data_MPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventVarGABA$age <= 16,], aes(x = Gamma_Event_Number_Variability, y = GABA.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(GABA.Cr ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA[EEG_MRS_data_MPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventVarGABA$age <= 16,])
  summary(lm.model)
  
  
  # MPFC Glutamate 
  
  MPFC_Glu <- MRS_glu %>% filter(roi == 8)
  MPFC_Glu <- MPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_Glu <- MPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  MPFC_Glu$Subject <- MPFC_Glu$ld8
  MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_Glu, by = "ID")
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  
  # Gamma number of events vs Glu
  EEG_MRS_data_MPFC_cleanEventGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGlu[EEG_MRS_data_MPFC_cleanEventGlu$age >= 10 & EEG_MRS_data_MPFC_cleanEventGlu$age <= 16,], aes(x = Gamma_Event_Number, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGlu[EEG_MRS_data_MPFC_cleanEventGlu$age >= 10 & EEG_MRS_data_MPFC_cleanEventGlu$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs Glu
  EEG_MRS_data_MPFC_cleanEventVarGlu <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGlu[EEG_MRS_data_MPFC_cleanEventVarGlu$age >= 10 & EEG_MRS_data_MPFC_cleanEventVarGlu$age <= 16,], aes(x = Gamma_Event_Number_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGlu[EEG_MRS_data_MPFC_cleanEventVarGlu$age >= 10 & EEG_MRS_data_MPFC_cleanEventVarGlu$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power vs Glu
  EEG_MRS_data_MPFC_cleanPowerGlu <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGlu[EEG_MRS_data_MPFC_cleanPowerGlu$age >= 10 & EEG_MRS_data_MPFC_cleanPowerGlu$age <= 16,], aes(x = log_Gamma_Power, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGlu[EEG_MRS_data_MPFC_cleanPowerGlu$age >= 10 & EEG_MRS_data_MPFC_cleanPowerGlu$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power variability vs Glu
  EEG_MRS_data_MPFC_cleanPowerVarGlu <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(Glu.Cr - mean(Glu.Cr, na.rm= T)) < (sd(Glu.Cr, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGlu[EEG_MRS_data_MPFC_cleanPowerVarGlu$age >= 10 & EEG_MRS_data_MPFC_cleanPowerVarGlu$age <= 16,], aes(x = log_Gamma_Power_Variability, y = Glu.Cr)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(Glu.Cr ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGlu[EEG_MRS_data_MPFC_cleanPowerVarGlu$age >= 10 & EEG_MRS_data_MPFC_cleanPowerVarGlu$age <= 16,])
  summary(lm.model)
  
  
  # Looking at glu/gaba ratio 
  
  # Prep MRS GABA data 
  MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
  MPFC_GABA <- MPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_GABA <- MPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  MPFC_GABA$Subject <- MPFC_GABA$ld8
  MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_GABA, by = "ID")
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  # MPFC Glutamate 
  
  MPFC_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
  MPFC_Glu <- MPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_Glu <- MPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  MPFC_Glu$Subject <- MPFC_Glu$ld8
  MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)
  
  EEG_MRS_data_MPFC <- merge(EEG_MRS_data_MPFC, Avg_MPFC_Glu, by = "ID")
  
  EEG_MRS_data_MPFC$gluGABAratio <- EEG_MRS_data_MPFC$Glu.Cr / EEG_MRS_data_MPFC$GABA.Cr
  
  # Gamma power of events vs Glu/GABA
  EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGABA, aes(x = log_Gamma_Power, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGABA)
  summary(lm.model)
  
  
  # Gamma power of events variability vs Glu/GABA
  EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA, aes(x = log_Gamma_Power_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA)
  summary(lm.model)
  
  
  # Gamma number of events vs Glu/GABA
  EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA, aes(x = Gamma_Event_Number, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA)
  summary(lm.model)
  
  # Gamma number of events variability vs Glu/GABA
  EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(gluGABAratio - mean(gluGABAratio, na.rm= T)) < (sd(gluGABAratio, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA, aes(x = Gamma_Event_Number_Variability, y = gluGABAratio)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(gluGABAratio ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA)
  summary(lm.model)
  
  
  # MPFC EEG vs. Gaba Glu Residuals 
  # Prep MRS GABA data 
  MPFC_GABA <- MRS_GABA %>% filter(roi == 8)
  MPFC_GABA <- MPFC_GABA%>%
    mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_GABA <- MPFC_GABA %>% 
    filter(abs(zscore_GABA) <= z_thres)
  MPFC_GABA$Subject <- MPFC_GABA$ld8
  MPFC_GABA <- MPFC_GABA %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_GABA <- aggregate(GABA.Cr ~ ID, data = MPFC_GABA, mean)
  
  EEG_MRS_data_MPFC <- merge(Gamma_MPFC_Gamma_SubLevel, Avg_MPFC_GABA, by = "ID")
  EEG_MRS_data_MPFC$Group <- factor(EEG_MRS_data_MPFC$Group)
  
  # MPFC Glutamate 
  
  MPFC_Glu <- MRS_glu %>% filter(roi == 8)
  MPFC_Glu <- MPFC_Glu%>%
    mutate(zscore_Glu = scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
  MPFC_Glu <- MPFC_Glu %>% 
    filter(abs(zscore_Glu) <= z_thres)
  MPFC_Glu$Subject <- MPFC_Glu$ld8
  MPFC_Glu <- MPFC_Glu %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  Avg_MPFC_Glu <- aggregate(Glu.Cr ~ ID, data = MPFC_Glu, mean)
  
  EEG_MRS_data_MPFC <- merge(EEG_MRS_data_MPFC, Avg_MPFC_Glu, by = "ID")
  
  # extract residuals 
  model <- lm(data = EEG_MRS_data_MPFC, GABA.Cr ~ Glu.Cr * age) 
  EEG_MRS_data_MPFC$residuals <- residuals(model) 
  
  # Gamma power of events vs GABA Glu residuals 
  EEG_MRS_data_MPFC_cleanPowerGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power - mean(log_Gamma_Power, na.rm= T)) < (sd(log_Gamma_Power, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerGABA[EEG_MRS_data_MPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerGABA$age <= 16,], aes(x = log_Gamma_Power, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ log_Gamma_Power + age, data = EEG_MRS_data_MPFC_cleanPowerGABA[EEG_MRS_data_MPFC_cleanPowerGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma power of events variability vs GABA Glu residuals
  EEG_MRS_data_MPFC_cleanPowerVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(log_Gamma_Power_Variability - mean(log_Gamma_Power_Variability, na.rm= T)) < (sd(log_Gamma_Power_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanPowerVarGABA[EEG_MRS_data_MPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerVarGABA$age <= 16,], aes(x = log_Gamma_Power_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ log_Gamma_Power_Variability + age, data = EEG_MRS_data_MPFC_cleanPowerVarGABA[EEG_MRS_data_MPFC_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanPowerVarGABA$age <= 16,])
  summary(lm.model)
  
  
  # Gamma number of events vs GABA Glu residuals
  EEG_MRS_data_MPFC_cleanEventGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number - mean(Gamma_Event_Number, na.rm= T)) < (sd(Gamma_Event_Number, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventGABA[EEG_MRS_data_MPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventGABA$age <= 16,], aes(x = Gamma_Event_Number, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ Gamma_Event_Number + age, data = EEG_MRS_data_MPFC_cleanEventGABA[EEG_MRS_data_MPFC_cleanEventGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventGABA$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs GABA Glu residuals
  EEG_MRS_data_MPFC_cleanEventVarGABA <- EEG_MRS_data_MPFC %>% filter(abs(Gamma_Event_Number_Variability - mean(Gamma_Event_Number_Variability, na.rm= T)) < (sd(Gamma_Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(residuals - mean(residuals, na.rm= T)) < (sd(residuals, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_MPFC_cleanEventVarGABA[EEG_MRS_data_MPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventVarGABA$age <= 16,], aes(x = Gamma_Event_Number_Variability, y = residuals)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(residuals ~ Gamma_Event_Number_Variability + age, data = EEG_MRS_data_MPFC_cleanEventVarGABA[EEG_MRS_data_MPFC_cleanEventVarGABA$age >= 10 & EEG_MRS_data_MPFC_cleanEventVarGABA$age <= 16,])
  summary(lm.model)
  
}


WholeBrainEEG_vs_PCA_MRS <- function () {
  
  all_df <- DelayOnly_Sublevel()
  alldata_delayOnly <- all_df$alldata_delayOnly
  alldata_delayOnly$Subject <- as.character(alldata_delayOnly$Subject)
  alldata_delayOnly <- alldata_delayOnly %>% separate(Subject, c("ID", "Date"), sep = "_" )
  
  
  library(missMDA) 
  library(FactoMineR)
  #### PCA ####
  # Glutamate 
  Glu_wide_PCA <- pivot_wider(MRS_glu, id_cols=ld8, names_from=label, values_from=Glu.Cr)
  Glu_wide_PCA <- dplyr::select(Glu_wide_PCA, -`R STS`, -`L STS`, -`R Thalamus`, -`R Caudate`, -`L Caudate`, -`R Posterior Insula`, -`L Posterior Insula`)
  nb_glu_kfold <- estim_ncpPCA(Glu_wide_PCA %>% dplyr::select(-ld8),method.cv = "Kfold", verbose = FALSE) 
  nb_glu_kfold$ncp
  nb_glu_kfold$ncp <- 2
  
  plot(0:5, nb_glu_kfold$criterion, xlab = "nb dim", ylab = "MSEP")
  
  res_comp_glu <- imputePCA(Glu_wide_PCA %>% dplyr::select(-ld8) %>% as.data.frame, ncp = nb_glu_kfold$ncp) # iterativePCA algorithm
  
  glu_pca <- PCA(res_comp_glu, ncp = nb_glu_kfold$ncp)
  plot(glu_pca)
  plot(glu_pca, choix="var")
  Glu_wide_PCA$dim1 <- glu_pca$ind$coord[,1]
  dim1_lai <- lm(data=Glu_wide_PCA, `MPFC` ~ dim1)
  summary(dim1_lai)
  ggplot(Glu_wide_PCA) + aes(y=dim1, x=`L Anterior Insula`) +
    geom_point() + geom_smooth(method="lm") + 
    theme_classic(base_size = 13) +xlab("L Anterior Insula") + ylab("dim1") + 
    theme(legend.key = element_rect(fill = "white", colour = "black"))
  
  
  # GABA
  GABA_wide_PCA <- pivot_wider(MRS_GABA, id_cols=ld8, names_from=label, values_from=GABA.Cr)
  GABA_wide_PCA <- dplyr::select(GABA_wide_PCA, -`R STS`, -`L STS`, -`R Thalamus`, -`R Caudate`, -`L Caudate`, -`R Posterior Insula`, -`L Posterior Insula`)
  
  nb_GABA_kfold <- estim_ncpPCA(GABA_wide_PCA %>% dplyr::select(-ld8), method.cv = "Kfold", verbose = FALSE) 
  nb_GABA_kfold$ncp
  nb_GABA_kfold$ncp <- 2
  
  plot(0:5, nb_GABA_kfold$criterion, xlab = "nb dim", ylab = "MSEP")
  
  res_comp_GABA <- imputePCA(GABA_wide_PCA %>% dplyr::select(-ld8) %>% as.data.frame, ncp = nb_GABA_kfold$ncp) # iterativePCA algorithm
  
  GABA_pca <- PCA(res_comp_GABA, ncp = nb_GABA_kfold$ncp)
  plot(GABA_pca)
  plot(GABA_pca, choix="var")
  GABA_wide_PCA$dim1 <- GABA_pca$ind$coord[,1]
  dim1_lai <- lm(data=GABA_wide_PCA, `MPFC` ~ dim1)
  summary(dim1_lai)
  ggplot(GABA_wide_PCA) + aes(y=dim1, x=`L Anterior Insula`) +
    geom_point() + geom_smooth(method="lm") + 
    theme_classic(base_size = 13) +xlab("L Anterior Insula") + ylab("dim1") + 
    theme(legend.key = element_rect(fill = "white", colour = "black"))
  
  
  PCA_measures <- merge(Glu_wide_PCA, GABA_wide_PCA, by = "ld8", suffixes = c("_Glu", "_GABA"))
  PCA_measures <- PCA_measures %>% separate(ld8, c("ID", "Date"), sep = "_" )
  
  EEG_PCA <- merge(alldata_delayOnly, PCA_measures, by = "ID")  
  
  
  # Gamma number of events vs GABA PCA
  EEG_MRS_data_cleanEventGABA <- EEG_PCA %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanEventGABA[EEG_MRS_data_cleanEventGABA$age >= 10 & EEG_MRS_data_cleanEventGABA$age <= 16,], aes(x = Gamma.Event_Number, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_GABA ~ Gamma.Event_Number + age, data = EEG_MRS_data_cleanEventGABA[EEG_MRS_data_cleanEventGABA$age >= 10 & EEG_MRS_data_cleanEventGABA$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs GABA PCA
  EEG_MRS_data_cleanEventVarGABA <- EEG_PCA %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanEventVarGABA[EEG_MRS_data_cleanEventVarGABA$age >= 10 & EEG_MRS_data_cleanEventVarGABA$age <= 16,], aes(x = Gamma.Event_Number_Variability, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_GABA ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_cleanEventVarGABA[EEG_MRS_data_cleanEventVarGABA$age >= 10 & EEG_MRS_data_cleanEventVarGABA$age <= 16,])
  summary(lm.model)
  
  # Gamma power vs GABA PCA
  EEG_MRS_data_cleanPowerGABA <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanPowerGABA[EEG_MRS_data_cleanPowerGABA$age >= 10 & EEG_MRS_data_cleanPowerGABA$age <= 16,], aes(x = Gamma.log1p_Trial_Power, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_GABA ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_cleanPowerGABA[EEG_MRS_data_cleanPowerGABA$age >= 10 & EEG_MRS_data_cleanPowerGABA$age <= 16,])
  summary(lm.model)
  
  # Gamma power variability vs GABA PCA
  EEG_MRS_data_cleanPowerVarGABA <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_GABA - mean(dim1_GABA, na.rm= T)) < (sd(dim1_GABA, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanPowerVarGABA[EEG_MRS_data_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_cleanPowerVarGABA$age <= 16,], aes(x = Gamma.log1p_Trial_Power_Variability, y = dim1_GABA)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_GABA ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_cleanPowerVarGABA[EEG_MRS_data_cleanPowerVarGABA$age >= 10 & EEG_MRS_data_cleanPowerVarGABA$age <= 16,])
  summary(lm.model)
  
  
  #glutamate
  # Gamma number of events vs Glu PCA
  EEG_MRS_data_cleanEventGlu <- EEG_PCA %>% filter(abs(Gamma.Event_Number - mean(Gamma.Event_Number, na.rm= T)) < (sd(Gamma.Event_Number, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanEventGlu[EEG_MRS_data_cleanEventGlu$age >= 10 & EEG_MRS_data_cleanEventGlu$age <= 16,], aes(x = Gamma.Event_Number, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_Glu ~ Gamma.Event_Number + age, data = EEG_MRS_data_cleanEventGlu[EEG_MRS_data_cleanEventGlu$age >= 10 & EEG_MRS_data_cleanEventGlu$age <= 16,])
  summary(lm.model)
  
  # Gamma number of events variability vs Glu PCA
  EEG_MRS_data_cleanEventVarGlu <- EEG_PCA %>% filter(abs(Gamma.Event_Number_Variability - mean(Gamma.Event_Number_Variability, na.rm= T)) < (sd(Gamma.Event_Number_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanEventVarGlu[EEG_MRS_data_cleanEventVarGlu$age >= 10 & EEG_MRS_data_cleanEventVarGlu$age <= 16,], aes(x = Gamma.Event_Number_Variability, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_Glu ~ Gamma.Event_Number_Variability + age, data = EEG_MRS_data_cleanEventVarGlu[EEG_MRS_data_cleanEventVarGlu$age >= 10 & EEG_MRS_data_cleanEventVarGlu$age <= 16,])
  summary(lm.model)
  
  # Gamma power vs Glu PCA
  EEG_MRS_data_cleanPowerGlu <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power - mean(Gamma.log1p_Trial_Power, na.rm= T)) < (sd(Gamma.log1p_Trial_Power, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanPowerGlu[EEG_MRS_data_cleanPowerGlu$age >= 10 & EEG_MRS_data_cleanPowerGlu$age <= 16,], aes(x = Gamma.log1p_Trial_Power, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_Glu ~ Gamma.log1p_Trial_Power + age, data = EEG_MRS_data_cleanPowerGlu[EEG_MRS_data_cleanPowerGlu$age >= 10 & EEG_MRS_data_cleanPowerGlu$age <= 16,])
  summary(lm.model)
  
  # Gamma power variability vs Glu PCA
  EEG_MRS_data_cleanPowerVarGlu <- EEG_PCA %>% filter(abs(Gamma.log1p_Trial_Power_Variability - mean(Gamma.log1p_Trial_Power_Variability, na.rm= T)) < (sd(Gamma.log1p_Trial_Power_Variability, na.rm= T) * 2)) %>% filter(abs(dim1_Glu - mean(dim1_Glu, na.rm= T)) < (sd(dim1_Glu, na.rm= T) * 2))
  
  ggplot(data = EEG_MRS_data_cleanPowerVarGlu[EEG_MRS_data_cleanPowerVarGlu$age >= 10 & EEG_MRS_data_cleanPowerVarGlu$age <= 16,], aes(x = Gamma.log1p_Trial_Power_Variability, y = dim1_Glu)) + geom_point() + stat_smooth(method = "lm")
  
  lm.model <- lm(dim1_Glu ~ Gamma.log1p_Trial_Power_Variability + age, data = EEG_MRS_data_cleanPowerVarGlu[EEG_MRS_data_cleanPowerVarGlu$age >= 10 & EEG_MRS_data_cleanPowerVarGlu$age <= 16,])
  summary(lm.model)
  
}

