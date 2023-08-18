
# Maria Code
# Clean the data

library(tidyverse)
library(ggplot2)
library(readxl)
library(Hmisc)
library(lmerTest)
library(RColorBrewer)
library(data.table)
library(LNCDR)

MRS <- read_csv('H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/13MP20200207_LCMv2fixidx.csv')

# Step 1 Outlier Detection - visual inspection of LCModel fits/spectra
# create a list of who to remove and remove them
lcm <- read_excel("H:/Projects/7TBrainMech/scripts/mri/MRSI_roi/txt/2_lcm.xlsx", col_names = FALSE)

lcm <- separate(lcm, "...1", c("ld8", "junk","y","x"),extra="merge", sep = "[-.]")
lcm <- select(lcm, -junk)
lcm$bad <- TRUE
MRS <- MRS %>% mutate(x=216+1-x,y=216+1-y)
MRS <- merge(MRS, lcm, by=c("ld8", "x", "y"), all=T) 
MRS <- filter(MRS, is.na(bad))
MRS <- select(MRS, -bad)

#keep people's correct coordinates
MRS <- MRS %>% filter(!is.na(roi))

# get rid of junk data noticed recently 
MRS<- MRS %>% filter(Glu.Cr != 0)

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

#### Make dataframes #### 
MRS_glu <- MRS %>%
  filter(Glu.SD <= 20) %>%
  group_by(roi) %>%
  mutate(zscore=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup 

MRS_gaba <- MRS %>%
  filter(GABA.SD <= 20) %>%
  group_by(roi) %>%
  mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup 


MRS_gaba <- MRS_gaba %>% separate(ld8,c("luna","vdate"), remove=F) %>% mutate(visdate = (as.Date(vdate, format = '%Y%m%d')))

MRS_glu <- MRS_glu %>% separate(ld8,c("luna","vdate"), remove=F) %>% mutate(visdate = (as.Date(vdate, format = '%Y%m%d'))) 


lunaize(ggplot(data = MRS_glu[MRS_glu$roi == 13,], aes(x = visdate, y = GMrat)) + geom_point() + stat_smooth(method = 'loess') + geom_line(aes(group=luna), size = 1) + ggtitle("R Thalamus")) + theme(plot.title = element_text(hjust = 0.5))
  
lm.model <- lmer(GMrat ~ visdate + (1 | luna), data = MRS_glu[MRS_glu$roi == 1,])
summary(lm.model)
  
