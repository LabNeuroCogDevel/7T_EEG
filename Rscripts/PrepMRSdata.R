

LoadMRS <- function() {

# Get data and remove bad quality data #
MRS <- read.csv("H:/Projects/7TBrainMech/scripts/eeg/Shane/MRSdata/13MP20200207_LCMv2fixidx.csv")

# Step 1 Outlier Detection - visual inspection of LCModel fits/spectra
# create a list of who to remove and remove them
lcm <- read_excel("H:/Projects/7TBrainMech/scripts/eeg/Shane/MRSdata/2_lcm.xlsx", col_names = FALSE)
lcm <- separate(lcm, "...1", c("ld8", "junk","y","x"),extra="merge", sep = "[-.]")
lcm <- dplyr::select(lcm, -junk)
lcm$bad <- TRUE
MRS <- MRS %>% mutate(x=216+1-x,y=216+1-y)
MRS <- merge(MRS, lcm, by=c("ld8", "x", "y"), all=T) 
MRS <- filter(MRS, is.na(bad))
MRS <- dplyr::select(MRS, -bad)




#keep people's correct coordinates
MRS <- MRS %>% filter(!is.na(roi))
# get rid of junk data
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

z_thres = 3

return(MRS)

}


cleanGlutamate <- function(MRS, roinumber) {
  
  
  library(tidyverse)
  library(gamm4)
  library(readxl)
  library(ggplot2)
  library(lubridate)
  library(mgcViz)
  library(tinytex)
  library(ggpubr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(gamm4)
  library(dplyr)
  library(checkmate)
  
  # good glu dataframe
  z_thres = 4
  
  MRS_glu <- MRS %>%
    filter(Glu.SD <= 20) %>%
    filter(roi == roinumber) %>%
    mutate(zscore_glu=scale(Glu.Cr, center=T, scale=T)) %>%
    filter(abs(zscore_glu) <= z_thres) %>% ungroup %>%  
    mutate(agegrp = cut(age, breaks=c(0,18,Inf), labels=c("10-18","19-30")))
  MRS_glu$vdate <-ymd(gsub(".*_","", MRS_glu$ld8))
  MRS_glu$dateNumeric <- as.numeric(as.POSIXct(MRS_glu$vdate, format="%m-%d-%Y"))
  MRS_glu$relDate <- pmax(0, MRS_glu$dateNumeric - 1560000000)
  
  
  # good gaba dataframe
  MRS_glu <- MRS_glu %>%
    filter(GABA.SD <= 20) %>%
    filter(roi == roinumber)  %>%
    mutate(zscore_gaba=scale(GABA.Cr, center=T, scale=T)) %>%
    filter(abs(zscore_gaba) <= z_thres) %>% ungroup  %>% 
    mutate(agegrp = cut(age, breaks=c(0,18,Inf), labels=c("10-18","19-30")))
  
  

return(MRS_glu)

}

allMRSRegions <- function() {
  ## Look at MRS data
  
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/PrepMRSdata.R", envir = knitr::knit_global(), chdir = TRUE)
  sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)
  MRS <- LoadMRS()
  
  ROIs <- data.frame(MRS$roi, MRS$label)
  
  
  
    ## ACC
  
    roinumber = 7
    
    z_thres = 2
    
    ACC_Glu <- cleanGlutamate(MRS, roinumber)
    
    
    
    ACC_Glu_invage <- lm(data=ACC_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=ACC_Glu, na.action=na.exclude)
    ACC_Glu$Glu_ageResids <- (residuals(var_Glu))
    ACC_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    
    ACC_Glu$Subject <- ACC_Glu$ld8
    
    ACC_Glu <- ACC_Glu %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    ACC_Glu$zscoredGABAMinusGlu <- (ACC_Glu$zscore_gaba) - (ACC_Glu$zscore_glu)
    ACC_Glu$zscoredGABAMinusGlu_abs <- abs(ACC_Glu$zscoredGABAMinusGlu)
    
    ACC_Glu$GABAMinusGlu <- (ACC_Glu$GABA.Cr) - (ACC_Glu$Glu.Cr)
    ACC_Glu$GABAMinusGlu_abs <- abs(ACC_Glu$GABAMinusGlu)
    
    
    
   ## Right posterior Insula
    
    roinumber = 3
    RPI_Glu <- cleanGlutamate(MRS, roinumber)
    
    
    RPI_Glu_invage <- lm(data=RPI_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=RPI_Glu, na.action=na.exclude)
    RPI_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    RPI_Glu$Glu_ageResids <- (residuals(var_Glu))
    
    RPI_Glu$Subject <- RPI_Glu$ld8
    
    RPI_Glu <- RPI_Glu %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    RPI_Glu$zscoredGABAMinusGlu <- (RPI_Glu$zscore_gaba) - (RPI_Glu$zscore_glu)
    RPI_Glu$zscoredGABAMinusGlu_abs <- abs(RPI_Glu$zscoredGABAMinusGlu)
    RPI_Glu$GABAMinusGlu <- (RPI_Glu$GABA.Cr) - (RPI_Glu$Glu.Cr)
    RPI_Glu$GABAMinusGlu_abs <- abs(RPI_Glu$GABAMinusGlu)
    
    
 
    
    ## Left posterior insula
    
    roinumber = 4
    LPI_Glu <- cleanGlutamate(MRS, roinumber)
    
    LPI_Glu_invage <- lm(data=LPI_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=LPI_Glu, na.action=na.exclude)
    LPI_Glu$Glu_ageResids <- (residuals(var_Glu))
    LPI_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    
    LPI_Glu$Subject <- LPI_Glu$ld8
    
    LPI_Glu <- LPI_Glu %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    LPI_Glu$zscoredGABAMinusGlu <- (LPI_Glu$zscore_gaba) - (LPI_Glu$zscore_glu)
    LPI_Glu$zscoredGABAMinusGlu_abs <- abs(LPI_Glu$zscoredGABAMinusGlu)
    LPI_Glu$GABAMinusGlu <- (LPI_Glu$GABA.Cr) - (LPI_Glu$Glu.Cr)
    LPI_Glu$GABAMinusGlu_abs <- abs(LPI_Glu$GABAMinusGlu)
    
    
  
    ## Right anterior insula 
    roinumber = 1
    RAI_Glu <- cleanGlutamate(MRS, roinumber)
    
    
    RAI_Glu_invage <- lm(data=RAI_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=RAI_Glu, na.action=na.exclude)
    RAI_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    RAI_Glu$Glu_ageResids <- (residuals(var_Glu))
    
    RAI_Glu$Subject <- RAI_Glu$ld8
    
    RAI_Glu <- RAI_Glu %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    RAI_Glu$zscoredGABAMinusGlu <- (RAI_Glu$zscore_gaba) - (RAI_Glu$zscore_glu)
    RAI_Glu$zscoredGABAMinusGlu_abs <- abs(RAI_Glu$zscoredGABAMinusGlu)
    RAI_Glu$GABAMinusGlu <- (RAI_Glu$GABA.Cr) - (RAI_Glu$Glu.Cr)
    RAI_Glu$GABAMinusGlu_abs <- abs(RAI_Glu$GABAMinusGlu)
    

    
    ## left anterior insula
    
    roinumber = 2
    LAI_Glu <- cleanGlutamate(MRS, roinumber)
    
    
    LAI_Glu_invage <- lm(data=LAI_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=LAI_Glu, na.action=na.exclude)
    LAI_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    LAI_Glu$Glu_ageResids <- (residuals(var_Glu))
    
    LAI_Glu$Subject <- LAI_Glu$ld8
    
    LAI_Glu <- LAI_Glu %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    LAI_Glu$zscoredGABAMinusGlu <- (LAI_Glu$zscore_gaba) - (LAI_Glu$zscore_glu)
    LAI_Glu$zscoredGABAMinusGlu_abs <- abs(LAI_Glu$zscoredGABAMinusGlu)
    LAI_Glu$GABAMinusGlu <- (LAI_Glu$GABA.Cr) - (LAI_Glu$Glu.Cr)
    LAI_Glu$GABAMinusGlu_abs <- abs(LAI_Glu$GABAMinusGlu)
    
   
 
    
    ## right caudate
    
    roinumber = 5
    RC_Glu <- cleanGlutamate(MRS, roinumber)
    
    RC_Glu_invage <- lm(data=RC_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=RC_Glu, na.action=na.exclude)
    RC_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    RC_Glu$Glu_ageResids <- (residuals(var_Glu))
    
    RC_Glu$Subject <- RC_Glu$ld8
    
    RC_Glu <- RC_Glu %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    RC_Glu$zscoredGABAMinusGlu <- (RC_Glu$zscore_gaba) - (RC_Glu$zscore_glu)
    RC_Glu$zscoredGABAMinusGlu_abs <- abs(RC_Glu$zscoredGABAMinusGlu)
    RC_Glu$GABAMinusGlu <- (RC_Glu$GABA.Cr) - (RC_Glu$Glu.Cr)
    RC_Glu$GABAMinusGlu_abs <- abs(RC_Glu$GABAMinusGlu)
    
   
    ## left caudate
    
    roinumber = 6
    LC_Glu <- cleanGlutamate(MRS, roinumber)
    
    LC_Glu_invage <- lm(data=LC_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=LC_Glu, na.action=na.exclude)
    LC_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    LC_Glu$Glu_ageResids <- (residuals(var_Glu))
    
    LC_Glu$Subject <- LC_Glu$ld8
    
    LC_Glu <- LC_Glu %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    LC_Glu$zscoredGABAMinusGlu <- (LC_Glu$zscore_gaba) - (LC_Glu$zscore_glu)
    LC_Glu$zscoredGABAMinusGlu_abs <- abs(LC_Glu$zscoredGABAMinusGlu)
    LC_Glu$GABAMinusGlu <- (LC_Glu$GABA.Cr) - (LC_Glu$Glu.Cr)
    LC_Glu$GABAMinusGlu_abs <- abs(LC_Glu$GABAMinusGlu)
    
  
    
     ## MPFC
    
    roinumber = 8
    MPFC_Glu <- cleanGlutamate(MRS, roinumber)
    
    
    MPFC_Glu_invage <- lm(data=MPFC_Glu, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=MPFC_Glu, na.action=na.exclude)
    MPFC_Glu$Glu_ageResidsABS <- abs(residuals(var_Glu))
    MPFC_Glu$Glu_ageResids <- (residuals(var_Glu))
    
    MPFC_Glu$Subject <- MPFC_Glu$ld8
    
    MPFC_Glu <- MPFC_Glu %>% separate(ld8,c("luna","vdate"), remove=F)

    
    # looking at GABA - Glu
    MPFC_Glu$zscoredGABAMinusGlu <- (MPFC_Glu$zscore_gaba) - (MPFC_Glu$zscore_glu)
    MPFC_Glu$zscoredGABAMinusGlu_abs <- abs(MPFC_Glu$zscoredGABAMinusGlu)
    MPFC_Glu$GABAMinusGlu <- (MPFC_Glu$GABA.Cr) - (MPFC_Glu$Glu.Cr)
    MPFC_Glu$GABAMinusGlu_abs <- abs(MPFC_Glu$GABAMinusGlu)
    
  
    
 ## right DLPFC
    
    roinumber = 9
    z_thres = 2
    RDLPFC <- cleanGlutamate(MRS, roinumber)
    
    
    RDLPFC_invage <- lm(data=RDLPFC, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=RDLPFC, na.action=na.exclude)
    
    RDLPFC$Glu_ageResidsABS <- abs(residuals(var_Glu))
    RDLPFC$Glu_ageResids <- (residuals(var_Glu))
    
    RDLPFC$Subject <- RDLPFC$ld8
    
    # GABA - Glu 
    RDLPFC$zscoredGABAMinusGlu<- RDLPFC$zscore_gaba - RDLPFC$zscore_glu
    RDLPFC$zscoredGABAMinusGlu_abs <- abs(RDLPFC$zscoredGABAMinusGlu)
    RDLPFC$GABAMinusGlu <- (RDLPFC$GABA.Cr) - (RDLPFC$Glu.Cr)
    RDLPFC$GABAMinusGlu_abs <- abs(RDLPFC$GABAMinusGlu)
    
    RDLPFC<- RDLPFC %>% separate(Subject,c("luna","vdate"), remove=F)
    
    #GABA - Glu vs age
    
    RDLPFC$zscoredGABAMinusGlu <- RDLPFC$zscore_gaba - RDLPFC$zscore_glu
    
    RDLPFC$zscoredGABAMinusGlu_abs <- abs(RDLPFC$zscoredGABAMinusGlu)
    

    
     ## left DLPFC
    
    roinumber = 10
    LDLPFC <- cleanGlutamate(MRS, roinumber)
    
    
    LDLPFC_invage <- lm(data=LDLPFC, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=LDLPFC, na.action=na.exclude)
    
    LDLPFC$Glu_ageResidsABS <- abs(residuals(var_Glu))
    LDLPFC$Glu_ageResids <- (residuals(var_Glu))
    
    LDLPFC$Subject <- LDLPFC$ld8
    
    # GABA - Glu 
    LDLPFC$zscoredGABAMinusGlu<- LDLPFC$zscore_gaba - LDLPFC$zscore_glu
    LDLPFC$zscoredGABAMinusGlu_abs <- abs(LDLPFC$zscoredGABAMinusGlu)
    LDLPFC$GABAMinusGlu <- (LDLPFC$GABA.Cr) - (LDLPFC$Glu.Cr)
    LDLPFC$GABAMinusGlu_abs <- abs(LDLPFC$GABAMinusGlu)
    
    LDLPFC<- LDLPFC %>% separate(Subject,c("luna","vdate"), remove=F)
  
  
     ## right thalamus
    
    roinumber = 13
    Rthalamus <- cleanGlutamate(MRS, roinumber)
    
    
    Rthalamus_invage <- lm(data=Rthalamus, Glu.Cr ~ invage + sex + GMrat)
    #summary(ROI_Glu_invage)
    
    var_Glu <- lm(Glu.Cr ~ invage + sex + GMrat, data=Rthalamus, na.action=na.exclude)
    Rthalamus$Glu_ageResidsABS <- abs(residuals(var_Glu))
    Rthalamus$Glu_ageResids <- (residuals(var_Glu))
    
    Rthalamus$Subject <- Rthalamus$ld8
    
    Rthalamus <- Rthalamus %>% separate(ld8,c("luna","vdate"), remove=F)
    
    # looking at GABA - Glu
    Rthalamus$zscoredGABAMinusGlu <- (Rthalamus$zscore_gaba) - (Rthalamus$zscore_glu)
    Rthalamus$zscoredGABAMinusGlu_abs <- abs(Rthalamus$zscoredGABAMinusGlu)
    Rthalamus$GABAMinusGlu <- (Rthalamus$GABA.Cr) - (Rthalamus$Glu.Cr)
    Rthalamus$GABAMinusGlu_abs <- abs(Rthalamus$GABAMinusGlu)
    
    
    lunaize(ggplot(data = Rthalamus , aes(y = zscore_glu, x = age)) + stat_smooth(method='gam') + xlab("Age") + ylab("abs(RPI GABA-Glu)") + geom_point() )
    
    gam.model <- gam(zscore_glu ~ s(age, k=4),
                     data = Rthalamus  )
    
    summary(gam.model)
    
    vars <- c("Subject", "luna", "visitnum", "label", "age", "Glu.Cr", "GABA.Cr", "zscore_glu", "zscore_gaba", "Glu_ageResids", "Glu_ageResidsABS", "zscoredGABAMinusGlu", "zscoredGABAMinusGlu_abs", "GABAMinusGlu", "GABAMinusGlu_abs")
    ACCdf <- ACC_Glu[vars]
    MPFCdf <- MPFC_Glu[vars]
    RDLPFCdf <- RDLPFC[vars]
    LDLPFCdf <- LDLPFC[vars]
    RPIdf <- RPI_Glu[vars]
    LPIdf <- LPI_Glu[vars]
    RAIdf <- RAI_Glu[vars]
    LAIdf <- LAI_Glu[vars]
    LCdf <- LC_Glu[vars]
    RCdf <- RC_Glu[vars]
    Rthalamusdf <- Rthalamus[vars]
    
    
    allMRSregions <- rbind(ACCdf, MPFCdf) %>% rbind(., RDLPFCdf) %>% rbind(., LDLPFCdf) %>% rbind(., Rthalamusdf) %>% mutate(visitno = visitnum)
    
  return(allMRSregions)
}
