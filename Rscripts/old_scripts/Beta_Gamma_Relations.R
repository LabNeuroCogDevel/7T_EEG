

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2,lme4, multcomp, tidyr, tidyverse,wesanderson, reshape2, corrplot, ggExtra, cowplot, ggsci, mgcv, ggpointdensity, devtools, ggpubr, itsadug, useful, rmarkdown, tinytex, gratia, gamm4, car, GGally, CCA)

#devtools::install_github('LabneuroCogDevel/LNCDR')
#library(devtools)
#devtools::install_github("strengejacke/sjPlot")
#library(LNCDR)
library(data.table)
require(dplyr)
library(factoextra)
library(tidyverse)
require(knitr)


# load in files
allBeta <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Delay1_2/Beta_all_data.csv')
allGamma <-  read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Delay1_2/Gamma_all_data.csv')

