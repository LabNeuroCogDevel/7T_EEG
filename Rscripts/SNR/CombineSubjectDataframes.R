# Libraries ----

library(LNCDR)
library(data.table)
library(dplyr)
library(factoextra)
library(ggplot2)
library(e1071)
library(caret)
attach(mtcars)
library(grid)
library(gridExtra)
library(plotrix)
library(mgcv)
library(readxl)
library(lme4)
library(lubridate)
library(checkmate)
library(lmerTest)
library(tidyr)
library(jtools)

# Initialize directory and empty dataframes ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/totalEvokedInduced_indivSubjects/")

# List all CSV files in the directory
csv_files <- list.files(pattern = ".csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it ----
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}


write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsTotalEvokedInduced_40hz.csv', row.names = F)


# Prep the file for merge 7t eeg ----

# rename channel to urchan to merge with channel labels and coordinates 
names(combined_data)[names(combined_data) == "Channel"] <- "urchan"


write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsSNR.csv', row.names = F)
