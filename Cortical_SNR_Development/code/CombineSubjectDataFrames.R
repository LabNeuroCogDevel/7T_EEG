#!/usr/bin/env Rscript

# Initialize directory and empty dataframes ----

combineSubjectDataframes <- function (freq){
oldpwd <- getwd()
# Set your working directory to the folder containing your CSV files
setwd(paste0("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/results/individualSubFiles/", freq,"hz/"))

# List all CSV files in the directory
csv_files <- list.files(pattern = paste0(freq,"Hz.csv"))

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it ----
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

# Prep the file for merge 7t eeg ----

# rename channel to urchan to merge with channel labels and coordinates 
names(combined_data)[names(combined_data) == "Channel"] <- "urchan"

write.csv(combined_data, file = paste0(oldpwd, '/results/allSubjectsSNR_allChans_',freq,'Hz.csv'), row.names = F)

}
