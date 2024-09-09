

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
library(eegUtils)
library(tvem)
library(interactions)
library(akima)
library(mice)


entropy <- read.csv('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE Lab/entropy/allSubjectsEntropy.csv')

nonEpContacts <- read.csv('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE Lab/nonEPcontacts.csv')

# List of subjects
subjects <-intersect(unique(nonEpContacts$Subject), unique(entropy$Subject))


# Initialize an empty data frame to store the results
result <- data.frame()

# Loop through subjects
for (subject in subjects) {
  # Select contacts for the current subject
  selectContacts <- nonEpContacts %>%
    filter(Subject == subject) %>%
    select(nonEpContacts)
  
  # Create a formatted string
  formattedString <- paste(selectContacts$nonEpContacts, collapse = '|')
  
  # Filter rows in entropy for the current subject
  filteredData <- entropy[grep(formattedString, entropy$label), ] %>%
    subset(., !grepl('ET', label)) %>%
    subset(., !grepl('EEG', label)) %>%
    filter(Subject == subject)
  
  # Append the result to the overall result data frame
  subjectResult <- bind_rows(result, filteredData)
  
  result <- rbind(result, subjectResult)
}



