# Preprocessing all auditory state data for the cortical SNR development project
# Pull in raw data from hera('Raw/EEG/7TBrainMech')
# run_preprocessing_pipeline.m (/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions)
# Set task as 'SNR'
# This only downsamples the data from 1024 Hz to 512 Hz (as opposed to other studies that are downsampled to 150 Hz)
# Fully preprocessed data will be in Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize

matlab -nodesktop -r "addpath(genpath('../Preprocessing_Functions/')); run_preprocessing_pipeline('SNR')"

