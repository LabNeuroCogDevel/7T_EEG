<br>

# MultiScale Entropy Throughout Adolescence 


### Project Lead
Shane D. McKeon

### Faculty Lead
Beatriz Luna 

### Project Start Date
July 2024

### Current Project Status
In progress 

### Datasets
LNCD 7T

### Github Repository
https://github.com/LabNeuroCogDevel/7T_EEG/tree/main/Entropy

## Code Documentation
**Preprocessing:**

Preprocessing can be run using [01_Entropy_Preprocessing.sh](/LabNeuroCogDevel/7T_EEG/blob/main/Entropy/01_Entropy_Preprocessing.sh)

Note this is the same preprocessing as the resting state data in the [Aperiodic EEG Project](https://labneurocogdevel.github.io/7T_EEG/fooofMRS.html)

  ```matlab -nodesktop -r "addpath(genpath('../Preprocessing_Functions/')); run_preprocessing_pipeline('Resting_State')" ```
  
* Initial preprocessing was done using matlab code run_preprocessing_pipeline.m (../Preprocessing_Functions) which first pulls in raw data from 'Raw/EEG/7TBrainMech'.
* Set the task as 'Resting_State' to select the resting state data.
* Bandpass filter between 0.5 Hz and 70 Hz
* Downsamples the data from 1024Hz to 150 Hz
* Removes bad channels (the following criterion is used)
  - arg_flatline: 8
    - Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal.
  - arg_highpass: [0.25 0.75]
    - Transition band for the initial high-pass filter in Hz. This is formatted as [transition-start, transition-end]
  - arg_channel: 0.7
    - Minimum channel correlation. If a channel is correlated at less than this value to a reconstruction of it based on other channels, it is considered abnormal in the given time window. This method requires that channel locations are available and roughly correct; otherwise a fallback criterion will be used.
  - arg_noisy: 5
    - If a channel has more line noise relative to its signal than this value, in standard deviations based on the total channel population, it is considered abnormal.
  - arg_burst: 15
    - Standard deviation cutoff for removal of bursts (via ASR). Data portions whose variance is larger than this threshold relative to the calibration data are considered missing data and will be removed. 
  - arg_window: 0.3
    - Criterion for removing time windows that were not repaired completely. This may happen if the artifact in a window was composed of too many simultaneous uncorrelated sources (for example, extreme movements such as jumps). This is the maximum fraction of contaminated channels that are tolerated in the final output data for each considered window.
* Interpolates missing channels
  - Dataset includes a few subjects that used a 128 cap as opposed to a 64 channel cap. The code removes the 4 channels that are found in 128 but not 64 and reinterpolates the missing channels that were removed from above
* Run ICA to identify eye movements and blinks
* Homogenize Channel locations
  - Read in the channel locations and make sure all files have the correct locations, especially the few subjects who were ran using a 128 channel cap
* Filter out 60 Hz artifact from line noise

<br> 
