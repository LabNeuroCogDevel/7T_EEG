Instructions to run preprocessing on 7T EEG data

Note: you need the [hera.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/hera.m) function for the pipeline to run

1. [run_preprocessing_pipeline.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/preprocessing/run_preprocessing_pipeline.m)
   - [remark.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/preprocessing/remark.m) function pulls the raw data from the database and remarks the trigger vales to single digits corresponding to the appropriate epoch
2. [preprocessing_pipeline.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/preprocessing/preprocessing_pipeline.m)
    - Loads in the desired subject
    - Loads in the file containing the coordinates of the cap electrodes
    - Checks to see if the subject has already been preprocessed. If yes, skips and moves on to the next
      - [file_locs.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/preprocessing/file_locs.m) will load in all the path and file names 
    - References the data to the electrodes corresponding to the mastoids
    - Filters out the data below 0.5 Hz (low frequency drifts) and 70 Hz
    - Resample the data to 150 Hz (all tasks besides SNR) or 512 Hz (for SNR task)
    - Removes external channels 
    - Imports channel coordinates into the EEG file
    - Bad channel rejection
      - See [clean_rawdata.m](https://github.com/sccn/clean_rawdata) EEGLAB plugin
   - Interpolate missing channels
      - See [pop_interp.m](https://github.com/INCF/p3-validator/blob/master/trunk/lib/eeglab9_0_4_5s/functions/popfunc/pop_interp.m)
   - Reference data to average for ICA
3. [runICAs](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/preprocessing/runICAs.m) for eye movement
4. Remove bad epochs via [epochclean.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/preprocessing/epochlean.m)
   - Uses [pop_epoch.m](https://github.com/sccn/eeglab/blob/develop/functions/popfunc/pop_epoch.m)
5. Move channel locations if subject used a 128 channel cap so that the coordinates are the same as the 64 channel cap
   - See [homogenizeChanLoc.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/preprocessing/homogenizeChanLoc.m)
6. Filter out the 60 Hz artifact 
