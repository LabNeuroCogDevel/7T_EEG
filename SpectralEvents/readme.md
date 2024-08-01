To run the spectral events pipeline on all channels averaged together

1. [RunSpectralEvents.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/spectral_event_functions/RunSpectralEvents.m)
   - Set desired task (Resting_State, MGS, SNR)
   - Set epoch (Delay or Fix; only is running resting state)
   - Set seconds; only if running delay epoch
   - Set band to analyze (gamma, beta, theta, alpha)
   - Pull subjects from the preprocessed folder
   - Create the X matrix ([createXmatrix.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/spectral_event_functions/createXmatrix.m)) that organizes the data in the correct format that the spectral event code expects
   - Creates an array of values corresponding to the trigger values, called classLabels
   - Defines event band and frequency values based on which band you're running
   - Runs the spectral event toolbox, see [Spectral Event Toolbox](https://github.com/jonescompneurolab/SpectralEvents)
   - Saves out large structs of the subjects spectral event data

2. [SpectralEvents_Struct_Extraction.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/spectral_event_functions/SpectralEvents_Struct_Extraction.m)
   - Initialize task, epoch, and band
   - Load in the structs
   - Extracts out the Event Number, Power, Duration, Max Power, and Max Frequency and generates the average of each
   - Crates a large table with all the subjects data


To run the spectral events pipeline on each individual channel 

1. [RunSpectralEvents_IndividualChannels.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/spectral_event_functions/RunSpectralEvents_IndividualChannels.m)
   - Same as above but on every channel
   - Calls [createXmatrix_individualChannels.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/spectral_event_functions/createXmatrix_individualChannels.m)
2. [SpectralEvents_Struct_Extraction_IndividualChannels.m](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Functions/spectral_event_functions/SpectralEvents_Struct_Extraction_IndividualChannels.m)
   - Same as above 

