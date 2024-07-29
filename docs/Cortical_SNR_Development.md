# Cortical SNR Development

## Reproducing
Code starts on 
[`01_Cortical_SNR_Preprocessing.sh`](https://github.com/LabNeuroCogDevel/7T_EEG/blob/main/Cortical_SNR_Development/01_Cortical_SNR_Preprocessing.sh)

```matlab
addpath(genpath('../Preprocessing_Functions/'));  % (1)!
run_preprocessing_pipeline('SNR');
```

1. These are shared preprocessing functions to process raw `.bdf` files into ica corrected `.set` fieldtrip files


