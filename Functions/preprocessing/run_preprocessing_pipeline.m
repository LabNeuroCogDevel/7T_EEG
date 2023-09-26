
% load in paths to feildtrip and eeglab 
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions/preprocessing'));
addpath(genpath('/resources/Euge/'))
run('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1/eeglab.m');

% Outpath
maindir = hera('Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data');

% which task do you want to run. CHANGE ACCORDING TO MGS, REST, SNR, etc. 
task = 'SNR'; 
  
disp(["i am running" task]) 
   
taskdirectory = [maindir, '/', task]; 

% initial values
lowBP = 0.5;
topBP = 70;
FLAG = 1;

%% settings
only128 = 0; % 0==do all, 1==only 128 channel subjects
condition = 1; %0 - if you want to overwrite an already existing file; 1- if you want it to skip subjects who have already been run through singlesubject

remark(task, taskdirectory); % change the trigger values to be single digit 


% gather all file paths for all subjects remarked data 
setfilesDir = [taskdirectory, '/remarked/1*_20*.set'];    
setfiles = all_remarked_set(setfilesDir);

n = size(setfiles,1); %number of EEG sets to preprocess

% loop through every subject to preprocess
for i = 1:n
    inputfile = setfiles{i};
    try
      preprocessing_pipeline(inputfile, taskdirectory, lowBP, topBP, FLAG, condition, task)
   catch e
      fprintf('Error processing "%s": %s\n',inputfile, e.message)
      for s=e.stack
         disp(s)
      end
   end
end


