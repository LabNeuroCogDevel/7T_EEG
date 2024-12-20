function [] = run_preprocessing_pipeline(task)

% load in paths to feildtrip and eeglab 
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/'));
addpath(genpath('/resources/Euge/'))
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
% addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20220104/')
%run('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1/eeglab.m');
run('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/eeglab.m')


% Outpath
% maindir = hera('Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data');
maindir = hera('Abby/preprocessed_data');

disp(["i am running" task]) 
   
taskdirectory = [maindir '/' task]; 

% initial values
lowBP = 0.5;
topBP = 70;
FLAG = 1;

%% settings
only128 = 0; % 0==do all, 1==only 128 channel subjects
condition = 1; %0 - if you want to overwrite an already existing file; 1- if you want it to skip subjects who have already been run through singlesubject
dryrun = 1; 

remark(task, taskdirectory, dryrun); % change the trigger values to be single digit 
remarkAS(taskdirectory,dryrun) % remarking antisaccade trials

% gather all file paths for all subjects remarked data 
setfilesDir = [taskdirectory, '/remarked/1*_20*.set'];    
setfiles = all_remarked_set(setfilesDir);
n = size(setfiles,1); %number of EEG sets to preprocess

% some participants are left out from lines 33 and 34 for AS
if task=="anti"
    setfilesDir = [taskdirectory, '/remarked/*.set'];
    setfiles = all_remarked_setAS(setfilesDir);
    n = length(setfiles);
end

% loop through every subject to preprocess
for i = 1:n
    inputfile = setfiles{i};
    try
        if task ~= "anti"
            preprocessing_pipeline(inputfile, taskdirectory, lowBP, topBP, FLAG, condition, task)
        elseif task == "anti"
            preprocessing_pipelineAS(inputfile,taskdirectory,lowBP,topBP,FLAG,condition,task)
        end
   catch e
      fprintf('Error processing "%s": %s\n',inputfile, e.message)
      for s=e.stack
         disp(s)
      end
   end
end


%% select ICA values to reject

ICA_Path = [taskdirectory '/ICAwhole'];
CleanICApath = [taskdirectory '/AfterWhole/ICAwholeClean/'];


EEGfileNames = dir([ICA_Path '/*.set']);

for fidx = 1:length(EEGfileNames)
    filename = EEGfileNames(fidx).name;
    locs = file_locs(fullfile(ICA_Path,filename), taskdirectory, task);
    if exist(locs.ICAwholeClean, 'file')
        fprintf('skipping; already created %s\n', locs.ICAwholeClean);
        continue
    end

    selectcompICA
end


%% Homogenize Chanloc
datapath = [taskdirectory '/AfterWhole/ICAwholeClean'];
savepath = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize'];

setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};
for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
end

correction_cap_location = hera('Projects/7TBrainMech/scripts/eeg/Shane/resources/ELchanLoc.ced');
for i = 1:length(setfiles)
    homogenizeChanLoc(setfiles{i},correction_cap_location,savepath, taskdirectory, task)
end


%% filter out the 60hz artifact from electronics 
datapath = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize'];

setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};
for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
end

for i = 1:length(setfiles)
    inputfile = setfiles{i};
    [filepath,filename ,ext] =  fileparts((setfiles{i}));

    EEG = pop_loadset(inputfile);
    EEG = pop_eegfiltnew(EEG, 59, 61, [], 1, [], 0);
    EEG = pop_saveset(EEG, 'filename', filename, 'filepath', datapath);
end


if task == "MGS" || task == "anti"
    %% Clean epochs to remove
    path_data = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize/'];
    epoch_path = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize/'];
    epoch_folder = [taskdirectory '/AfterWhole/epoch/'];
    epoch_rj_marked_folder = [taskdirectory '/AfterWhole/epochclean/'];

    EEGfileNames = dir([path_data, '/*_icapru.set']);

    revisar = {};
    for currentEEG = 1:size(EEGfileNames,1)
        if currentEEG == 284
            continue
        end
        filename = [EEGfileNames(currentEEG).name];
        inputfile = [epoch_path,filename];
        revisar{currentEEG} = epochlean(inputfile,epoch_folder,epoch_rj_marked_folder,task,taskdirectory);
    end
end

