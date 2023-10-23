
set(0,'DefaultFigureVisible','off'); %set figure visibility to off
%% set needed paths
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Functions/SNR_processing')));

%% set initial values 
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize');
triggerValue = '4'; 
% channelValues = [4,5,6,36,37,38]; if you only want to run DLPFC
outpath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR');

%% load in all the data files
setfiles0 = dir([datapath,'/*icapru*.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    idvalues{j} = (setfiles0(j).name(1:14));
end

numSubj = length(idvalues);
errorSubjects = cell(1,numSubj);

%% run evoked and induced activity function 
for i = 1:numSubj
    subject = idvalues{i};
    inputfile = setfiles{i};
    savePath = [outpath '/' subject '_SNRdata.csv'];

    if ~exist(savePath, 'file')
        fprintf('The file %s does not exist, running SNR code\n', savePath)

       [errorSubjects] = evokedInducedActivity(i, inputfile, triggerValue,outpath,subject,errorSubjects);
    end

    disp(i);
end

