%% set needed paths
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Functions/SNR_processing')));

%% set initial values 
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize');
triggerValue = '4'; 
channelValues = [4,5,6,36,37,38];
errorSubjects = [];

%% load in all the data files
setfiles0 = dir([datapath,'/*icapru*.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    idvalues(j,:) = (setfiles0(j).name(1:14));
end

numSubj = length(idvalues);


%% run evoked activity function 
for i = 1:numSubj
    ersp = zeros(121,200); % clear previous subjects data by zeroing
    itc = zeros(121,200); % clear previous subjects data by zeroing
    powbase = zeros(1,121); % clear previous subjects data by zeroing

    for c = 1:length(channelValues)
        inputfile = setfiles{i};
        [ersp,itc,powbase,times,freqs, errorSubjects] = evokedActivity(i, inputfile, triggerValue, channelValues(c),errorSubjects);
        channelERSP{c} = ersp; 
        channelITC{c} = itc; 
        channelPowbase{c} = powbase; 
    end

    subjectERSP{i} = channelERSP;
    subjectITC{i} = channelITC;
    subjectPowbase{i} = channelPowbase;

    disp(i);
end
