
set(0,'DefaultFigureVisible','off'); %set figure visibility to off

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
    idvalues{j} = (setfiles0(j).name(1:14));
end

numSubj = length(idvalues);
subjectERSP = cell(1,50);
subjectITC = cell(1,50);
subjectPowbase= cell(1,50);

%% run evoked activity function 
for i = 251:numSubj
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

save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectERPs_251_end.mat','subjectERSP')
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectITC_251_end.mat', 'subjectITC')
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectPowbase_251_end.mat', 'subjectPowbase')


subjecttfData = cell(1,50);
%% run induced activity function 
for i = 251:numSubj
    tfdata = zeros(121,200,150); % clear previous subjects data by zeroing

    for c = 1:length(channelValues)
        inputfile = setfiles{i};
        [errorSubjects, tfdata] = inducedActivity(i, inputfile, triggerValue, channelValues(c),errorSubjects);
        channeltfdata{c} = tfdata;
    end

    subjecttfData{i} = channeltfdata;
    
    disp(i);
end

save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/trialTFdata_251_end.mat','subjecttfData', '-v7.3')



