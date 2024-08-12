
%% set needed paths
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Functions/SNR_processing')));
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions/preprocessing'));

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

%% run evoked and induced activity function
for i = 1:numSubj
    subject = idvalues{i};
    inputfile = setfiles{i};
    savePath = [outpath '/' subject '_SNRdata.csv'];

    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset(inputfile); % load in eeg file

    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );

    EEG.times = EEG.times + 500; % shift the time so that you are moving the evoked response into the prestim baseline so that we can see how the newtimef function calculates it
    for l = 1:length(EEG.event)
        EEG.event(l).latency = EEG.event(l).latency + 500;
        EEG.urevent(l).latency = EEG.urevent(l).latency + 500;
    end

    EEG = pop_rmdat( EEG, {triggerValue},[-0.2 0.8] ,0); % select events with trigger value you want
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
    EEG = eeg_checkset( EEG );
    
   
    EEG = pop_epoch( EEG, {triggerValue}, [-0.2 0.8], 'epochinfo', 'yes'); % create epochs using selected events

    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
    EEG = eeg_checkset( EEG );

    parpool('local', 4);
    parfor c = 1:64
        [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = pop_newtimef( EEG, 1, c, [-199  799], [2  15] , 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'baseline',[0], 'freqs', [10 70], 'plotphase', 'off', 'padratio', 4, 'winsize', 102);
        channelERSP{c} = ersp;
        channelITC{c} = itc;
        channelPowbase{c} = powbase;
        channeltfdata{c} = tfdata;
        channelfreqs{c} = freqs;
    end
    delete(gcp('nocreate')); %stop the current parallel pool before you start a new one

idx = find(channelfreqs{1,1}>=36 & channelfreqs{1,1}<= 45); % only going to look at the values between the freq of interest

    % baseline  power
for j = 1:length(channelPowbase)
    channel = channelPowbase{1,j};
    selectfreq = channel(:,idx);
    avgChannel = mean(selectfreq,2);
    subjectAvgPowBase500{j} = avgChannel;
end


end
