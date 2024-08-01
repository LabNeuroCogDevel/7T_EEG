function [] = RunSpectralEvents()

%% For Terminal: addpath(('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))

addpath(('/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))
addpath(genpath('Functions/spectral_event_functions/'));
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20220104'))
ft_defaults

% values to initalize
task = 'Resting_State'; % which task do you want to run. CHANGE ACCORDING TO MGS, REST, SNR, etc. *****
epoch = 'Delay'; % which epoch do you want to run, if REST there is no epoch
band = 'Gamma';
seconds = '3_4';


% In path
inpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data');

% Outpath
outpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/');

% task directory
taskdirectory = [inpath, '/', task];

datapath = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize'];

%load in all the delay files
setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    idvalues(j,:) = (setfiles0(j).name(1:14));
end

numSubj = length(idvalues);
x = cell(1,numSubj);
classLabels = cell(1,numSubj);

% check to see if subject has already been run
[num,txt,raw] = xlsread(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allVariables.csv'));
alreadyRun = txt(:,1);
save(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/alreadyRun.mat'), 'alreadyRun')

%for terminal:
% load(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/alreadyRun.mat'))

%% Create X matrix for Spectral Event Toolbox

[x,subject] = createXmatix(setfiles, task, epoch);

if task == 'MGS' && epoch == 'Delay'
    for i = 1:length(x)
        classLabels{i} = 4+zeros(1,size(x{i},2));
    end

elseif task == 'MGS' && epoch == "Fix"
    for i = 1:length(x)
        classLabels{i} = 2+zeros(1,size(x{i},2));
    end

elseif task == 'Resting_State'
    for i = 1:length(x)
        classLabels{i} = 16130+zeros(1,size(x{i},2));
    end

end

%% Run spectral event toolbox
if band == 'Gamma'
    eventBand = [35,65]; %Frequency range of spectral events
    fVec = (1:70); %Vector of fequency values over which to calculate TFR

elseif band == 'Beta'
    eventBand = [13,30]; %Frequency range of spectral events
    fVec = (1:30); %Vector of fequency values over which to calculate TFR

elseif alpha == 'Alpha'
    eventBand = [8,12]; %Frequency range of spectral events
    fVec = (1:30); %Vector of fequency values over which to calculate TFR

elseif band == 'Theta'
    eventBand = [4,7]; %Frequency range of spectral events
    fVec = (1:30); %Vector of fequency values over which to calculate TFR

end

Fs = 150; % Sampling rate of time-series
findMethod = 1; % Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; % Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);

% %remove emptry cells (from subjects that had already been run)
% x = x(~cellfun('isempty',x));
% Subjects = Subjects(~cellfun('isempty',Subjects));

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis, x, classLabels);

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    save([savefolder band '_SpecEvents_' seconds '.mat'], 'specEv_struct')
    save([savefolder band '_TFR_'seconds '.mat'], 'TFRs');
    save([savefolder band '_subs_' seconds '.mat'], 'Subjects')
    save([savefolder band '_Xmatrix_' seconds '.mat'], 'x')


elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    save([savefolder band '_SpecEvents_FIX.mat'], 'specEv_struct')
    save([savefolder band '_TFR_FIX.mat'], 'TFRs')
    save([savefolder band '_subs_FIX.mat'], 'Subjects')
    save([savefolder band '_Xmatrix_FIX.mat'], 'x')

elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    save([savefolder band '_SpecEvents_RS.mat'], 'specEv_struct')
    save([savefolder band '_TFR_RS.mat'], 'TFRs')
    save([savefolder band '_subs_RS.mat'], 'Subjects')
    save([savefolder band '_Xmatrix_RS.mat'], 'x')

end


end

