
%% For Terminal: addpath(('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))
%if class labels doesnt work:

addpath(('/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))
addpath(genpath('Functions'));
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20220104'))
ft_defaults

datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Resting_State/AfterWhole/ICAwholeClean_homogenize');

%load in all the delay files
setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    
    idvalues(j,:) = (setfiles0(j).name(1:14));
end

numSubj = length(idvalues);
x = cell(1,numSubj);
classLabels = cell(1,numSubj);

% check to see if subject has already been run
[num,txt,raw] = xlsread(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_allVariables.csv'));
alreadyRun = txt(:,1);
save(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/alreadyRun.mat'), 'alreadyRun')

%for terminal:
% load(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/alreadyRun.mat'))


%% Create X matrix for Spectral Event Toolbox
clear x
clear Subjects

for i = 1:numSubj
    inputfile = setfiles{i};
    subject = inputfile(93:106);
    
    %     if any(strcmp(alreadyRun, subject))
    %         warning('%s already complete', subject)
    %         continue
    %
    %     else
    
    inputfile = setfiles{i};
    
    clear events
    clear data
    clear hdr
    clear avgData
    clear classLabels
    clear allData
    
    hdr = ft_read_header(inputfile);
    data = ft_read_data(inputfile, 'header', hdr);
    events = ft_read_event(inputfile, 'header', hdr);
    
    cfg = [];
    cfg.dataset = inputfile;
    cfg.headerfile = inputfile;
    cfg.channel =   {'all', '-POz'};
    cfg.trialdef.eventtype = 'trigger';
    cfg.trialdef.eventvalue = '16130';
    cfg.trialdef.prestim = 0;
    cfg.trialdef.poststim = 4;
    cfg.trialfun = 'ft_trialfun_general';
    cfg.trialdef.ntrials = 139;
    cfg.event = events;
    
    if ~ischar(cfg.event(1).value)
        
        newevents = cfg.event;
        for j = 1:length(newevents)
            newevents(j).value = num2str(newevents(j).value);
        end
        cfg.event = newevents;
        
    end
    
    try
        
        [cfg]= ft_definetrial(cfg);
        
    catch
        warning('%s wont run through feildtrip', subject)
        wontRun{i} = subject;
        continue;
    end
    
    
    [data] = ft_preprocessing(cfg);
    
    
    %redefine the trails to be between 2-3 of each eyes open trigger
    cfg.trl = [];
    cfg = rmfield(cfg, 'trl');
    cfg.toilim = [2 3];
    data = ft_redefinetrial(cfg, data);
    
%     for k = 1:size(data.label,1) %channels
        for j = 1:length(data.trial) %trial
            trialData = (data.trial{1,j});
            avgTrialData = mean(trialData, 1); 
            allData(j,:) = avgTrialData;
        end
        
%     end
    
%     avgData = squeeze(mean(allData, 1));
    
    x{i} = allData';
    Subjects{i} = subject;
end

%% Gamma Analysis 
eventBand = [30,75]; %Frequency range of spectral events
fVec = (30:75); %Vector of fequency values over which to calculate TFR
Fs = 150; %Sampling rate of time-series
findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; %Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);
% 
% %remove emptry cells (from subjects that had already been run)
 x = x(~cellfun('isempty',x));
% Subjects = Subjects(~cellfun('isempty',Subjects));

for i = 1:length(x)
    classLabels{i} = 16130+zeros(1,size(x{i},2));
end

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  x ,classLabels);


save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_SpecEvents_RS.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Resting_State/Gamma_TFR_RS.mat'), 'TFRs')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Resting_State_Subs_RS.mat'), 'Subjects')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Resting_State_xMatrix.mat'), 'x')


%% Theta Analysis
eventBand = [4,7]; %Frequency range of spectral events
fVec = (1:30); %Vector of fequency values over which to calculate TFR
Fs = 150; %Sampling rate of time-series
findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; %Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);
%             classLabels = 4+zeros(1,size(avgData,1));

%remove emptry cells (from subjects that had already been run)
% x = x(~cellfun('isempty',x));
% Subjects = Subjects(~cellfun('isempty',Subjects));

for i = 1:length(x)
    classLabels{i} = 4+zeros(1,size(x{i},2));
end


[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  x ,classLabels);

save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_SpecEvents_RS.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_TFR_RS.mat'), 'TFRs')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Resting_State/Theta_RS.mat'), 'Subjects')


%% Beta Analysis
eventBand = [13,30]; %Frequency range of spectral events
fVec = (1:30); %Vector of fequency values over which to calculate TFR
Fs = 150; %Sampling rate of time-series
findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; %Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);
%             classLabels = 4+zeros(1,size(avgData,1));

%remove emptry cells (from subjects that had already been run)
x = x(~cellfun('isempty',x));
Subjects = Subjects(~cellfun('isempty',Subjects));

for i = 1:length(x)
    classLabels{i} = 2+zeros(1,size(x{i},2));
end


[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  x ,classLabels);

save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_SpecEvents_RS.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_TFR_RS.mat'), 'TFRs')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Resting_State/Beta_newSubs_RS.mat'), 'Subjects')



%% Alpha Analysis


eventBand = [8,12]; %Frequency range of spectral events
fVec = (1:30); %Vector of fequency values over which to calculate TFR
Fs = 150; %Sampling rate of time-series
findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; %Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);
%         classLabels = 4+zeros(1,size(avgData,1));

%remove emptry cells (from subjects that had already been run)
x = x(~cellfun('isempty',x));
Subjects = Subjects(~cellfun('isempty',Subjects));

for i = 1:length(x)
    classLabels{i} = 2+zeros(1,size(x{i},2));
end

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  x ,classLabels);

save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_SpecEvents_RS.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_TFR_RS.mat'), 'TFRs')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Resting_State/Alpha_newSubs_RS.mat'), 'Subjects')




