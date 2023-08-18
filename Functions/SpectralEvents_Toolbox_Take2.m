
%% For Terminal: addpath(('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))
%if class labels doesnt work:
for i = 1:length(x)
    classLabels{i} = 2+zeros(1,size(x{i},2));
end


%%


addpath(('/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))
addpath(genpath('Functions'));
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20191127'))
ft_defaults

datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/ICAwholeClean_homogenize');

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

%% check to see if subject has already been run
% Gamma Analysis 

clear x
clear Subjects

[num,txt,raw] = xlsread(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allVariables.csv'));
alreadyRun = txt(:,1);
save(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/alreadyRun.mat'), 'alreadyRun')

%for terminal:
% load(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/alreadyRun.mat'))


for i = 1:numSubj
    inputfile = setfiles{i};
    subject = inputfile(84:97);
    
    if any(strcmp(alreadyRun, subject))
        warning('%s already complete', subject)
        continue
        
    else
    
        inputfile = setfiles{i};
        
        clear events
        clear data
        clear hdr
        clear avgData
        clear classLabels
        
        hdr = ft_read_header(inputfile);
        data = ft_read_data(inputfile, 'header', hdr);
        events = ft_read_event(inputfile, 'header', hdr);
        
        cfg = [];
        cfg.dataset = inputfile;
        cfg.headerfile = inputfile;
        cfg.channel =   {'all', '-POz'};
        cfg.trialdef.eventtype = 'trigger';
        cfg.trialdef.eventvalue = '2';
        cfg.trialdef.prestim = 0;
        cfg.trialdef.poststim = 1;
        cfg.trialfun = 'ft_trialfun_general';
        cfg.trialdef.ntrials = 95;
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
        
        
        %redefine the trails to be between 1-2 seconds of the delay period
        cfg.trl = [];
        cfg = rmfield(cfg, 'trl');
        cfg.toilim = [0 1];
        data = ft_redefinetrial(cfg, data);
        
        for k = 1:size(data.label,1) %channels
            for j = 1:length(data.trial) %trial
                trialData = (data.trial{1,j});
                channel = trialData(k,:);
                avgData(j,:) = channel;
            end
            
            x{i} = avgData';
            Subjects{i} = subject;
            
        end
    end
end
%%

eventBand = [30,75]; %Frequency range of spectral events
fVec = (30:75); %Vector of fequency values over which to calculate TFR
Fs = 150; %Sampling rate of time-series
findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; %Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);

%remove emptry cells (from subjects that had already been run)
x = x(~cellfun('isempty',x));
Subjects = Subjects(~cellfun('isempty',Subjects));

for i = 1:length(x)
    classLabels{i} = 2+zeros(1,size(x{i},2));
end

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  x ,classLabels);


save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_SpecEvents_newSubs_FIX.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_TFR_newSubs_FIX.mat'), 'TFRs')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Fix/Gamma_newSubs_FIX.mat'), 'Subjects')




% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/GammaAvgEventNumber_allevents.mat'), 'AvgEventNumber')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/GammaAvgEventDuration_allevents.mat'), 'AvgEventDuration')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/GammaAvgPower_allevents.mat'), 'AvgPower')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/GammaPeakPower_allevents.mat'), 'peakPower')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/GammaSDPower_allevents.mat'), 'SDPower')


addpath(('Z:/DB_SQL/'));

all_ages = db_query('select id ||''_''|| to_char(vtimestamp,''YYYYMMDD'') as IDvalues, age, visitno from visit natural join enroll where etype like ''LunaID''  and vtype ilike ''eeg''  order by id, age ');
all_ages.idvalues= char(all_ages.idvalues);
writetable(all_ages, hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/agefile.csv'));


for j = 1 : length(setfiles0)
    
    idvalues(j,:) = (setfiles0(j).name(1:14));
end

MeanAvgEventNumber = mean(AvgEventNumber,1)';
MeanAvgEventDuration = mean(AvgEventDuration,1)';
MeanAvgPower = mean(AvgPower, 1)';
MeanpeakPower = mean(peakPower, 1)';


T = table(idvalues, MeanAvgEventNumber, MeanAvgEventDuration, MeanAvgPower, MeanpeakPower);
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_allevents.mat'), 'T')

all_info = innerjoin(T, all_ages);



all_info = all_info(all_info.MeanAvgPower ~= 0, :);

writetable(all_info, hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_Spectral_Analysis_Table_all_allevents.csv'));



DLPFC_AvgEventNumber = AvgEventNumber(6,:)';
DLPFC_AvgEventDuration = AvgEventDuration(6,:)';
DLPFC_AvgPower = AvgPower(6,:)';
DLPFC_peakPower = peakPower(6,:)';

DLPFC_T = table(idvalues, DLPFC_AvgEventNumber, DLPFC_AvgEventDuration, DLPFC_AvgPower, DLPFC_peakPower);
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Gamma_Spectral_Analysis_Table_allevents.mat'), 'T')


all_info = innerjoin(DLPFC_T, all_ages);


all_info = all_info(all_info.DLPFC_AvgEventNumber ~= 0, :);

writetable(all_info, hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Gamma_Spectral_Analysis_Table_all_allevents.csv'));


%% Theta Analysis
numSubj = length(idvalues);
x = cell(1,numSubj);
classLabels = cell(1,numSubj);
% SE_Subs_Channel = cell(length(setfiles), 63);

% clear AvgEventDuration
% clear AvgEventFSpan
% clear AvgEventNumber
% clear AvgPower
% clear peakPower

% check to see if subject has already been run

clear x
clear Subjects

[num,txt,raw] = xlsread(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_allVariables.csv'));
alreadyRun = txt(:,1);
save(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/alreadyRun.mat'), 'alreadyRun')

%for terminal:
% load(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/alreadyRun.mat'))



for i = 1 : numSubj
    inputfile = setfiles{i};
    subject = inputfile(84:97);
    
    if any(strcmp(alreadyRun, subject))
        warning('%s already complete', subject)
        continue
        
    else
        
        %     clear events
        %     clear data
        %     clear hdr
        %     clear avgData
        %     clear classLabels
        
        hdr = ft_read_header(inputfile);
        data = ft_read_data(inputfile, 'header', hdr);
        events = ft_read_event(inputfile, 'header', hdr);
        
        cfg = [];
        cfg.dataset = inputfile;
        cfg.headerfile = inputfile;
        cfg.channel = {'all', '-POz'};
        cfg.trialdef.eventtype = 'trigger';
        cfg.trialdef.eventvalue = '2';
        cfg.trialdef.prestim = 0;
        cfg.trialdef.poststim = 1;
        cfg.trialfun = 'ft_trialfun_general';
        cfg.trialdef.ntrials = 95;
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
            continue;
            warning('%s wont run through feildtrip', subject)
            wontRun{i} = subject;
        end
        
        
        [data] = ft_preprocessing(cfg);
        
        
        %redefine the trails to be between 1-2 seconds of the delay period
        cfg.trl = [];
        cfg = rmfield(cfg, 'trl');
        cfg.toilim = [0 1];
        data = ft_redefinetrial(cfg, data);
        
        for k = 1:size(data.label,1)
            for j = 1:length(data.trial)
                trialData = (data.trial{1,j});
                channel = trialData(k,:);
                avgData(j,:) = channel;
            end
            
            x{i} = avgData';
            Subjects{i} = subject;
        end
    end
end

eventBand = [4,7]; %Frequency range of spectral events
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
%
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_SpecEvents_newSubs_FIX.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_TFR_newSubs_FIX.mat'), 'TFRs')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Fix/Theta_newSubs_FIX.mat'), 'Subjects')



% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/ThetaAvgEventNumber_allevents.mat'), 'AvgEventNumber')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/ThetaAvgEventDuration_allevents.mat'), 'AvgEventDuration')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/ThetaAvgPower_allevents.mat'), 'AvgPower')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/ThetaPeakPower_allevents.mat'), 'peakPower')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/ThetaSDPower_allevents.mat'), 'SDPower')

%
AgeFile = xlsread(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/ERPs/EdaDI.xlsx'));

addpath(('Z:/DB_SQL/'));

all_ages = db_query('select id ||''_''|| to_char(vtimestamp,''YYYYMMDD'') as IDvalues, age, visitno from visit natural join enroll where etype like ''LunaID''  and vtype ilike ''eeg''  order by id, age ');
all_ages.idvalues= char(all_ages.idvalues);

for j = 1 : length(setfiles0)
    
    idvalues(j,:) = (setfiles0(j).name(1:14));
end
%
MeanAvgEventNumber = mean(AvgEventNumber,1)';
MeanAvgEventDuration = mean(AvgEventDuration,1)';
MeanAvgPower = mean(AvgPower, 1)';
MeanpeakPower = mean(peakPower, 1)';

T = table(idvalues, MeanAvgEventNumber, MeanAvgEventDuration, MeanAvgPower, MeanpeakPower);
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_allevents.mat'), 'T')

all_info = innerjoin(T, all_ages);



all_info = all_info(all_info.MeanAvgEventNumber ~= 0, :);

writetable(all_info, hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Theta/Theta_Spectral_Analysis_Table_all_allevents.csv'));



%theta DLPFC
%pull out DLPFC
DLPFC_AvgEventNumber = AvgEventNumber(6,:)';
DLPFC_AvgEventDuration = AvgEventDuration(6,:)';
DLPFC_AvgPower = AvgPower(6,:)';
DLPFC_peakPower = peakPower(6,:)';

DLPFC_T = table(idvalues, DLPFC_AvgEventNumber, DLPFC_AvgEventDuration, DLPFC_AvgPower, DLPFC_peakPower);
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Theta_Spectral_Analysis_Table_allevents.mat'), 'T')


all_info = innerjoin(DLPFC_T, all_ages);


all_info = all_info(all_info.DLPFC_AvgEventNumber ~= 0, :);

writetable(all_info, hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Theta_Spectral_Analysis_Table_all_allevents.csv'));


%
%
%% Beta Analysis
clear x
clear Subjects
clear avgData

numSubj = length(idvalues);
x = cell(1,numSubj);
classLabels = cell(1,numSubj);

[num,txt,raw] = xlsread(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_allVariables.csv'));
alreadyRun = txt(:,1);
save(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/alreadyRun.mat'), 'alreadyRun')

%for terminal:
% load(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/alreadyRun.mat'))

for i = 1 : numSubj
    
    inputfile = setfiles{i};
    subject = inputfile(84:97);
    
    if any(strcmp(alreadyRun, subject))
        warning('%s already complete', subject)
        continue
        
    else
        
        %     clear events
        %     clear data
        %     clear hdr
        %     clear avgData
        %     clear classLabels
        %
        hdr = ft_read_header(inputfile);
        data = ft_read_data(inputfile, 'header', hdr);
        events = ft_read_event(inputfile, 'header', hdr);
        
        cfg = [];
        cfg.dataset = inputfile;
        cfg.headerfile = inputfile;
        cfg.channel = {'all', '-POz'};
        cfg.trialdef.eventtype = 'trigger';
        cfg.trialdef.eventvalue = '2';
        cfg.trialdef.prestim = 0;
        cfg.trialdef.poststim = 1;
        cfg.trialfun = 'ft_trialfun_general';
        cfg.trialdef.ntrials = 95;
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
            continue;
        end
        
        
        [data] = ft_preprocessing(cfg);
        
        
        %redefine the trails to be between 1-2 seconds of the delay period
        cfg.trl = [];
        cfg = rmfield(cfg, 'trl');
        cfg.toilim = [0 1];
        data = ft_redefinetrial(cfg, data);
        
        for k = 1:size(data.label,1)
            for j = 1:length(data.trial)
                trialData = (data.trial{1,j});
                channel = trialData(k,:);
                avgData(j,:) = channel;
            end
            
            x{i} = avgData';
            Subjects{i} = subject;
            
        end
    end
    
end

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

save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_SpecEvents_newSubs_FIX.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_TFR_newSubs_FIX.mat'), 'TFRs')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Fix/Beta_newSubs_FIX.mat'), 'Subjects')



%
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/BetaAvgEventNumber_allevents.mat'), 'AvgEventNumber')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/BetaAvgEventDuration_allevents.mat'), 'AvgEventDuration')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/BetaAvgPower_allevents.mat'), 'AvgPower')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/BetaPeakPower_allevents.mat'), 'peakPower')
% save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/BetaSDPower_allevents.mat'), 'SDPower')


AgeFile = xlsread(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/ERPs/EdaDI.xlsx'));

addpath(('Z:/DB_SQL/'));

all_ages = db_query('select id ||''_''|| to_char(vtimestamp,''YYYYMMDD'') as IDvalues, age, visitno from visit natural join enroll where etype like ''LunaID''  and vtype ilike ''eeg''  order by id, age ');
all_ages.idvalues= char(all_ages.idvalues);

for j = 1 : length(setfiles0)
    
    idvalues(j,:) = (setfiles0(j).name(1:14));
end


MeanAvgEventNumber = mean(AvgEventNumber,1)';
MeanAvgEventDuration = mean(AvgEventDuration,1)';
MeanAvgPower = mean(AvgPower, 1)';
MeanpeakPower = mean(peakPower, 1)';

T = table(idvalues, MeanAvgEventNumber, MeanAvgEventDuration, MeanAvgPower, MeanpeakPower);
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_allevents.mat'), 'T')

all_info = innerjoin(T, all_ages);



all_info = all_info(all_info.MeanAvgEventNumber ~= 0, :);

writetable(all_info, hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Beta/Beta_Spectral_Analysis_Table_all_allevents.csv'));
%


%Beta DLPFC
%pull out DLPFC
DLPFC_AvgEventNumber = AvgEventNumber(6,:)';
DLPFC_AvgEventDuration = AvgEventDuration(6,:)';
DLPFC_AvgPower = AvgPower(6,:)';
DLPFC_peakPower = peakPower(6,:)';

DLPFC_T = table(idvalues, DLPFC_AvgEventNumber, DLPFC_AvgEventDuration, DLPFC_AvgPower, DLPFC_peakPower);
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Beta_Spectral_Analysis_Table_allevents.mat'), 'T')


all_info = innerjoin(DLPFC_T, all_ages);


all_info = all_info(all_info.DLPFC_AvgEventNumber ~= 0, :);

writetable(all_info, hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/DLPFC/DLPFC_Beta_Spectral_Analysis_Table_all_allevents.csv'));




SE_Subs_Channel = cell(length(setfiles), 63);

%% Alpha Analysis

numSubj = length(idvalues);
x = cell(1,numSubj);
classLabels = cell(1,numSubj);

clear x
clear Subjects
clear channel
clear avgData

[num,txt,raw] = xlsread(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Alpha_allVariables.csv'));
alreadyRun = txt(:,1);
save(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/alreadyRun.mat'), 'alreadyRun')

% clear AvgEventDuration
% clear AvgEventFSpan
% clear AvgEventNumber
% clear AvgPower
% clear peakPower

for i = 1 : numSubj
    inputfile = setfiles{i};
    subject = inputfile(84:97);
    
    if any(strcmp(alreadyRun, subject))
        warning('%s already complete', subject)
        continue
        
    else
        
        %     clear events
        %     clear data
        %     clear hdr
        %     clear avgData
        %     clear classLabels
        
        hdr = ft_read_header(inputfile);
        data = ft_read_data(inputfile, 'header', hdr);
        events = ft_read_event(inputfile, 'header', hdr);
        
        cfg = [];
        cfg.dataset = inputfile;
        cfg.headerfile = inputfile;
        cfg.channel =   {'all', '-POz'};
        cfg.trialdef.eventtype = 'trigger';
        cfg.trialdef.eventvalue = '2';
        cfg.trialdef.prestim = 0;
        cfg.trialdef.poststim = 1;
        cfg.trialfun = 'ft_trialfun_general';
        cfg.trialdef.ntrials = 95;
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
            continue;
        end
        
        
        [data] = ft_preprocessing(cfg);
        
        
        %redefine the trails to be between 1-2 seconds of the delay period
        cfg.trl = [];
        cfg = rmfield(cfg, 'trl');
        cfg.toilim = [0 1];
        data = ft_redefinetrial(cfg, data);
        
        for k = 1:size(data.label,1)
            for j = 1:length(data.trial)
                trialData = (data.trial{1,j});
                channel = trialData(k,:);
                avgData(j,:) = channel;
            end
            
            x{i} = avgData';
            Subjects{i} = subject;
            
        end
        
    end
    
end

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


save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_SpecEvents_newSubs_FIX.mat'), 'specEv_struct')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_TFR_newSubs_FIX.mat'), 'TFRs')
save(hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Alpha/Fix/Alpha_newSubs_FIX.mat'), 'Subjects')




