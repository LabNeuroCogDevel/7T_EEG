
%% Extract our data from the large struct

addpath(('/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))
addpath((genpath('Functions/spectral_event_functions/')));
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20191127'))
ft_defaults

% values to initalize
task = 'Resting_State'; % which task do you want to run. CHANGE ACCORDING TO MGS, REST, SNR, etc. *****
epoch = 'Delay'; % which epoch do you want to run, if REST there is no epoch
band = 'Gamma';
seconds = '3_4';


if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    specEv_struct = load([savefolder band '_specEV_DLPFC_Delay_' seconds '.mat']);
    specEv_struct = specEv_struct.specEv_struct; 

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    specEv_struct = load([savefolder band '_specEV_DLPFC_Fix.mat']);
    specEv_struct = specEv_struct.specEv_struct; 
    
elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    specEv_struct = load([savefolder band '_specEV_DLPFC_RS.mat']);
    specEv_struct = specEv_struct.specEv_struct; 

end


EventNumber = cell(1, length(specEv_struct));
EventDuration = cell(1, length(specEv_struct));
Power = cell(1, length(specEv_struct));
EventMaxPower = cell(1, length(specEv_struct));
EventMaxFreq = cell(1, length(specEv_struct));

for i = 1:length(specEv_struct)
    
    subject = specEv_struct(i);
    if isempty(subject.TrialSummary)
        continue
    end
    
    EventNumber{1,i} = subject.TrialSummary.TrialSummary.eventnumber;
    EventDuration{1,i} = subject.TrialSummary.TrialSummary.meaneventduration;
    Power{1,i} = subject.TrialSummary.TrialSummary.meaneventpower;
    EventMaxPower{1,i} = subject.Events.Events.maximapower;
    EventMaxFreq{1,i} = subject.Events.Events.maximafreq;
end


if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    save([savefolder band '_EventNumber' seconds '.mat'], 'EventNumber')
    save([savefolder band '_EventDuration' seconds '.mat'], 'EventDuration')
    save([savefolder band '_Power' seconds '.mat'], 'Power')
    save([savefolder band '_EventMaxPower' seconds '.mat'], 'EventMaxPower')
    save([savefolder band '_EventMaxFreq' seconds '.mat'], 'EventMaxFreq')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    save([savefolder band '_EventNumber' epoch '.mat'], 'EventNumber')
    save([savefolder band '_EventDuration' epoch '.mat'], 'EventDuration')
    save([savefolder band '_Power' epoch '.mat'], 'Power')
    save([savefolder band '_EventMaxPower' epoch '.mat'], 'EventMaxPower')
    save([savefolder band '_EventMaxFreq' epoch '.mat'], 'EventMaxFreq')
    
elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    save([savefolder band '_EventNumberRS.mat'], 'EventNumber')
    save([savefolder band '_EventDurationRS.mat'], 'EventDuration')
    save([savefolder band '_PowerRS.mat'], 'Power')
    save([savefolder band '_EventMaxPowerRS.mat'], 'EventMaxPower')
    save([savefolder band '_EventMaxFreqRS.mat'], 'EventMaxFreq')

end


for i = 1:length(specEv_struct)
    
    AvgPower_PerSubject(1,i) = mean(Power{1,i});
    AvgEventNumber_PerSubject(1,i) = mean(EventNumber{1,i});
    AvgEventDuration_PerSubject(1,i) = mean(EventDuration{1,i});
    AvgEventMaxPower_PerSubject(1,i) = mean(EventMaxPower{1,i});
    AvgEventMaxFreq_PerSubject(1,i) = mean(EventMaxFreq{1,i});
    SDPower_PerSubject(1,i) = std(Power{1,i});
    
end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    save([savefolder band '_AvgEventNumber' seconds '.mat'], 'AvgEventNumber_PerSubject')
    save([savefolder band '_AvgEventDuration' seconds '.mat'], 'AvgEventDuration_PerSubject')
    save([savefolder band '_AvgPower' seconds '.mat'], 'AvgPower_PerSubject')
    save([savefolder band '_AvgEventMaxPower' seconds '.mat'], 'AvgEventMaxPower_PerSubject')
    save([savefolder band '_AvgEventMaxFreq' seconds '.mat'], 'AvgEventMaxFreq_PerSubject')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    save([savefolder band '_AvgEventNumber' epoch '.mat'], 'AvgEventNumber_PerSubject')
    save([savefolder band '_AvgEventDuration' epoch '.mat'], 'AvgEventDuration_PerSubject')
    save([savefolder band '_AvgPower' epoch '.mat'], 'AvgPower_PerSubject')
    save([savefolder band '_AvgEventMaxPower' epoch '.mat'], 'AvgEventMaxPower_PerSubject')
    save([savefolder band '_AvgEventMaxFreq' epoch '.mat'], 'AvgEventMaxFreq_PerSubject')
    
elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    save([savefolder band '_AvgEventNumberRS.mat'], 'AvgEventNumber_PerSubject')
    save([savefolder band '_AvgEventDurationRS.mat'], 'AvgEventDuration_PerSubject')
    save([savefolder band '_AvgPowerRS.mat'], 'AvgPower_PerSubject')
    save([savefolder band '_AvgEventMaxPowerRS.mat'], 'AvgEventMaxPower_PerSubject')
    save([savefolder band '_AvgEventMaxFreqRS.mat'], 'AvgEventMaxFreq_PerSubject')

end

largeT = table('Size', [0, 5], 'VariableTypes', {'string', 'double', 'double','double','double'});
largeT.Properties.VariableNames = {'Subject', 'Trial' ,[band '_Trial_Power'],[band '_Event_Number'], [band '_Event_Duration']};

%Extract info trial by trial
for i = 1:length(specEv_struct) %Subjects
    clear Trial_power
    clear Trial_eventDuration
    clear Trial_eventNumber
    
    T = table('Size', [length(Power{1,i}), 5], 'VariableTypes', {'string', 'double', 'double','double','double'});
    T.Properties.VariableNames = {'Subject', 'Trial', 'Band',[band '_Trial_Power'], [band '_Event_Number'], [band '_Event_Duration']};
    
    if isempty(Power{1,i}) || isempty(EventNumber{1,i}) || isempty(EventDuration{1,i})
        continue
    end
    

    Subject_power = Power{1,i};

    Subject_channel_eventNumber = EventNumber{1,i};

    Subject_channel_eventDuration= EventDuration{1,i};


    if isempty(Power{1,i}) || isempty(EventNumber{1,i}) || isempty(EventDuration{1,i})
        continue
    end
    
    
    T.Trial_Power = Subject_power;
    T.Event_Number = Subject_channel_eventNumber;
    T.Event_Duration = Subject_channel_eventDuration;
    T.Trial = [1:length(T.Trial_Power)]';
    T.Subject(1:length(T.Trial_Power)) = idvalues(i,:);
    T.Band(1:length(T.Trial_Power)) = band;
    largeT = [largeT; T];
    
    
end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    writetable(largeT, [savefolder band '_subs_' seconds '.csv']);

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    writetable(largeT, [savefolder band '_subs_Fix.csv']);

elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    writetable(largeT, [savefolder band '_subs_RS.csv']);

end


%% peak frequency and power

largeT = table('Size', [0, 4], 'VariableTypes', {'string', 'double', 'double','double'});
largeT.Properties.VariableNames = {'Subject', 'Trial', 'Band', [band 'Peak_Frequency'], [band 'Peak_Power']};

% %Extract info trial by trial
for i = 1:length(specEv_struct) %Subjects
    
    
    T = table('Size', [length(EventMaxFreq{1,i}), 4], 'VariableTypes', {'string', 'double', 'double','double'});
    T.Properties.VariableNames = {'Subject', 'Trial', [band 'Peak_Frequency'], [band 'Peak_Power']};
    
    if isempty(EventMaxFreq{1,i}) || isempty(EventMaxPower{1,i})
        continue
    end
    
    
    Subject_peakFrequency = EventMaxFreq{1,i};
    
    Subject_peakPower = EventMaxPower{1,i};
    
    
    
    if isempty( EventMaxFreq(1,i)) || isempty(EventMaxPower(1,i))
        continue
    end
    
    
    T.Peak_Frequency = Subject_peakFrequency;
    T.Peak_Power = Subject_peakPower;
    T.Trial = [1:length(T.Peak_Frequency)]';
    T.Subject(1:length(T.Peak_Frequency)) = Subjects(i);
    T.Band(1:length(T.Peak_Frequency)) = band;
    largeT = [largeT; T];
    
end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    writetable(largeT, [savefolder band '_subs_PeakFreq_Power' seconds '.csv']);

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    writetable(largeT, [savefolder band '_subs_PeakFreq_Power_Fix.csv']);

elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    writetable(largeT, [savefolder band '_subs_PeakFreq_Power_RS.csv']);

end




