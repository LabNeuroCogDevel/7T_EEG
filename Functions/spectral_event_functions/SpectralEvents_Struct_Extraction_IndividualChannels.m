
%% Extract our data from the large struct
addpath(genpath('Functions/spectral_event_functions/'));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/')))
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20191213'))
ft_defaults

% values to initalize
task = 'Resting_State'; % which task do you want to run. CHANGE ACCORDING TO MGS, REST, SNR, etc. *****
epoch = 'Delay'; % which epoch do you want to run, if REST there is no epoch
band = 'Gamma';
seconds = '3_4';
region = 'selectiveFrontal';

% In path
inpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data');

% Outpath
outpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/');


if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    specEv_struct = load([savefolder band '_specEV_' region '_Delay_' seconds '.mat']);
    specEv_struct = specEv_struct.specEv_struct;

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    specEv_struct = load([savefolder band '_specEV_' region '_Fix.mat']);
    specEv_struct = specEv_struct.specEv_struct;

elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    specEv_struct = load([savefolder band '_specEV_' region '_RS.mat']);
    specEv_struct = specEv_struct.specEv_struct;

end


largeT = table('Size', [0, 6], 'VariableTypes', {'string', 'double','double', 'double','double','double'});
largeT.Properties.VariableNames = {'Subject','Channel', 'Trial','Band',[band 'Trial_Power'],[band 'Event_Number'],[band 'Event_Duration']};

for j = 1:length(specEV_struct_allChannels)
    specEv_struct = specEV_struct_allChannels{1,j};

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


    for i = 1:length(specEv_struct)

        AvgPower_PerSubject(1,i) = mean(Power{1,i});
        AvgEventNumber_PerSubject(1,i) = mean(EventNumber{1,i});
        AvgEventDuration_PerSubject(1,i) = mean(EventDuration{1,i});
        AvgEventMaxPower_PerSubject(1,i) = mean(EventMaxPower{1,i});
        AvgEventMaxFreq_PerSubject(1,i) = mean(EventMaxFreq{1,i});
        SDPower_PerSubject(1,i) = std(Power{1,i});

    end

    %Extract info trial by trial
    for i = 1:length(specEv_struct) %Subjects
        clear Trial_power
        clear Trial_eventDuration
        clear Trial_eventNumber

        T = table('Size', [length(Power{1,i}), 6], 'VariableTypes', {'string', 'double','double', 'double','double','double'});
        T.Properties.VariableNames = {'Subject','Channel', 'Trial' ,'Band',[band 'Trial_Power'],[band 'Event_Number'],[band 'Event_Duration']};

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
        T.Subject(1:length(T.Trial_Power)) = SubjectsRan(:,i);
        T.Channel = j+zeros(size(Subject_power,1),1);
        T.Band(1:length(T.Trial_Power)) = band;

        largeT = [largeT; T];

    end
end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    writetable(largeT, [savefolder band '_allChannels_TrialLevel_Delay' seconds '.csv']);

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    writetable(largeT, [savefolder band '_allChannels_TrialLevel_Fix.csv']);

elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    writetable(largeT, [savefolder band '_allChannels_TrialLevel_RS.csv']);

end


%% peak frequency and power

largeT = table('Size', [0, 6], 'VariableTypes', {'string', 'double','double', 'double','double','double'});
largeT.Properties.VariableNames = {'Subject','Channel', 'Trial', 'Band',[band 'Trial_Power'],[band 'Event_Number'],[band 'Event_Duration']};

% Extract info trial by trial
for i = 1:length(specEv_struct) %Subjects


    T = table('Size', [length(EventMaxFreq{1,i}), 4], 'VariableTypes', {'string', 'double', 'double','double'});
    T.Properties.VariableNames = {'Subject', 'Trial', 'Band', [band 'Peak_Frequency'], [band 'Peak_Power']};

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
    T.Band(1:length(T.Peak_Frequency)) = band;
    T.Subject(1:length(T.Peak_Frequency)) = Subjects(i);
    largeT = [largeT; T];

end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/' epoch seconds '/'];
    writetable(largeT, [savefolder band '_allChannels_TrialLevel_PeakFreq_Delay' seconds '.csv']);

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/' epoch ];
    writetable(largeT, [savefolder band '_allChannels_TrialLevel_PeakFreq_Fix.csv']);

elseif task == 'Resting_State'
    savefolder = [outpath band '/' task '/'];
    writetable(largeT, [savefolder band '_allChannels_TrialLevel_PeakFreq_RS.csv']);

end





