function [errorSubjects] = evokedInducedActivity(i, inputfile, triggerValue, outpath, subject,errorSubjects)

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset(inputfile); % load in eeg file

if EEG.srate ~= 512 %check that the subject has been resampled to 512 Hz, and if not
    EEG = pop_resample(EEG, 512, 0.8, 0.4);

end

[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );

EEG = pop_rmdat( EEG, {triggerValue},[-0.2 0.5] ,0); % select events with trigger value you want
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
EEG = eeg_checkset( EEG );

try
    EEG = pop_epoch( EEG, {triggerValue}, [-0.2 0.5], 'epochinfo', 'yes'); % create epochs using selected events
catch
    disp(['subject doesnt have correct epoch: ' inputfile(112:125)])
    errorSubjects{i} = [inputfile(112:125)];
    ersp = [];
    itc = [];
    powbase = [];
    times = [];
    freqs = [];
    tfdata = [];
    return;
end

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
EEG = eeg_checkset( EEG );

try
    EEG = pop_rmbase( EEG, [-.2 0] ,[]); % remove baseline
catch
    disp(['subject doesnt have correct baseline values: ' inputfile(112:125)])
    errorSubjects{i} = [inputfile(112:125)];
    ersp = [];
    itc = [];
    powbase = [];
    times = [];
    freqs = [];
    tfdata = [];
    return;
end

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off');
EEG = eeg_checkset( EEG );


% compute the time frequency spectrum
% [-200  499]: time bin to analyze , [2  15]: 2 wavelets at low frequencies and 15 at high frequencies
% baseline 0: use entire section before 0 for the baseline and use divisive strategy for baseline correction
% freqs, [10 70]: look at frequencies 10 - 70 Hz
% padratio 16: increases frequency resolution
% powbase is the baseline power spectrum

% tfdata returns the time frequency array for each trial

parpool('local', 6);
parfor c = 1:64
    [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = pop_newtimef(EEG, 1, c, [-200 499], [2 15], 'baseline',[0], 'freqs', [10 70], 'topovec', c, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'plotphase', 'off', 'padratio', 16);

    channelERSP{c} = ersp;
    channelITC{c} = itc;
    channelPowbase{c} = powbase;
    channeltfdata{c} = tfdata;
end
delete(gcp('nocreate')); %stop the current parallel pool before you start a new one


% ERSP
for j = 1:length(channelERSP)
    channel = channelERSP{1,j};
    avgChannel = mean(mean(channel));
    subjectAvgERSP{j} = avgChannel;
end

% ITC
for j = 1:length(channelITC)
    channel = channelITC{1,j};
    avgChannel = abs(mean(mean(channel)));
    subjectAvgITC{j} = avgChannel;
end

% baseline  power
for j = 1:length(channelPowbase)
    channel = channelPowbase{1,j};
    avgChannel = mean(channel,2);
    subjectAvgPowBase{j} = avgChannel;
end

% TF
for j = 1:length(channeltfdata)
    channel = channeltfdata{1,j};
    avgChannel = mean(squeeze(abs(mean(mean(channel))))); % average across time and frequency for each trial and then average across trials
    subjectAvgTF{j} = avgChannel;
end

%% create a table with all the info

T = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double','double','double','double'});
T.Properties.VariableNames = {'Subject', 'Channel', 'ERSP', 'ITC','BaselinePower','Induced'};

for c = 1:64
    T.ERSP(c) = subjectAvgERSP{c};
    T.ITC(c) = subjectAvgITC{c};
    T.BaselinePower(c) = subjectAvgPowBase{c};
    T.Induced(c) = subjectAvgTF{c};
    T.Channel(c) = c;
end
T.Subject(1:64) = subject;

savePath = [outpath '/' subject '_SNRdata.csv'];
writetable(T, savePath);

errorSubjects{i} = [];
eeglab redraw;
close all;

