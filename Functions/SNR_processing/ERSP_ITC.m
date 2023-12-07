function [errorSubjects, subjectAvgERSP, subjectAvgITC, subjectAvgPowBase, subjectAvgTF] = ERSP_ITC(i, inputfile, triggerValue, outpath, subject,errorSubjects)


% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
try
    EEG = pop_loadset(inputfile); % load in eeg file
catch
    disp(['subject missing fdt: ' subject])
    errorSubjects{i} = subject;
    subjectAvgERSP = [];
    subjectAvgITC = [];
    subjectAvgPowBase = [];
    subjectAvgTF = [];
    return;
end

% if max([EEG.event.type]) ~= 4 || max([EEG.event.type]) ~= 121
%     errorSubjects{i} = subject;
%     ersp = [];
%     itc = [];
%     powbase = [];
%     times = [];
%     freqs = [];
%     tfdata = [];
%     return;
% end

if EEG.srate ~= 512 %check that the subject has been resampled to 512 Hz, and if not
    try
        EEG = pop_resample(EEG, 512, 0.8, 0.4);
    catch
        disp(['subject doesnt have correct epoch: ' subject])
        errorSubjects{i} = subject;
        subjectAvgERSP = [];
        subjectAvgITC = [];
        subjectAvgPowBase = [];
        subjectAvgTF = [];
        return;
    end
end

[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );

EEG = pop_rmdat( EEG, {triggerValue},[-0.2 0.8] ,0); % select events with trigger value you want
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
EEG = eeg_checkset( EEG );

try
    EEG = pop_epoch( EEG, {triggerValue}, [-0.2 0.8], 'epochinfo', 'yes'); % create epochs using selected events
catch
    disp(['subject doesnt have correct epoch: ' subject])
    errorSubjects{i} = subject;
    subjectAvgERSP = [];
    subjectAvgITC = [];
    subjectAvgPowBase = [];
    subjectAvgTF = [];
    return;
end

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
EEG = eeg_checkset( EEG );
%
% try
%     EEG = pop_rmbase( EEG, [-.2 0] ,[]); % remove baseline
% catch
%     disp(['subject doesnt have correct baseline values: ' subject])
%     errorSubjects{i} = subject;
%     ersp = [];
%     itc = [];
%     powbase = [];
%     times = [];
%     freqs = [];
%     tfdata = [];
%     return;
% end
%
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off');
% EEG = eeg_checkset( EEG );


% compute the time frequency spectrum
% [-200  499]: time bin to analyze , [2  15]: 2 wavelets at low frequencies and 15 at high frequencies
% baseline 0: use entire section before 0 for the baseline and use divisive strategy for baseline correction
% freqs, [36 45]: look at frequencies 10 - 70 Hz
% padratio 16: increases frequency resolution
% powbase is the baseline power spectrum

% tfdata returns the time frequency array for each trial

parpool('local', 4);
parfor c = 1:64
    figure;
    [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = pop_newtimef( EEG, 1, c, [-199  799], [2  15] , 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'baseline',[0], 'freqs', [10 70], 'plotphase', 'off', 'padratio', 4, 'winsize', 102);
    channelERSP{c} = ersp;
    channelITC{c} = itc;
    channelPowbase{c} = powbase;
    channeltfdata{c} = tfdata;
    channelfreqs{c} = freqs;
end
delete(gcp('nocreate')); %stop the current parallel pool before you start a new one

idx = find(channelfreqs{1,1}>=36 & channelfreqs{1,1}<= 45); % only going to look at the values between the freq of interest

% ERSP
for j = 1:length(channelERSP)
    channel = channelERSP{1,j};
    selectfreq = channel(idx,:);
    avgChannel = mean(mean(selectfreq));
    subjectAvgERSP{j} = avgChannel;
end

% ITC
for j = 1:length(channelITC)
    channel = channelITC{1,j};
    selectfreq = channel(idx,:);
    avgChannel = mean(mean(abs(selectfreq)));
    subjectAvgITC{j} = avgChannel;
end

% baseline  power
for j = 1:length(channelPowbase)
    channel = channelPowbase{1,j};
    selectfreq = channel(:,idx);
    avgChannel = mean(selectfreq,2);
    subjectAvgPowBase{j} = avgChannel;
end

% TF
for j = 1:length(channeltfdata)
    channel = channeltfdata{1,j};
    selectfreq = channel(idx,:);
    avgChannel = mean(squeeze(mean(mean(abs(selectfreq))))); % average across time and frequency for each trial and then average across trials
    subjectAvgTF{j} = avgChannel;
end


errorSubjects{i} = [];
eeglab redraw;
close all;

