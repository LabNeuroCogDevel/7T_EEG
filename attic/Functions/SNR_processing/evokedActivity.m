function [errorSubjects, channelPower] = evokedActivity(i, inputfile, triggerValue, outpath, subject,errorSubjects)

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
try
    EEG = pop_loadset(inputfile); % load in eeg file
catch
    disp(['subject missing fdt: ' subject])
    errorSubjects{i} = subject;
    channelPower = [];
    return;
end


% if max([EEG.event.type]) ~= 4
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
        disp(['subject doesnt have correct epoch: ' subject]);
        errorSubjects{i} = subject;
        channelPower = [];

        return;
    end
end

[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );

EEG = pop_rmdat( EEG, {triggerValue},[0 0.8] ,0); % select events with trigger value you want
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
EEG = eeg_checkset( EEG );

try
    EEG = pop_epoch( EEG, {triggerValue}, [0 0.8], 'epochinfo', 'yes'); % create epochs using selected events
catch
    disp(['subject doesnt have correct epoch: ' subject])
    errorSubjects{i} = subject;
    channelPower = [];
   
    return;
end

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
EEG = eeg_checkset( EEG );

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

%remove the prestim data now that its been baselined
% try
%     EEG = pop_epoch( EEG, {triggerValue}, [0 0.5], 'epochinfo', 'yes'); % create epochs using selected events
% catch
%     disp(['subject doesnt have correct epoch: ' subject])
%     errorSubjects{i} = subject;
%     ersp = [];
%     itc = [];
%     powbase = [];
%     times = [];
%     freqs = [];
%     tfdata = [];
%     return;
% end

averagedData = mean(EEG.data,3); %average across trials. left with channels x time points
hz = linspace(0,EEG.srate/2, floor(EEG.pnts/2)+1);

for c = 1:64
    channelPower_allFreqs(c,:) = abs(fft(averagedData(c,:))).^2;
    
    idx = find(hz>=36 & hz<= 45);
    channelPower(c,:) = mean(channelPower_allFreqs(c,idx),2);
    
end


