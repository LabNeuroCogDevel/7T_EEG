%% set needed paths
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Functions/SNR_processing')));
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions/preprocessing'));

%% set initial values
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize');
triggerValue = '4';
% channelValues = [4,5,6,36,37,38]; if you only want to run DLPFC
outpath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/');

%% load in all the data files
setfiles0 = dir([datapath,'/*icapru*.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); 
    % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    idvalues{j} = (setfiles0(j).name(1:14));
end

idvalues = unique(idvalues);
numSubj = length(idvalues);
errorSubjects = cell(1,numSubj);

%% run evoked and induced activity function
for i = 1:numSubj
    disp(i);
    subject = idvalues{i};
    inputfile = setfiles{i};

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
        continue;
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
            continue;
        end
    end
    
try
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );

    EEG = pop_rmdat( EEG, {triggerValue},[-0.2 0.8] ,0); % select events with trigger value you want
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
    EEG = eeg_checkset( EEG );
catch
       disp(['subject wont run through pop_rmdat: ' subject])
        errorSubjects{i} = subject;
        subjectAvgERSP = [];
        subjectAvgITC = [];
        subjectAvgPowBase = [];
        subjectAvgTF = [];
        continue;
        
end

    try
        EEG = pop_epoch( EEG, {triggerValue}, [-0.2 0.8], 'epochinfo', 'yes'); % create epochs using selected events
    catch
        disp(['subject doesnt have correct epoch: ' subject])
        errorSubjects{i} = subject;
        subjectAvgERSP = [];
        subjectAvgITC = [];
        subjectAvgPowBase = [];
        subjectAvgTF = [];
        continue;
    end

    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
    EEG = eeg_checkset( EEG );

    time = EEG.times;
    baselineIdx = find(time > -199 & time < 0);
    baselinetime = time(:,baselineIdx);
    
    for c = 51
        chanData = squeeze(EEG.data(c,baselineIdx,:));
        avgData = mean(chanData,2); %avg all the epochs in the time domain then take fft to find evoked
        % window_percentage = 10;
        % window_size = round(length(avgData) * (window_percentage / 100));
        % hanning_window = hanning(length(avgData));
        % windowed_signal = avgData .* hanning_window;
        
        % Perform FFT
        N = 2^nextpow2(length(avgData));   % Next power of 2 from signal length for FFT
        fft_result = fft(avgData, N)/length(avgData);
        evokedPower = (abs(fft_result(1:N/2+1)).^2); % Square the magnitude, turn to db
        
        frequencies = 512/2*linspace(0,1,N/2+1);
        
        % % Plot the original signal and the FFT result
        % subplot(2,1,1);
        % plot(baselinetime, avgData);
        % title('Original Signal');
        % xlabel('Time (s)');
        % ylabel('Amplitude');
        % 
        % subplot(2,1,2);
        % plot(frequencies, evokedPower);
        % xlim([0,70]);
        % title('FFT');
        % xlabel('Frequency (Hz)');
        % ylabel('Power');
        % 
        % to find power on each individual trial
        for t = 1:size(chanData, 2)
            trialData = chanData(:,t);
            % Apply a 10% Hanning window
            % window_percentage = 10;
            % window_size = round(length(trialData) * (window_percentage / 100));
            % hanning_window = hanning(length(trialData));
            % windowed_signal = trialData .* hanning_window;
            
            % Perform FFT
            N = 2^nextpow2(length(trialData));   % Next power of 2 from signal length for FFT
            fft_result = fft(trialData, N)/length(trialData);
            trialPower_spectrum(:,t) = (abs(fft_result(1:N/2+1)).^2); % Square the magnitude, turn to db

            frequencies = 512/2*linspace(0,1,N/2+1);
            
            % % Plot the original signal and the FFT result
            % subplot(2,1,1);
            % plot(evokedtime, mean(chanData,2));
            % title('Original Signal');
            % xlabel('Time (s)');
            % ylabel('Amplitude');
            % 
            % subplot(2,1,2);
            % plot(frequencies, mean(trialPower_spectrum,2));
            % xlim([0,70]);
            % title('FFT');
            % xlabel('Frequency (Hz)');
            % ylabel('Power');

        end
        
        totalPower = mean(trialPower_spectrum,2);
        inducedPower = totalPower - evokedPower;
        figure;
        plot(frequencies, inducedPower);
        xlim([0,70]);
        title('FFT ');
        xlabel('Frequency (Hz)');
        ylabel('Power');
        hold on; 
        plot(frequencies, totalPower, 'LineWidth',4, 'Color','r');

    end
    subjectTotalPower(i) = mean(totalPower(isfinite(totalPower))); %omit -Inf values that resulted from 10*log10(0)
    subjectEvokedPower(i) = mean(evokedPower(isfinite(evokedPower))); %omit -Inf values that resulted from 10*log10(0)
    subjectInducedPower(i) = mean(inducedPower(isfinite(inducedPower))); %omit -Inf values that resulted from 10*log10(0)


end



 T = table('Size', [0, 4], 'VariableTypes', {'string', 'double','double','double'});
 T.Properties.VariableNames = {'Subject', 'TotalPower', 'EvokedPower','InducedPower'};
 
 for i = 1:length(idvalues)
     T.Subject(i) = idvalues{i};
     T.TotalPower(i) = subjectTotalPower(i); 
     T.EvokedPower(i) = subjectEvokedPower(i); 
     T.InducedPower(i) = subjectInducedPower(i); 

 end
 
 savePath = [outpath '/SNRdata_40Hz_totalEvokedInduced_baseline.csv'];
 writetable(T, savePath)

