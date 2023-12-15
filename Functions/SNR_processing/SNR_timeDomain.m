
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
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
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
    timeidx = find(time < 0);
    evokedtime = find(time > 50 & time < 500);

    for c = 51
        chanData = EEG.data(c,:,:);
        baseline = squeeze(chanData(:,timeidx,:));
        trialdata = squeeze(chanData(:,evokedtime,:));
        avgTrialdata = mean(trialdata, 2);
        avgBaseline = mean(mean(baseline, 2));

        for t = 1:length(avgTrialdata)
            baselineRemoved(t) = avgTrialdata(t,1) - avgBaseline;

        end

        evoked = (baselineRemoved.^2);
        noise = avgTrialdata.^2;
        SNR = 10*log10(evoked./noise');
        avgSNR(i) = mean(SNR);

        trialVar = var(trialdata,0,1);

%         avgTrialVariance = mean(trialVar);
%         sumEvoked = sum(evoked);
%         avgEvoked = mean(evoked);
%         varianceBasedSNR(i) = 10*log10(avgEvoked/(avgTrialVariance));


        
        % to plot an average SNR using windows 
        % % Parameters
        % window_size = 30;
        % overlap = 0.5;
        % 
        % % Calculate step size and number of windows
        % step_size = round(window_size * (1 - overlap));
        % num_windows = floor((length(SNR) - window_size) / step_size) + 1;
        % 
        % % Initialize array for averaged values
        % averages = zeros(1, num_windows);
        % 
        % % Calculate averages for each window
        % for i = 1:num_windows
        %     start_index = (i - 1) * step_size + 1;
        %     end_index = start_index + window_size - 1;
        % 
        %     window = SNR(start_index:end_index);
        %     averages(i) = mean(window);
        % 
        %      timewindow = evokedtime(start_index:end_index);
        %     timeaverages(i) = mean(timewindow);
        % end

    end

end

if triggerValue == '2'
    
    savePath = [outpath '/' 'SNR_20Hz.mat'];
    save(savePath, 'avgSNR');
    
elseif triggerValue =='3'
    
    savePath = [outpath '/' 'SNR_30Hz.mat'];
    save(savePath, 'avgSNR');
    
elseif triggerValue =='4'
    
    savePath = [outpath '/' 'SNR_40Hz.mat'];
    save(savePath, 'avgSNR');
    writematrix(varianceBasedSNR, [outpath '/SNR_40Hz.csv']);
end


 T = table('Size', [0, 2], 'VariableTypes', {'string', 'double'});
 T.Properties.VariableNames = {'Subject', 'SNR'};
 
 for i = 1:length(idvalues)
     T.Subject(i) = idvalues{i};
     T.SNR(i) = avgSNR(:,i); 
 end
 
 savePath = [outpath '/SNRdata_40Hz_method2.csv'];
 writetable(T, savePath)

