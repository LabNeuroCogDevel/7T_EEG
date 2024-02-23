function [] = ITC_allChans(triggerValue)

%% set needed paths
cd '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions'
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Functions/SNR_processing')));
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions/preprocessing'));

%% set initial values
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize');
% triggerValue = '4';
% channelValues = [4,5,6,36,37,38]; if you only want to run DLPFC
outpath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/ITC_indivSubs_allChans');

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

    if triggerValue == '2'
        savePath = [outpath '/' subject '_ITC_20Hz.csv'];
    elseif triggerValue =='3'
        savePath = [outpath '/' subject '_ITC_30Hz.csv'];
    elseif triggerValue == '4'
        savePath = [outpath '/' subject '_ITC_40Hz.csv'];
    end

    if ~exist(savePath, 'file')
        fprintf('The file %s does not exist, running SNR code\n', savePath)

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

        for c = 1:64
            [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = pop_newtimef( EEG, 1, c, [-199  799], [2  15] , 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'baseline',[0], 'freqs', [10 70], 'plotphase', 'off', 'padratio', 4, 'winsize', 102);
            subjectITC(:,:,c) = itc;
            close all;
        end







        allchans = [];
        % create a giant array for all subs and all TF values
        for c = 1:64
            chanTFarray = abs(subjectITC(:,:,c));
            % Split the char variable into two parts using '_'

            [timeGrid, freqGrid] = meshgrid(times, freqs);
            % Reshape the matrix and grids into a 3-column array
            subDataArray = [timeGrid(:), freqGrid(:), chanTFarray(:), repmat(c, 9600, 1)];

            allchans = [allchans; subDataArray];
        end


        T = table('Size', [0, 5], 'VariableTypes', {'string', 'double', 'double','double','double'});
        T.Properties.VariableNames = {'Subject', 'Channel','time','freq','ITC'};

        T.Subject(1:length(allchans)) = idvalues{i};
        T.Channel = allchans(:,4);
        T.freq = allchans(:,2);
        T.time = allchans(:,1);
        T.ITC = allchans(:,3);

        if triggerValue == '2'

            savePath = [outpath '/' subject '_ITC_20Hz.csv'];
            writetable(T, savePath);


        elseif triggerValue =='3'
            savePath = [outpath '/' subject '_ITC_30Hz.csv'];
            writetable(T, savePath);

        elseif triggerValue == '4'
            savePath = [outpath '/' subject '_ITC_40Hz.csv'];
            writetable(T, savePath);

        end
    end

end
end
