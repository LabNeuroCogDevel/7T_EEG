
%% set needed paths
cd '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions'
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Functions/SNR_processing')));
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions/preprocessing'));

%% set initial values
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize');
triggerValue = '3';
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





    % if triggerValue == '2'
    %
    %     savePath = [outpath '/' 'allsubjectsITCmatries_20Hz.mat'];
    %     save(savePath, 'subjectITC');
    %
    % elseif triggerValue =='3'
    %
    %     savePath = [outpath '/' 'allsubjectsITCmatries_30Hz.mat'];
    %     save(savePath, 'subjectITC');
    %
    % elseif triggerValue =='4'
    %
    %     savePath = [outpath '/' 'allsubjectsITCmatries_40Hz.mat'];
    %     save(savePath, 'subjectITC');
    % end
    %


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


            savePath = [outpath '/' subject '_ITC_40Hz.csv'];
            writetable(T, savePath);













    % 
    % % Get the size of the original array
    % [m, n, p] = size(subjectITC);
    % 
    % % Create a cell array to store the 2D arrays
    % cellArray = cell(p, 2);
    % 
    % % Loop through the 3rd dimension and extract 2D arrays into the cell array
    % for i = 1:p
    %     cellArray{i,2} = subjectITC(:, :, i);
    %     cellArray{i,1} = idvalues{1,i};
    % end
    % 
    % % Extract the common subject IDs
    % commonSubjectIDs = intersect(cellArray(:, 1), ages(:, 1));
    % 
    % % Initialize the merged cell array
    % mergedCellArray = [];
    % 
    % % Loop through common subject IDs and concatenate corresponding rows
    % for i = 1:length(commonSubjectIDs)
    %     subjectID = commonSubjectIDs{i};
    % 
    %     % Find rows with the common subject ID in each cell array
    %     idx1 = ismember(cellArray(:, 1), subjectID);
    %     idx2 = ismember(ages(:, 1), subjectID);
    % 
    %     % Concatenate corresponding rows
    %     mergedRow = [subjectID, {cellArray{idx1, 2}}, ages{idx2, 2}];
    % 
    %     % Append the merged row to the result
    %     mergedCellArray = [mergedCellArray; mergedRow];
    % end


    % for f = 1:size(subjectITC, 1)
    %     for t = 1:size(subjectITC,2)
    %         for s = 1:length(mergedCellArray)
    %             subTFarray = abs(mergedCellArray{s,2});
    %
    %             % Split the char variable into two parts using '_'
    %             splitValues = strsplit(mergedCellArray{s,1}, '_');
    %
    %             % Convert the split values to numeric
    %             timeFreqValue(s,1) = str2double(splitValues{1});
    %             timeFreqValue(s,2) = str2double(splitValues{2});
    %             timeFreqValue(s,3) = subTFarray(f,t);
    %             timeFreqValue(s,4) = mergedCellArray{s,3}; % array of everyones (f,t) freq time value and age
    %         end
    %
    %         % Find rows where value in first column is 1
    %         rowsWithZeroInColumn1 = find(timeFreqValue(:, 3) == 0);
    %
    %         % Remove rows with all values 0
    %         timeFreqValue(rowsWithZeroInColumn1, :) = [];
    %
    %         % Mean-center the age vector
    %         meanAge = mean(timeFreqValue(:,4));
    %         timeFreqValue(:,5) = timeFreqValue(:,4) - meanAge;
    %
    %         TFtable = array2table(timeFreqValue, 'VariableNames', {'lunaID','visitDate','TFvalue', 'age', 'meanCenteredAge'});
    %
    %         lme = fitlme(TFtable,'TFvalue~meanCenteredAge+(1|lunaID)');
    %
    % %        model = fitlm(timeFreqValue(:,2), timeFreqValue(:,1), 'linear');
    %        ageEstimate(f,t) = lme.Coefficients{2,2};
    %        ageTvalue(f,t) = lme.Coefficients{2,4};
    %        agePvalue(f,t) = lme.Coefficients{2,6};
    %
    %        intEstimate(f,t) = lme.Coefficients{1,2};
    %        intTvalue(f,t) = lme.Coefficients{1,4};
    %        intPvalue(f,t) = lme.Coefficients{1,6};
    %
    %     end
    % end
    %
    % if triggerValue == '2'
    %     save([outpath '/ageEstimate_20hz.mat'], 'ageEstimate');
    %     save([outpath '/ageTvalue_20hz.mat'], 'ageTvalue');
    %     save([outpath '/agePvalue_20hz.mat'], 'agePvalue');
    %     save([outpath '/intEstimate_20hz.mat'], 'intEstimate');
    %     save([outpath '/intTvalue_20hz.mat'], 'intTvalue');
    %     save([outpath '/intPvalue_20hz.mat'], 'intPvalue');
    %
    %     writematrix(ageEstimate, [outpath '/ageEstimate_20hz.csv']);
    %     writematrix(ageTvalue, [outpath '/ageTvalue_20hz.csv']);
    %     writematrix(agePvalue, [outpath '/agePvalue_20hz.csv']);
    %     writematrix(intEstimate, [outpath '/intEstimate_20hz.csv']);
    %     writematrix(intTvalue, [outpath '/intTvalue_20hz.csv']);
    %     writematrix(intPvalue, [outpath '/intPvalue_20hz.csv']);
    %
    % elseif triggerValue =='3'
    %     save([outpath '/ageEstimate_30hz.mat'], 'ageEstimate');
    %     save([outpath '/ageTvalue_30hz.mat'], 'ageTvalue');
    %     save([outpath '/agePvalue_30hz.mat'], 'agePvalue');
    %     save([outpath '/intEstimate_30hz.mat'], 'intEstimate');
    %     save([outpath '/intTvalue_30hz.mat'], 'intTvalue');
    %     save([outpath '/intPvalue_30hz.mat'], 'intPvalue');
    %
    %     writematrix(ageEstimate, [outpath '/ageEstimate_30hz.csv']);
    %     writematrix(ageTvalue, [outpath '/ageTvalue_30hz.csv']);
    %     writematrix(agePvalue, [outpath '/agePvalue_30hz.csv']);
    %     writematrix(intEstimate, [outpath '/intEstimate_30hz.csv']);
    %     writematrix(intTvalue, [outpath '/intTvalue_30hz.csv']);
    %     writematrix(intPvalue, [outpath '/intPvalue_30hz.csv']);
    %
    % elseif triggerValue =='4'
    %     save([outpath '/ageEstimate_40hz.mat'], 'ageEstimate');
    %     save([outpath '/ageTvalue_40hz.mat'], 'ageTvalue');
    %     save([outpath '/agePvalue_40hz.mat'], 'agePvalue');
    %     save([outpath '/intEstimate_40hz.mat'], 'intEstimate');
    %     save([outpath '/intTvalue_40hz.mat'], 'intTvalue');
    %     save([outpath '/intPvalue_40hz.mat'], 'intPvalue');
    %
    %     writematrix(ageEstimate, [outpath '/ageEstimate_40hz.csv']);
    %     writematrix(ageTvalue, [outpath '/ageTvalue_40hz.csv']);
    %     writematrix(agePvalue, [outpath '/agePvalue_40hz.csv']);
    %     writematrix(intEstimate, [outpath '/intEstimate_40hz.csv']);
    %     writematrix(intTvalue, [outpath '/intTvalue_40hz.csv']);
    %     writematrix(intPvalue, [outpath '/intPvalue_40hz.csv']);
    % end
    %

    % sigEstimates = ageEstimate.*(agePvalue<0.05);
    %
    % % Create a heatmap
    % imagesc(times, freqs, ageEstimate);
    % % Customize the plot as needed
    % colormap('jet'); % Adjust the colormap as needed
    % colorbar;
    % title('Time-Frequency by age linear estimates');
    % xlabel('Time');
    % ylabel('Frequency');
    %
    %
    % % Create a heatmap
    % imagesc(times, freqs, sigEstimates);
    % % Customize the plot as needed
    % colormap('jet'); % Adjust the colormap as needed
    % colorbar;
    % title('Time-Frequency by age sig estimates');
    % xlabel('Time');
    % ylabel('Frequency');
    %
    %
    % % Create a heatmap
    % imagesc(times, freqs, -log10(agePvalue));
    % % Customize the plot as needed
    % colormap('jet'); % Adjust the colormap as needed
    % colorbar;
    % title('Time-Frequency by age linear p values');
    % xlabel('Time');
    % ylabel('Frequency');
    %
    %
    % % Create a heatmap
    % imagesc(times, freqs, intEstimate);
    % % Customize the plot as needed
    % colormap('jet'); % Adjust the colormap as needed
    % colorbar;
    % title('Time-Frequency by intercept linear estimates');
    % xlabel('Time');
    % ylabel('Frequency');


    allSubs = [];
    % create a giant array for all subs and all TF values
    for c = 1:size(mergedCellArray,1)
        subTFarray = abs(mergedCellArray{c,2});
        % Split the char variable into two parts using '_'
        splitValues = strsplit(mergedCellArray{s,1}, '_');

        [timeGrid, freqGrid] = meshgrid(times, freqs);
        % Reshape the matrix and grids into a 3-column array
        subDataArray = [repmat(str2double(splitValues), 9600, 1),timeGrid(:), freqGrid(:), subTFarray(:), ];

        allSubs = [allSubs; subDataArray];
    end

    allSubstable = array2table(allSubs, 'VariableNames', {'lunaID','visitDate','age', 'time','freqs', 'ITC'});

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
