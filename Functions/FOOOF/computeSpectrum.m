
addpath(genpath('/resources/Euge/'))
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20220104'))
ft_defaults
run('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1/eeglab.m');

datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize');

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

channelValues = [1:64];

n = size(setfiles,1); %number of EEG sets to preprocess
 
allT = table('Size', [0, 4], 'VariableTypes', {'string', 'double', 'double','double'});
allT.Properties.VariableNames = {'Subject', 'Channel', 'freqs','power'};

for i = 1:n
    inputfile = setfiles{i};
    subject = inputfile(122:135);
      
    clear events
    clear data
    clear hdr
    clear avgData
    clear classLabels
    clear allData
    
    hdr = ft_read_header(inputfile);
    data = ft_read_data(inputfile, 'header', hdr);
    events = ft_read_event(inputfile, 'header', hdr);
    
    [psd, freqs] = pwelch(data', hamming(500), 128, [], 150); %(x,window,noverlap,f,fs)


    subT = table('Size', [0, 4], 'VariableTypes', {'string', 'double', 'double','double'});
    subT.Properties.VariableNames = {'Subject', 'Channel', 'freqs','power'};

    for c = 1:64

        chanT = table('Size', [0, 4], 'VariableTypes', {'string', 'double', 'double','double'});
        chanT.Properties.VariableNames = {'Subject', 'Channel', 'freqs','power'};

        channelPSD = psd(:,c);

        chanT.Subject(1:length(freqs),:)= subject;
        chanT.Channel(1:length(freqs),:)= c;
        chanT.freqs(1:length(freqs),:)= freqs;
        chanT.power(1:length(freqs),:)= channelPSD;

        subT = [subT ; chanT];


    end

    allT = [allT ; subT];
    
end

writetable(allT, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/PSDtable_allElectrodes.csv');


% merge7t = readtable('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv');
% ages = (merge7t(:, [1, 4,5])); % Select the first and third columns
% 
% % Define a custom function to concatenate the values with an underscore
% mergeColumns = @(col1, col2) strcat(num2str(col1), '_', num2str(col2));
% 
% % Apply the custom function to each row and create a new merged column
% ages.Subject = arrayfun(mergeColumns, ages.lunaid, ages.eeg_date, 'UniformOutput', false);
% 
% % Convert the cell array in the 'MergedColumn' to a cell array of strings
% ages.Subject = cellstr(ages.Subject);
% 
% mergedTable = innerjoin(ages, PSDtable);
% 
% 
% writetable(mergedTable, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/PSDtable_allElectrodes.csv');
% 
% % Find subjects within the specified age range
% youngest = table2cell(mergedTable(mergedTable.eeg_age >= 10 & mergedTable.eeg_age < 16,:));
% middle = table2cell(mergedTable(mergedTable.eeg_age >= 16 & mergedTable.eeg_age < 22,:));
% oldest = table2cell(mergedTable(mergedTable.eeg_age > 22,:));
% 
% avgYoungest = mean(cat(2, youngest{:, 5}), 2);
% avgMiddle = mean(cat(2, middle{:, 5}), 2);
% avgOldest = mean(cat(2,oldest{:,5}),2);
% 
% plot(freqs(15:175), log(avgYoungest(15:175)), 'DisplayName', '10-15');
% hold on;
% plot(freqs(15:175), log(avgMiddle(15:175)), 'DisplayName', '16-22');
% hold on;
% plot(freqs(15:175), log(avgOldest(15:175)), 'DisplayName', '22+'); 
% legend('Location', 'Best');  % You can adjust the 'Location' as needed
% 
% 
% % find the individual freq bands 
% alphaIdx = find(freqs >= 8 & freqs <= 12);
% betaIdx = find(freqs > 12 & freqs <= 30); 
% gammaIdx = find(freqs > 30 & freqs <= 50); 
% 
% for i = 1:length(youngest)
%     subjectPSD = youngest{i,5};
%     alphaPSDyoungest{i,:} = subjectPSD(alphaIdx,:);
%     betaPSDyoungest{i,:} = subjectPSD(betaIdx,:);
%     gammaPSDyoungest{i,:} = subjectPSD(gammaIdx,:);
% end
% 
% for i = 1:length(middle)
%     subjectPSD = middle{i,5};
%     alphaPSDmiddle{i,:} = subjectPSD(alphaIdx,:);
%     betaPSDmiddle{i,:} = subjectPSD(betaIdx,:);
%     gammaPSDmiddle{i,:} = subjectPSD(gammaIdx,:);
% 
% end
% 
% for i = 1:length(oldest)
%     subjectPSD = oldest{i,5};
%     alphaPSDoldest{i,:} = subjectPSD(alphaIdx,:);
%     betaPSDoldest{i,:} = subjectPSD(betaIdx,:);
%     gammaPSDoldest{i,:} = subjectPSD(gammaIdx,:);
% 
% end
% 
% avgalphaYoungest = mean(cat(2, alphaPSDyoungest{:,1}), 1); 
% avgbetaYoungest = mean(cat(2, betaPSDyoungest{:,1}), 1); 
% avggammaYoungest = mean(cat(2, gammaPSDyoungest{:,1}), 1); 
% 
% avgalphaMiddle = mean(cat(2, alphaPSDmiddle{:,1}), 1); 
% avgbetaMiddle = mean(cat(2, betaPSDmiddle{:,1}), 1); 
% avggammaMiddle = mean(cat(2, gammaPSDmiddle{:,1}), 1); 
% 
% avgalphaOldest = mean(cat(2, alphaPSDoldest{:,1}), 1); 
% avgbetaOldest = mean(cat(2, betaPSDoldest{:,1}), 1); 
% avggammaOldest = mean(cat(2, gammaPSDoldest{:,1}), 1); 
% 
% 
% 
% 
