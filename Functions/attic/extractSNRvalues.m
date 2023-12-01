
%% Evoked activity

%% put all files into one ERSP
first50 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectERPs_first50.mat');
next51to100 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectERPs_51_100.mat');
next101to150 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectERPs_101_150.mat');
next151to200 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectERPs_151_200.mat');
next201to250 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectERPs_201_250.mat');
next251end = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectERPs_251_end.mat');

first50 = first50.subjectERSP;
next51to100 = next51to100.subjectERSP;
next101to150 = next101to150.subjectERSP;
next151to200 = next151to200.subjectERSP;
next201to250 = next201to250.subjectERSP;
next251end = next251end.subjectERSP;

% Create a logical index to identify empty cells
isEmptynext51to100 = cellfun('isempty', next51to100);
isEmptynext101to150 = cellfun('isempty', next101to150);
isEmptynext151to200 = cellfun('isempty', next151to200);
isEmptynext201to250 = cellfun('isempty', next201to250);
isEmptynext251end = cellfun('isempty', next251end);

% Use logical indexing to extract non-empty cells
notEmptyNext51to100 = next51to100(~isEmptynext51to100);
notEmptyNext101to150 = next101to150(~isEmptynext101to150);
notEmptyNext151to200 = next151to200(~isEmptynext151to200);
notEmptyNext201to250 = next201to250(~isEmptynext201to250);
notEmptyNext251end = next251end(~isEmptynext251end);


%concatenate the arrays together
allSubjectsERSP = [first50, notEmptyNext51to100, notEmptyNext101to150, notEmptyNext151to200, notEmptyNext201to250, notEmptyNext251end];
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsERP.mat','allSubjectsERSP')


%% put all files into one ITC
first50 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectITC_first50.mat');
next51to100 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectITC_51_100.mat');
next101to150 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectITC_101_150.mat');
next151to200 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectITC_151_200.mat');
next201to250 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectITC_201_250.mat');
next251end = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectITC_251_end.mat');

first50 = first50.subjectITC;
next51to100 = next51to100.subjectITC;
next101to150 = next101to150.subjectITC;
next151to200 = next151to200.subjectITC;
next201to250 = next201to250.subjectITC;
next251end = next251end.subjectITC;

% Create a logical index to identify empty cells
isEmptynext51to100 = cellfun('isempty', next51to100);
isEmptynext101to150 = cellfun('isempty', next101to150);
isEmptynext151to200 = cellfun('isempty', next151to200);
isEmptynext201to250 = cellfun('isempty', next201to250);
isEmptynext251end = cellfun('isempty', next251end);

% Use logical indexing to extract non-empty cells
notEmptyNext51to100 = next51to100(~isEmptynext51to100);
notEmptyNext101to150 = next101to150(~isEmptynext101to150);
notEmptyNext151to200 = next151to200(~isEmptynext151to200);
notEmptyNext201to250 = next201to250(~isEmptynext201to250);
notEmptyNext251end = next251end(~isEmptynext251end);


%concatenate the arrays together
allSubjectsITC = [first50, notEmptyNext51to100, notEmptyNext101to150, notEmptyNext151to200, notEmptyNext201to250, notEmptyNext251end];
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsITC.mat','allSubjectsITC')


%% put all files into one power base
first50 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectPowbase_first50.mat');
next51to100 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectPowbase_51_100.mat');
next101to150 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectPowbase_101_150.mat');
next151to200 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectPowbase_151_200.mat');
next201to250 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectPowbase_201_250.mat');
next251end = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/subjectPowbase_251_end.mat');

first50 = first50.subjectPowbase;
next51to100 = next51to100.subjectPowbase;
next101to150 = next101to150.subjectPowbase;
next151to200 = next151to200.subjectPowbase;
next201to250 = next201to250.subjectPowbase;
next251end = next251end.subjectPowbase;

% Create a logical index to identify empty cells
isEmptynext51to100 = cellfun('isempty', next51to100);
isEmptynext101to150 = cellfun('isempty', next101to150);
isEmptynext151to200 = cellfun('isempty', next151to200);
isEmptynext201to250 = cellfun('isempty', next201to250);
isEmptynext251end = cellfun('isempty', next251end);

% Use logical indexing to extract non-empty cells
notEmptyNext51to100 = next51to100(~isEmptynext51to100);
notEmptyNext101to150 = next101to150(~isEmptynext101to150);
notEmptyNext151to200 = next151to200(~isEmptynext151to200);
notEmptyNext201to250 = next201to250(~isEmptynext201to250);
notEmptyNext251end = next251end(~isEmptynext251end);


%concatenate the arrays together
allSubjectsSubjectPowbase = [first50, notEmptyNext51to100, notEmptyNext101to150, notEmptyNext151to200, notEmptyNext201to250, notEmptyNext251end];
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsPowbase.mat','allSubjectsSubjectPowbase')


%% average across time and frequency 
% baseline power
for i = 1:length(allSubjectsSubjectPowbase)
    for j = 1:6
        channel = allSubjectsSubjectPowbase{1,i}{1,j};
        avgChannel{j,1} = mean(channel,2);
    end
   subjectAvgPowbase{i} = avgChannel;

end
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsAvgPowbase.mat','subjectAvgPowbase')


% ERSP
for i = 1:length(allSubjectsERSP)
    for j = 1:6
        channel = allSubjectsERSP{1,i}{1,j};
        avgChannel{j,1} = mean(mean(channel));
    end
   subjectAvgERSP{i} = avgChannel;

end
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsAvgERSP.mat','subjectAvgERSP')


% ITC
for i = 1:length(allSubjectsITC)
    for j = 1:6
        channel = allSubjectsITC{1,i}{1,j};
        avgChannel{j,1} = abs(mean(mean(channel)));
    end
   subjectAvgITC{i} = avgChannel;

end
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsAvgITC.mat','subjectAvgITC')



%% create a table with all the info

largeT = table('Size', [0, 5], 'VariableTypes', {'string', 'double', 'double','double','double'});
largeT.Properties.VariableNames = {'Subject', 'Channel', 'ERSP', 'ITC','BaselinePower'};

for i = 1:length(idvalues)
    T = table('Size', [0, 5], 'VariableTypes', {'string', 'double', 'double','double','double'});
    T.Properties.VariableNames = {'Subject', 'Channel', 'ERSP', 'ITC','BaselinePower'};
    for c = 1:6
        T.ERSP(c) = subjectAvgERSP{:,i}{c,1};
        T.ITC(c) = subjectAvgITC{:,i}{c,1};
        if isempty(subjectAvgPowbase{:,i}{c,1})
            T.BaselinePower(c) = NaN;
        else
            T.BaselinePower(c) = subjectAvgPowbase{:,i}{c,1};
        end
        T.Channel(c) = channelValues(c);
    end
    T.Subject(1:6) = idvalues{i};

    largeT = [largeT; T];
end

writetable(largeT, '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/evokedActivityAllSubjects.csv');





%% Induced activity 

% put all files into one ERSP
first50 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/trialTFdata_1_50.mat');
next51to100 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/trialTFdata_51_100.mat');
next101to150 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/trialTFdata_101_150.mat');
next151to200 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/trialTFdata_151_200.mat');
next201to250 = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/trialTFdata_201_250.mat');
next251end = load('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/trialTFdata_251_end.mat');

first50 = first50.subjecttfData;
next51to100 = next51to100.subjecttfData;
next101to150 = next101to150.subjecttfData;
next151to200 = next151to200.subjecttfData;
next201to250 = next201to250.subjecttfData;
next251end = next251end.subjecttfData;

% Create a logical index to identify empty cells
isEmptynext51to100 = cellfun('isempty', next51to100);
isEmptynext101to150 = cellfun('isempty', next101to150);
isEmptynext151to200 = cellfun('isempty', next151to200);
isEmptynext201to250 = cellfun('isempty', next201to250);
isEmptynext251end = cellfun('isempty', next251end);

% Use logical indexing to extract non-empty cells
notEmptyNext51to100 = next51to100(~isEmptynext51to100);
notEmptyNext101to150 = next101to150(~isEmptynext101to150);
notEmptyNext151to200 = next151to200(~isEmptynext151to200);
notEmptyNext201to250 = next201to250(~isEmptynext201to250);
notEmptyNext251end = next251end(~isEmptynext251end);


%concatenate the arrays together
allSubjectsTF = [first50, notEmptyNext51to100, notEmptyNext101to150, notEmptyNext151to200, notEmptyNext201to250, notEmptyNext251end];
save('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR/allSubjectsTFdata.mat','allSubjectsTF', '-v7.3')


% need to take abs of the TF data to get power
for i = 1:length(allSubjectsTF)
    for j = 1:6
        channel = allSubjectsTF{1,i}{1,j};
        avgChannel{j,1} = abs(mean(mean(channel)));
    end
   subjectAvgITC{i} = avgChannel;

end

