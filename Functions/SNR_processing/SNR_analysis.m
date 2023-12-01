
set(0,'DefaultFigureVisible','on'); %set figure visibility to off
%% set needed paths
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Functions/SNR_processing')));
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions/preprocessing'));

%% set initial values 
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize');
triggerValue = '4'; 
% channelValues = [4,5,6,36,37,38]; if you only want to run DLPFC
outpath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/SNR');

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

numSubj = length(idvalues);
errorSubjects = cell(1,numSubj);

%% run evoked and induced activity function
for i = 1:numSubj
    
    subject = idvalues{i};
    inputfile = setfiles{i};
    savePath = [outpath '/' subject '_SNRdata.csv'];
    
    if ~exist(savePath, 'file')
        fprintf('The file %s does not exist, running SNR code\n', savePath)
        
        [errorSubjects, subjectAvgERSP, subjectAvgITC, subjectAvgPowBase, subjectAvgTF] = ERSP_ITC(i, inputfile, triggerValue, outpath, subject,errorSubjects);
        
        [errorSubjectsEvoked, channelPower] = evokedActivity(i, inputfile, triggerValue, outpath, subject,errorSubjects);
        
        
        %% create a table with all the info
        if ~isempty(subjectAvgERSP) || ~isempty(channelPower)
            T = table('Size', [0, 7], 'VariableTypes', {'string', 'double', 'double','double','double','double','double'});
            T.Properties.VariableNames = {'Subject', 'Channel', 'ERSP', 'ITC','BaselinePower','Induced', 'Evoked'};
            
            for c = 1:64
                T.ERSP(c) = subjectAvgERSP{c};
                T.ITC(c) = subjectAvgITC{c};
                T.BaselinePower(c) = subjectAvgPowBase{c};
                T.Induced(c) = subjectAvgTF{c};
                T.Evoked(c) = channelPower(c,1);
                
                T.Channel(c) = c;
            end
            T.Subject(1:64) = subject;
            
            savePath = [outpath '/' subject '_SNRdata.csv'];
            writetable(T, savePath);
            
            clear subjectAvgERSP
            clear subjectAvgITC
            clear subjectAvgPowBase
            clear subjectAvgTF
            clear channelPower
        end
        disp(i);
        
    end
end



%% addressing the failed subjects
% 
% % Find the indices of empty cells
% emptyCellIndices = cellfun('isempty', errorSubjects);
% 
% % Use logical indexing to remove empty cells
% errorSubjects = errorSubjects(~emptyCellIndices);
% 
% for i = 1:length(errorSubjects)
%     sub = errorSubjects{i};
%     rawPath_file = hera('Raw/EEG/7TBrainMech'); % path to raw data 
% 
%     ssmatches = dir(fullfile(rawPath_file, ['/' sub '/*_SS*bdf'])); 
%     rawFile = fullfile(ssmatches(1).folder, ssmatches(1).name);
% 
%     EEG = pop_biosig(rawFile); %load in the raw file
%     [micromed_time,mark]=make_photodiodevector(EEG); %find the trigger values
% 
%     newMark = fixMarks(mark);
% 
%     processedPath = hera('Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize/');
% 
%     ssmatches = dir(fullfile(processedPath, [sub '*_SS*u.set'])); 
%     preprocessedFile = fullfile(ssmatches(1).folder, ssmatches(1).name);
% 
%     preprocessedEEG = pop_loadset(preprocessedFile); %load in the raw file
% 
%     for j = 1:length(EEG.event)
%         preprocessedEEG.event(j).urevent = EEG.event(j).urevent;
%         preprocessedEEG.event(j).latency = EEG.event(j).latency;
%         preprocessedEEG.event(j).type = newMark(j);
%     end
% 
%    preprocessedEEG = pop_saveset(preprocessedEEG, preprocessedFile);
% 
% end
% 
% 
