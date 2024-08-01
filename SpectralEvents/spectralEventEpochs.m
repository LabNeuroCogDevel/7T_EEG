function [] = spectralEventEpochs()

%% Create Epochs for each condition
% Select epochs
outpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/PermutEpoch');
datapath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/ICAwholeClean_homogenize');
setfiles0 = dir([datapath,'/*.set']);
setfiles = {};
for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
end

mdirect{1}.mark = '1';mdirect{1}.outputdir = outpath;mdirect{1}.name = 'ITI';
mdirect{2}.mark =  '2';mdirect{2}.outputdir =  outpath;mdirect{2}.name = 'Fix';
mdirect{3}.mark =  '-3';mdirect{3}.outputdir = outpath;mdirect{3}.name = 'CueL';
mdirect{4}.mark =  '3';mdirect{4}.outputdir = outpath;mdirect{4}.name = 'CueR';
mdirect{5}.mark =  '-5';mdirect{5}.outputdir = outpath;mdirect{5}.name = 'MGSL';
mdirect{6}.mark =  '5';mdirect{6}.outputdir = outpath;mdirect{6}.name = 'MGSR';

%need to do delay separatly cause it is a longer epoch
mdirect{1}.mark =  '4'; mdirect{1}.outputdir = outpath; mdirect{1}.name = 'Delay';


wind = [2  4];
misssubj = cell(size(setfiles,1),size(mdirect,2)+1);
for currentEEG = 1:size(setfiles,1)
    missingbycond = extraepochs(setfiles{currentEEG},mdirect ,wind);

    [a,b,c] = fileparts(setfiles{currentEEG});
    misssubj{1,1} = b;

    misssubj{1,1} = missingbycond{1,1};
    save epocasXtipo misssubj
end

%% permutations
addpath(genpath('resources/permutaciones'));
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/PermutEpoch');
Resultados_folder = hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/ERPs');

% Firth part generates the grand average for each condition
condition{1} = 'ITI';
condition{2} = 'Fix';
condition{3} = 'CueR';
condition{4} = 'CueL';
condition{5} = 'Delay';
condition{6} = 'MGSR';
condition{7} = 'MGSL';

for i=1:length(condition)
    outputname = [condition{i},'_GAv'];
    Bins = grandaverage(datapath,outputname,condition{i});
end

%%
%Rois Parietal (R y L), Frontal lateral (R y L), partially invented
cfg.rois{1}=[5 6 7]; %left-frontal
cfg.rois{2}=[40 39 38]; %right-frontal
cfg.rois{3}=[22 23 31]; %left-parietal
cfg.rois{4}=[31 59 58]; %right-parietal


end