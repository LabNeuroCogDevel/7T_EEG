

datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/sEEG_rawData/');
savePath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/sEEG_rawData/');

%load in all the files
setfiles0 = dir([datapath, 'P*/Rest/rest*ns2*']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1:length(setfiles0)
    idvalues(j,:) = (setfiles0(j).folder(72:78));
    
    if isfile([savePath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient
        load([savePath, idvalues(j,:), '/Rest/data.mat']);
        montage = readtable([savePath, idvalues(j,:), '/Rest/montage', idvalues(j,:), '.xlsx']);

        
        
        
    end
    
end

        
   
    

