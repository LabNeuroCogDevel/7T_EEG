
addpath(genpath('Functions'));

addpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/fooof/fooof_mat-main/fooof_mat'))
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20220104'))
ft_defaults

addpath('C:/Users/sdm42/AppData/Local/Programs/Python')

datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Resting_State/AfterWhole/ICAwholeClean_homogenize');

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


for i = 1:numSubj
    inputfile = setfiles{i};
    subject = inputfile(93:106);
    
    %     if any(strcmp(alreadyRun, subject))
    %         warning('%s already complete', subject)
    %         continue
    %
    %     else
    
    inputfile = setfiles{i};
    
    clear events
    clear data
    clear hdr
    clear avgData
    clear classLabels
    clear allData
    
    hdr = ft_read_header(inputfile);
    data = ft_read_data(inputfile, 'header', hdr);
    events = ft_read_event(inputfile, 'header', hdr);
    
    [psd, freqs] = pwelch(data', 500, [], [], 150);
    
    freqs = freqs';
    psd = psd';

% FOOOF settings
settings = struct();  % Use defaults
f_range = [30, 75];

% Run FOOOF
fooof_results = fooof(freqs, psd, f_range, settings);

% Print out the FOOOF Results
fooof_results
    
end
