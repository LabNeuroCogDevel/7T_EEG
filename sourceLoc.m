
% attempting source localization for HNN 

addpath(genpath(hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20220104'))
ft_defaults

datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize');
entropyPath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sourceLocalized/');

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


eeglab                        % start eeglab
eeglabPath = fileparts(which('eeglab'));                 % save its location
bemPath = fullfile(eeglabPath, 'plugins', 'dipfit', 'standard_BEM');    % load the dipfit plugin
EEG = pop_loadset(fullfile(eeglabPath, 'sample_data', 'eeglab_data_epochs_ica.set'));   % load the sample eeglab epoched dataset
EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile(bemPath, 'standard_vol.mat'), ...
           'coordformat','MNI','mrifile',fullfile(bemPath, 'standard_mri.mat'), ...
           'chanfile',fullfile(bemPath, 'elec', 'standard_1005.elc'), ...
           'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] , ...
           'chansel',[1:32] );

%% Leadfield Matrix calculation
dataPre = eeglab2fieldtrip(EEG, 'preprocessing', 'dipfit');   % convert the EEG data structure to fieldtrip

cfg = [];
cfg.channel = {'all'};
cfg.reref = 'yes';
cfg.refchannel = {'all'};
dataPre = ft_preprocessing(cfg, dataPre);

vol = load('-mat', EEG.dipfit.hdmfile);

cfg            = [];
cfg.elec       = dataPre.elec;
cfg.headmodel  = vol.vol;
cfg.resolution = 10;   % use a 3-D grid with a 1 cm resolution
cfg.unit       = 'mm';
cfg.channel    = { 'all' };
[sourcemodel] = ft_prepare_leadfield(cfg);

