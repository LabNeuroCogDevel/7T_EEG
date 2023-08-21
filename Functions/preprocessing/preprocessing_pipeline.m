function [] = preprocessing_pipeline(inputfile, outpath, lowBP, topBP, FLAG, condition, task, varargin)

% what file are we using
if ~exist(inputfile,'file'), error('inputfile "%s" does not exist!', inputfile), end
[d, currentName, ext ] = fileparts(inputfile);

%% cap locations
cap_location = fullfile(eeglabpath,'/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
if ~exist(cap_location, 'file'), error('cannot find file for 128 channel cap: %s', cap_location), end
correction_cap_location = hera('Projects/7TBrainMech/scripts/eeg/Shane/resources/ChanLocMod128to64.ced');
if ~exist(correction_cap_location, 'file'), error('cannot find file for correction 128 channel cap: %s', correction_cap_location), end

%% Files
subj_files = file_locs(inputfile, outpath, task);

% to know how far your script is with running
fprintf('==========\n%s:\n\t Initial Preprocessing(%s,%f,%f,%s)\n',...
    currentName, inputfile, lowBP, topBP, outpath)


% where to save things
filter_folder = 'filtered';
chanrj_folder = 'channels_rejected';
epoch_folder = 'epoched';
icawholein_folder ='rerefwhole';
epoch_rj_marked_folder = 'marked_epochs';
epochrj_folder = 'rejected_epochs';
icaout = fullfile(outpath, 'ICA');
icawholeout = fullfile(outpath, 'ICAwhole');

% and what files will we create
rerefwhole_name = [currentName '_rerefwhole'];
chrm_name   = [currentName '_badchannelrj'];
epochrj_name = [currentName '_epochs_rj'];
% FIXME: these is not actually used!? but is recoreded in data_removed WF20190911
datarm_name = [currentName '_baddatarj'];
epochrm_name = [currentName '_badepochrj'];
icawholeout_name = [currentName '_rerefwhole_ICA'];
% TODO: collect other pop_saveset filenames here

commonPlus = {'AFz','C1','C2','C3','C4','C5','C6','CP1','CP2','CP3','CP4',...
    'CP5','CP6','CPz','Cz','F1','F2','F3','F4','F5','F6','F7','F8','FC1',...
    'FC2','FC3','FC4','FC5','FC6','FCz','Fp1','Fp2','FT10','FT9','Fz','I1',...
    'I2','O1','O2','Oz','P1','P10','P2','P3','P4','P5','P6','P7','P8','P9',...
    'PO10','PO3','PO4','PO7','PO8','PO9','POz','Pz','T7','T8',...
    'AF8','AF7','AF4','AF3'};

epochrj = fullfile(outpath, epochrj_folder, [epochrj_name '.set']);

if condition == 1
    icawholeoutFile = fullfile( icawholeout, [icawholeout_name '.set']);
    
else
    icawholeoutFile = 'no';
end

if exist(icawholeoutFile, 'file')
    warning('%s already complete (have "%s")! todo load from file', currentName, icawholeout_name)
    return
end

if condition == 1
    xEEG = load_if_exists(subj_files.filter);
end

if isstruct(xEEG)
    [ALLEEG EEG CURRENTSET] = pop_newset([], xEEG, 0,...
        'setname',currentName,...
        'gui','off');
else
    
    %% load EEG set and re- referance
    EEG = pop_loadset(inputfile);
    
    if size(EEG.data,1) < 100
        Flag128 = 0;
        EEG = pop_reref(EEG, [65 66]); %does this "restore" the 40dB of "lost SNR" ? -was it actually lost? ...this is potentially undone by PREP
        EEG = eeg_checkset(EEG);
        % find the cap to use
    else
        %[129 130] are the mastoid externals for the 128 electrode
        EEG = pop_reref(EEG, [129 130]); %does this "restore" the 40dB of "lost SNR" ? -was it actually lost? ...this is potentially undone by PREP
        EEG = eeg_checkset(EEG);
    end
    
    %stores EEG set in ALLEEG, give setname
    ALLEEG = [];
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,...
        'setname',currentName,...
        'gui','on');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    
    EEG.subject = currentName(1:findstr(currentName,'mgs')-2);
    EEG.condition =  currentName(findstr(currentName,'mgs'):end);
    
    %% Filtering
    
    %band-pass filter between low and max boundle (or 1 and 90 Hz)
    % TODO: params are
    %  EEG, locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz, minphase);
    % why 3380, 0, [], 0
    % > Warning: Transition band is wider than maximum stop-band width.
    % For better results a minimum filter order of 6760 is recommended.
    % Reported might deviate from effective -6dB cutoff frequency.
    % [FLAGexist] = checkdone(fullfile(outpath,filter_folder,[currentName '_filtered']));
    
    EEG = pop_eegfiltnew(EEG, lowBP, topBP, 3380, 0, [], 0);
    % filtorder = 3380 - filter order (filter length - 1). Mandatory even. performing 3381 point bandpass filtering.
    % pop_eegfiltnew() - transition band width: 0.5 Hz
    % pop_eegfiltnew() - passband edge(s): [0.5 90] Hz
    % pop_eegfiltnew() - cutoff frequency(ies) (-6 dB): [0.25 90.25] Hz
    % pop_eegfiltnew() - filtering the data (zero-phase)
    
    %give a new setname and overwrite unfiltered data
    EEG = pop_editset(EEG,'setname',[currentName '_bandpass_filtered']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    
    % %50 hz notch filter: 47.5-52.5
    % EEG = pop_eegfiltnew(EEG, 47.5,52.5,826,1,[],0);
    
    %% Resample data
    
    % Downsample the data to 150Hz using antialiasing filter
    EEG = pop_resample(EEG, 150, 0.8, 0.4); %0.8 is fc and 0.4 is dc. Default is .9 and .2. We dont know why Alethia changed them
    
    % Downsample the data to 512Hz
    % EEGb = pop_resample( EEG, 512);
    EEG = eeg_checkset(EEG);
    
    %change setname
    EEG = pop_editset(EEG,'setname',[currentName '_filtered']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    %save filtered data
    EEG = pop_saveset( EEG, 'filename',[currentName '_filtered'], ...
        'filepath',fullfile(outpath,filter_folder));
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end

%% CHANNELS
% remove external channels
EEG = pop_select( EEG,'nochannel',{'EX5' 'EX6' 'EX7' 'EX8' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp' 'FT7' 'FT8' 'TP7' 'TP8' 'TP9' 'TP10'});

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%import channel locations

%change this maybe
% eeglab should have already been added with addpath

EEG=pop_chanedit(EEG, 'lookup', cap_location);

if size(EEG.data,1) > 100
    EEG = pop_select( EEG,'channel',commonPlus);
    EEG = pop_chanedit(EEG, 'load', {correction_cap_location 'filetype' 'autodetect'});
    % 128    'AF8' --> 64    'AF6'
    % 128    'AF7' --> 64    'AF5'
    % 128    'AF4' --> 64    'AF2'
    % 128    'AF3' --> 64    'AF1'
    Flag128 = 1;
end
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% bad channel rejection

%bad channels are in general the channels that you can't save even if you
%reject 5-10% of your datapoints

%different options for channel rejection are displayed. Option 3 is ued.

% %1.look at standard deviation in bar plots and remove the channels
% with big std's --> manual
% stdData = std(EEG.data,0,2);
% figure(idx); bar(stdData)

% %2.kurtosis
% EEG = pop_rejchan(EEG, 'elec',[1:8] ,'threshold',5,'norm','on','measure','kurt');

%3.clean_rawdata
if condition == 1
    
    xEEG = load_if_exists(subj_files.chanrj);
end

originalEEG = EEG;
if isstruct(xEEG)
    [ALLEEG EEG] = eeg_store(ALLEEG, xEEG, CURRENTSET);
else
    EEG = clean_rawdata(EEG, 8, [0.25 0.75], 0.7, 5, 15, 0.3); % we can discuss that
    % vis_artifacts(EEG,originalEEG,'NewColor','red','OldColor','black');
    %change setname
    EEG = pop_editset(EEG,'setname', chrm_name);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG);
    
    EEG = pop_saveset( EEG,'filename', chrm_name, ...
        'filepath',sprintf('%s/%s/',outpath,chanrj_folder));
end

if condition == 1
    xEEG = load_if_exists(subj_files.rerefwhole_name);
end

if isstruct(xEEG)
    [ALLEEG EEG] = eeg_store(ALLEEG, xEEG, CURRENTSET);
else
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    if ~any(find(cellfun (@any,regexpi (fieldnames(EEG.etc), 'clean_channel_mask'))));
        EEG.etc.clean_channel_mask=42;
    else
    end
    
    %save the channels that were rejected in a variable
    channels_removed{1} = chrm_name; %setname
    channels_removed{2} = setdiff({originalEEG.chanlocs.labels},{EEG.chanlocs.labels}, 'stable');
    channels_removed{3} = find(EEG.etc.clean_channel_mask==0);
    %also save the channels that were rejected in the EEG struct
    EEG.channels_rj = channels_removed{2};
    EEG.channels_rj_nr = length(EEG.channels_rj);
    
    %save the proportion of the dataset that were rejected in a variable
    data_removed{1} = datarm_name; %setname
    data_removed{2} = length(find(EEG.etc.clean_sample_mask==0))/EEG.srate;%
    data_removed{3} = length(find(EEG.etc.clean_sample_mask==0))/length(EEG.etc.clean_sample_mask);%
    %also save the data that were rejected in the EEG struct
    EEG.data_rj    = data_removed{2};
    EEG.data_rj_nr = data_removed{3};
    
    %% interpolate channels
    % POSSIBLE PROBLEMS
    %  - injecting extra channels ontop of expected 64 (n>64)
    %  - 128 missing expected labels, adding too few back (n<64)
    if Flag128 == 1
        nchan = 64;
        ngood = length(EEG.chanlocs);
        %  128 cap doesn't have exactly the same postions as 64
        % remove 4 that are in the wrong place and reinterpret
        % AND interp any bad channels
        % do this by removing the 4 128weirdos
        % from the already trimmed (no bad channels) in EEG.chanlocs
        
        need_128interp = [2  3  35  36 ];
        % get the names of those to remove
        n128name = {originalEEG.chanlocs(need_128interp).labels};
        % should always be {'AF5','AF1','AF2','AF6'} ??
        
        % find where they are in current EEG files (if they haven't already been removed)
        n128here_idx = find(ismember({EEG.chanlocs.labels},n128name));
        % keep those that aren't the ones we matched
        % remove from chanlocs, data and update nbcan
        % WARNING -- who knows what else we should have changed to update the set info!
        keep_idx = setdiff(1:ngood, n128here_idx);
        EEG.chanlocs = EEG.chanlocs(keep_idx);
        EEG.data = EEG.data(keep_idx,:);
        EEG.nbchan = length(keep_idx);
        
        %EEG_i = pop_interp(EEG, interp_ch, 'spherical');
        fprintf('%d channels in orig; want to interpolate %d bad and move %d\n',...
            originalEEG.nbchan, nchan - ngood, length(need_128interp))
        EEG_i = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
        
        % could swap these channels (they're close, but not the same)
        % BUT WE DONT
        % 128    'AF7' --> 64    'AF5' In this point channel 2
        % 128    'AF3' --> 64    'AF1' In this point channel 3
        % 128    'AF4' --> 64    'AF2' In this point channel 35
        % 128    'AF8' --> 64    'AF6' In this point channel 36
        % lines above modify channel information and pocition in data to make
        %  it the same for 64 and 128 cap
        
        % need to do destructive swapping. need a copy
        EEG = EEG_i;
        EEG.chanlocs(2) = EEG_i.chanlocs(3);%EEG.chanlocs(2) must by 'AF1' in 64 cap
        EEG.chanlocs(3) = EEG_i.chanlocs(2);%EEG.chanlocs(3) must by 'AF5' in 64 cap
        %     EEG.chaninfo.filecontent(4,:) = EEG_i.chaninfo.filecontent(3,:); This
        %     is not necessary i think, but just in case...
        %     EEG.chaninfo.filecontent(4,1) = '3';
        %     EEG.chaninfo.filecontent(3,:) = EEG_i.chaninfo.filecontent(4,:);
        %     EEG.chaninfo.filecontent(3,1) = '2';
        EEG.data(2,:) = EEG_i.data(3,:);% ALERT ALERT Lines latelly added
        EEG.data(3,:) = EEG_i.data(2,:);% ATERT ALERT Lines latelly added
    else
        EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
    end
    % eeglab redraw
    %% re-reference: average reference
    if FLAG
        EEG = pop_reref( EEG, []);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,...
            'setname',[currentName '_avref'],...
            'gui','off');
        
    else
    end
    
    %save whole rereferenced data for ICA whole
    %save epochs rejected EEG data
    EEG = pop_editset(EEG, 'setname', rerefwhole_name);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    EEG = pop_saveset(EEG, 'filename', rerefwhole_name, 'filepath', fullfile(outpath,icawholein_folder));
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end

%% Whole data ICA run

icawholein = fullfile(outpath, icawholein_folder, [rerefwhole_name '.set']);

if ~exist(subj_files.icawhole, 'file')
    runICAs(icawholein, icawholeout)
else
    fprintf('have %s, not rerunning\n', subj_files.icawholeout)
end




























EEGfileNames = dir([path_data '/*.set']);


for fidx = 285:length(EEGfileNames)
    filename = EEGfileNames(fidx).name;
    locs = file_locs(fullfile(path_data,filename));
    if exist(locs.ICAwholeClean, 'file')
        fprintf('skipping; already created %s\n', locs.ICAwholeClean);
        continue
    end

    selectcompICA
end

% Clean epochs to remove
path_data = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/ICAwholeClean/');
epoch_folder =  hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/epoch/');
epoch_rj_marked_folder = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/epochclean/');

%% added to code to be pushed through manually 
% redo = [10173 10195 10644 10997 11299 11451 11543 11630 11631 11634 11640 11645 11661 11664 11669 11673];
% for i = 1:length(redo) 
%     file(i,:) = sprintf('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/ICAwholeClean/%d', redo(i));
%     EEGfileNames(i,:) = dir(hera([file(i,:) '*_icapru.set']));
% end

%%
EEGfileNames = dir([path_data, '/*_icapru.set']);

revisar = {};
for currentEEG = 1:size(EEGfileNames,1)
    filename = [EEGfileNames(currentEEG).name];
    inputfile = [path_data,filename];
revisar{currentEEG} = epochlean(inputfile,epoch_folder,epoch_rj_marked_folder);
end
% save revisar revisar

%% Epoching small 
clear all
close all
%% paths
addpath(genpath('Functions'));
% addpath(hera('Projects/7TBrainMech/scripts/eeg/toolbox/eeglab14_1_2b'))
%% files
folder = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/epochclean/');
%% added to code to run manually 
% redo = [10173 10195 10644 10997 11299 11451 11543 11630 11631 11634 11640 11645 11661 11664 11669 11673];
% for i = 1: length(redo) 
%     file(i,:) = sprintf('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/epochclean/%d', redo(i));
%     setfiles0(i,:) = dir(hera([file(i,:) '*ICA_icapru_epochs_rj.set']));
% end
%%

filter = [folder,'*ICA_icapru_epochs_rj.set'];
setfiles0 = dir(filter);
setfiles = {};
for ica = 1:length(setfiles0)
setfiles{ica,1} = fullfile(folder, setfiles0(ica).name); % cell array with EEG file names
% setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end


%% Homogenize Chanloc
datapath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/ICAwholeClean');
% datapath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/epochclean');
outpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/ICAwholeClean_homogenize');

%% added to code to run manually 
% redo = [10173 10195 10644 10997 11299 11451 11543 11630 11631 11634 11640 11645 11661 11664 11669 11673];
% 
% for i = 1: length(redo) 
%     file(i,:) = sprintf('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/epochclean/%d', redo(i));
%     setfiles0(i,:) = dir(hera([file(i,:) '*.set']));
% end
%%
setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};
for epo = 1:length(setfiles0)
setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
% setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

correction_cap_location = hera('Projects/7TBrainMech/scripts/eeg/Shane/resources/ELchanLoc.ced');
for i = 1:length(setfiles)
    homogenizeChanLoc(setfiles{i},correction_cap_location,outpath)
end

%% END HERE FOR SPECTRAL EVENTS 

%% Create Epochs for each condition
% Select epochs
outpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/PermutEpoch');
datapath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Prep/AfterWhole/ICAwholeClean_homogenize');
setfiles0 = dir([datapath,'/*.set']);
setfiles = {};
for epo = 1:length(setfiles0)
setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
% setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
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
%     misssubj{1,2} = missingbycond{1,1};
%     misssubj{1,3} = missingbycond{1,2};
%     misssubj{1,4} = missingbycond{1,3};
%     misssubj{1,5} = missingbycond{1,4};
%     misssubj{1,6} = missingbycond{1,5};
%     misssubj{1,6} = missingbycond{1,6};
%     
    
    misssubj{1,1} = missingbycond{1,1};
    save epocasXtipo misssubj
end

%% Permutaciones
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
% cfg.rois{5}=[cfg.rois{1} cfg.rois{2} cfg.rois{3} cfg.rois{4}]; % All
% cfg.rois{6} = [32 23 59 31 24 60 58];

%% Call function to make graphs 

CreateERPGraphs(condition, cfg, datapath, Bins, Resultados_folder);

