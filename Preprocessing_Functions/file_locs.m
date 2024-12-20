function loc = file_locs(filepath, savedir, task)
% From any file in the eeg pipeline, get all paths used
if task == "MGS"
    id_task_regexp='(?<id>\d{5}_\d{8})_(?<task>MGS)(?<extra>.*)';


elseif task == "SNR"
    id_task_regexp='(?<id>\d{5}_\d{8})_(?<task>ss)(?<extra>.*)';


elseif task == "Resting_State"
    id_task_regexp='(?<id>\d{5}_\d{8})_(?<task>REST)(?<extra>.*)';

elseif task == "anti"
    id_task_regexp='(?<id>\d{5}_\d{8})_(?<task>anti)(?<extra>.*)';

end

[filedir, filename, filext ] = fileparts(filepath);


% todo maybe this isn't a problem.
if ~exist(filepath, 'file'), error('file "%s" does not exist', filepath), end
if ~strmatch(filext, 'set'), error('file "%s" is not a .set file', filepath), end

% identify name
parts = regexp(filename, id_task_regexp, 'ignorecase','names');
if isempty(parts),error('%s failed to match id_task pattern (%s)', filename, id_task_regexp), end

subj_task = [parts(1).id '_' parts(1).task];

%% define file locations
mkname=@(f, ext) fullfile(savedir, f, [subj_task ext '.set']);
loc.remarked        = mkname('remarked',          '_Rem');
loc.filter          = mkname('filtered',          '_Rem_filtered');
loc.filterbp        = mkname('bandpass_filtered', '_bandpass_filtered') ;
loc.chanrj          = mkname('channels_rejected', '_Rem_badchannelrj');
loc.rerefwhole_name = mkname('rerefwhole',        '_Rem_rerefwhole');
loc.epoch           = mkname('epoched',           '_Rem_epochs');
loc.epoch_rj_marked = mkname('marked_epochs',     '_Rem_epochs_marked');
loc.epochrj         = mkname('rejected_epochs',   '_Rem_epochs_rj');
loc.icaout          = mkname('ICA',               '_Rem_epochs_rj_ICA');
loc.icawhole        = mkname('ICAwhole',          '_Rem_rerefwhole_ICA');
loc.SASICA          = mkname('ICA',               '_Rem_epochs_rj_ICA_SAS');
%TODO: add rerefwhole_name

%% stats and logic on files in/completed
loc.allfiles = fieldnames(loc);
files_exist = cellfun(@(x) exist(loc.(x),'file'), loc.allfiles);
loc.missing = loc.allfiles(~files_exist);
loc.ncomplete = sum(files_exist);
loc.is_finished = exist(loc.icaout, 'file') ~= 0;

% not single subject 
loc.epochClean      = mkname('AfterWhole/epochclean', '_Rem_rerefwhole_ICA_icapru_epochs_rj');
loc.ICAwholeClean      = mkname('AfterWhole/ICAwholeClean', '_Rem_rerefwhole_ICA_icapru');
loc.epochCleanHomongenize      = mkname('AfterWhole/epochclean_homogenize', '_Rem_rerefwhole_ICA_icapru_epochs_rj');
loc.ICACleanHomongenize      = mkname('AfterWhole/ICAwholeClean_homogenize', '_Rem_rerefwhole_ICA_icapru');


end
