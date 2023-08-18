function [EEG] = remarked_data_restingState(dryrun)
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Functions'));
run('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1/eeglab.m');

% should we actually run, or just say what we'd do (dry run)
% defaut to just printing, not actually running
if nargin < 1
    dryrun=0;
end
eeglab

path_file = hera('Raw/EEG/7TBrainMech');
outputpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Resting_State/remarked');

%directory of EEG data
[path,folder] = fileparts(path_file);
d =[path,'/',folder,'/'];
% d ='/Users/macbookpro/Documents/MBCS/BA/DMT_experiment/drive-download-20190719T191352Z-001/';
% names = dir([d,'*mgs.bdf']);

namesOri = dir([d,'*/*.bdf']);
restIDX = find (cellfun (@any,regexpi ( {namesOri.name}.', 'rest')));
% names = namesOri();
% 
% % names = dir([d,'*eyecal.bdf']);
% 
% names = {names(~[names.isdir]).name}; %cell array with EEG file names
% nr_eegsets = size(names,2); %number of EEG sets to preprocess

for idx = restIDX'
    
    currentName = namesOri(idx).name(1:end-4);
    d = [namesOri(idx).folder '/'];

    %to know how far your script is with running
    %disp(currentName);
 
    %% skip if we've already done
    finalfile=fullfile(outputpath, [currentName '_Rem.set']);
    if exist(finalfile,'file')
        fprintf('already have %s\n', finalfile)
        continue
    end
    if dryrun
        fprintf('want to run %s; set dryrun=0 to actually run\n', finalfile)
        continue
    end
    fprintf('making %s\n',finalfile);
  
    %% load EEG set
%     EEG = pop_biosig([d currentName '.bdf'],'ref',[65 66] );
    EEG = pop_biosig([d currentName '.bdf']);
    EEG.setname=[currentName 'Rem']; %name the EEGLAB set (this is not the set file itself)

    eeglab redraw
    
      [micromed_time,mark]=make_photodiodevector(EEG);
%     plot([micromed_time;micromed_time],[-100;0],'r')
%     figure(5);plot(micromed_time ,mark ,'r*')
    
    %%
    mark = mark - min(mark);
    mark(mark>65000) = 0;

    
    ending = mod(mark',10);
    

    for i=unique(mark)
        mmark=find(mark==i);
        if ~isempty(mmark)
            for j = 1:length(mmark)
                %             EEG.event(mmark).type = cond{i+1};
                EEG = pop_editeventvals(EEG,'changefield',{mmark(j) 'type' i});
            end
        end
    end
    %% GUARDO EL DATASET
    
EEG = pop_resample( EEG, 512);

EEG = pop_saveset( EEG, 'filename',[currentName '_Rem.set'],'filepath',outputpath);
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

end

clear all 
close all


end

