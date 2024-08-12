
% quick check that ICA_Path is valid
if ~ exist(ICA_Path,'dir'); error('no ICA_Path directory "%s"!', ICA_Path); end

% prefix ='_prep_allepoch_ica.set'
EEG = pop_loadset('filename',filename,'filepath',ICA_Path);
% EEG.setname = ['UG_' suj_n, '_',prefix(1:end-4)];
EEG = eeg_checkset( EEG );

% DB point on 28
[EEG, com] = pop_selectcomps(EEG, [1:30] );
pop_eegplot(EEG, 0, 1, 1);
pop_saveset( EEG, 'filename',[filename(1:end-4),'_pesos.set'],'filepath',CleanICApath);
% pop_eegplot( EEG, 0, 1, 1);
eeglab redraw
lowBP
% TEXTO = 'Lista la sel de comp? (If not put [0])';
% interp_user_def = input(TEXTO);
% 
% if interp_user_def~= 0
cmps = find(EEG.reject.gcompreject);
eeglab redraw
OUTEEG = pop_subcomp( EEG, [cmps], 1);
eeglab redraw
pop_saveset( OUTEEG, 'filename',[filename(1:end-4),'_icapru.set'],'filepath',CleanICApath);
% end
close all
