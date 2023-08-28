% EEGLAB history file generated on the 28-Aug-2023
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','11633_20180510_SS_Rem_rerefwhole_ICA_icapru.set','filepath','/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize/');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
EEG = pop_rmdat( EEG, {'4'},[-0.2 0.5] ,0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '4'  }, [-0.2         0.5], 'newname', '11633_20180510_SS_Rem_rerefwhole_ICA epochs', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-200.1953 0] ,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
EEG = eeg_checkset( EEG );
figure; pop_newtimef( EEG, 1, 6, [-200  499], [2  15] , 'topovec', 6, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'F3', 'baseline',[0], 'freqs', [10 70], 'plotphase', 'off', 'padratio', 4);
eeglab redraw;
