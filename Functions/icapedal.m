[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','11771_20190416_rest_epochs_rj_ICA.set','filepath','C:\\Users\\Amelie\\Documents\\LNCDpasantia\\miniBatchs\\Results\\ICA\\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = eeg_checkset( EEG );
pop_topoplot(EEG,0, [1:30] ,'11771_20190416_rest_epochs_rj_ICA',[5 6] ,0,'electrodes','on');
EEG = eeg_checkset( EEG );
pop_selectcomps(EEG, [1:30] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_subcomp( EEG, [2  4], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','11771_20190416_rest_epochs_rj_ICAprun','gui','off'); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',2,'study',0); 