% EEGLAB history file generated on the 19-Oct-2022
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','11882_20220826_run_893726.18.tsv_Rem_take2.set','filepath','H:\\Projects\\7TBrainMech\\scripts\\eeg\\Shane\\Habit\\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, [65 66] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',75,'plotfreqz',1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'lookup','H:\\Projects\\7TBrainMech\\scripts\\eeg\\Shane\\Functions\\resources\\eeglab2022.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '13'  '14'  '15'  }, [-0.5         0.5], 'newname', 'BDF file epochs', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-500 0] ,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'gui','off'); 
EEG = eeg_checkset( EEG );
figure; pop_erpimage(EEG,1, [27],[[]],'O1',10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [27] EEG.chanlocs EEG.chaninfo } );
EEG = eeg_checkset( EEG );
figure; pop_plottopo(EEG, [1:65] , 'BDF file epochs', 0, 'ydir',1);
EEG = eeg_checkset( EEG );
figure; pop_timtopo(EEG, [-500      498.0469], [NaN], 'ERP data and scalp maps of BDF file epochs');
eeglab redraw;
