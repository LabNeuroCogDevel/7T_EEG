% EEGLAB history file generated on the 17-Nov-2023
% ------------------------------------------------

EEG.etc.eeglabvers = '14.1.2'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = eeg_checkset( EEG );
EEG.etc.eeglabvers = '2022.1'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG.setname='10129_20180919_ss_Rem';
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG.setname='10129_20180919_ss_Rem_avref';
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG.setname='10129_20180919_ss_Rem_rerefwhole';
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '4'  }, [-0.2         0.8], 'newname', '10129_20180919_ss_Rem_rerefwhole_ICA pruned with ICA epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-199 0] ,[]);
EEG = eeg_checkset( EEG );
figure; pop_spectopo(EEG, 1, [-199.2188      798.8281], 'EEG' , 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
figure; pop_newtimef( EEG, 1, 17, [-199  799], [0] , 'topovec', 17, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'CP3', 'baseline',[NaN], 'plottype', 'curve', 'plotphase', 'off', 'padratio', 1, 'winsize', 102);
figure; pop_newtimef( EEG, 1, 17, [-199  799], [0] , 'topovec', 17, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'CP3', 'baseline',[NaN], 'freqs', [35 45], 'plotphase', 'off', 'nfreqs', 2);
figure; pop_newtimef( EEG, 1, 1, [-199  799], [3         0.8] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'Fp1', 'baseline',[0], 'plotphase', 'off', 'padratio', 1, 'winsize', 102);
figure; pop_newtimef( EEG, 1, 1, [-199  799], [3         0.8] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'Fp1', 'baseline',[NaN], 'plotphase', 'off', 'padratio', 1, 'winsize', 102);
figure; pop_newtimef( EEG, 1, 1, [-199  799], [0] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'Fp1', 'baseline',[0], 'plotphase', 'off', 'padratio', 1, 'winsize', 102);
figure; pop_newtimef( EEG, 1, 1, [-199  799], [2  15] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'Fp1', 'baseline',[0], 'plotphase', 'off', 'padratio', 1, 'winsize', 102);
figure; pop_newtimef( EEG, 1, 1, [-199  799], [2  15] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'Fp1', 'baseline',[0], 'plotphase', 'off', 'padratio', 1, 'winsize', 102);
figure; pop_newtimef( EEG, 1, 1, [-199  799], [2  15] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'Fp1', 'baseline',[0], 'freqs', [10 70], 'plotphase', 'off', 'padratio', 1);
figure; pop_newtimef( EEG, 1, 1, [-199  799], [2  15] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'Fp1', 'baseline',[0], 'plotphase', 'off', 'padratio', 4, 'winsize', 102);
