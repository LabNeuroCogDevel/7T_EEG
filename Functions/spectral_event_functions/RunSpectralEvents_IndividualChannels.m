
addpath(('/Projects/7TBrainMech/scripts/eeg/toolbox/SpectralEvents-master/'))
addpath(genpath('Functions/spectral_event_functions/'));
addpath(hera('/Projects/7TBrainMech/scripts/fieldtrip-20220104'))
ft_defaults

% values to initalize
task = 'Resting_State'; % which task do you want to run. CHANGE ACCORDING TO MGS, REST, SNR, etc. *****
epoch = 'Delay'; % which epoch do you want to run, if REST there is no epoch
band = 'Gamma';
seconds = '3_4';


% In path
inpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data');

% Outpath
outpath = hera('Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/');

% task directory
taskdirectory = [inpath, '/', task];

datapath = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize'];

%load in all the delay files
setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    idvalues(j,:) = (setfiles0(j).name(1:14));
end

numSubj = length(idvalues);
x = cell(1,numSubj);
subject = cell(1,numSubj);
additionalSubject = cell(1,numSubj);
classLabels = cell(1,numSubj);

[num,txt,raw] = xlsread(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/Gamma_allChannels_TrialLevel_Delay_3_4.xls'));
alreadyRun = txt(:,1);

%for terminal:
% load(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/Gamma/alreadyRun.mat'))

%% Create X matrix for Spectral Event Toolbox

[x,subject] = createXmatix_individualChannels(setfiles, task, epoch);

%% Run spectral event toolbox
if band == 'Gamma'
    eventBand = [35,65]; %Frequency range of spectral events
    fVec = (1:70); %Vector of fequency values over which to calculate TFR

elseif band == 'Beta'
    eventBand = [13,30]; %Frequency range of spectral events
    fVec = (1:30); %Vector of fequency values over which to calculate TFR

elseif alpha == 'Alpha'
    eventBand = [8,12]; %Frequency range of spectral events
    fVec = (1:30); %Vector of fequency values over which to calculate TFR

elseif band == 'Theta'
    eventBand = [4,7]; %Frequency range of spectral events
    fVec = (1:30); %Vector of fequency values over which to calculate TFR

end

Fs = 150; % Sampling rate of time-series
findMethod = 1; % Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; % Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);


% remove emptry cells (from subjects that had already been run)
x = x(~cellfun('isempty',x));
SubjectsRan = additionalSubject(~cellfun('isempty',additionalSubject));

%% this is to run the spectral event analysis on every individual channel
for c = 1:63 %channels
    for s = 1:length(x)%subject
        clear oneChannel_allTrials
        for t = 1:length(x{1,s})  %trials

            oneChannel_allTrials(t,:) =  x{1,s}{1,t}(c,:); %subject s, trial t, channel c

            if task == 'MGS' && epoch == 'Delay'
                classLabels{s} = 4+zeros(1,length(x{1,s}) ); %needs to be the same length as the trial number

            elseif task == 'MGS' && epoch == "Fix"
                classLabels{s} = 2+zeros(1,length(x{1,s}) ); %needs to be the same length as the trial number

            elseif task == 'Resting_State'
                classLabels{s} = 16130+zeros(1,length(x{1,s}) ); %needs to be the same length as the trial number

            end
        end

        oneChannel_allSubjects{s} = oneChannel_allTrials';

    end
    [specEv_struct, TFRs, X, tvec] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis, oneChannel_allSubjects, classLabels);
    spevEV_struct_allChannels{c} = specEv_struct;
    TFRs_allChannels{c} = TFRs;
    X_allChannels{c} = X;

    fprintf('c = %.0f', c);
end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_allChannels_Delay.mat'], 'spevEV_struct_allChannels','-v7.3')
    save([savefolder band '_TFR_allChannels_Delay.mat'], 'TFRs_allChannels','-v7.3');
    save([savefolder band '_Xmatrix_allChannels_Delay.mat'], 'x','-v7.3')


elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_allChannels_Fix.mat'], 'spevEV_struct_allChannels','-v7.3')
    save([savefolder band '_TFR_allChannels_Fix.mat'], 'TFRs_allChannels','-v7.3');
    save([savefolder band '_Xmatrix_allChannels_Fix.mat'], 'x','-v7.3')

elseif task == 'Resting_State'
     savefolder = [outpath band '/'];
    save([savefolder band '_specEV_allChannels_RS.mat'], 'spevEV_struct_allChannels','-v7.3')
    save([savefolder band '_TFR_allChannels_RS.mat'], 'TFRs_allChannels','-v7.3');
    save([savefolder band '_Xmatrix_allChannels_RS.mat'], 'x','-v7.3')

end


% Average all subject TRFs together for each channel

for i = 1:63
    channel = TFRs_allChannels{i};
    for t = 1:size(channel{3},3)
        for s = 1:length(channel)
            if size(channel{s},3) ~= 197
                continue;
            else
                subject = channel{s};
                trial(:,:,s) = subject(:,:,t);
            end
        end
        avgTrial(:,:,t) = mean(trial, 3);
    end
    avgTFRperChannel{i} = avgTrial;
end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/'];
    save([savefolder band '_avgTFRperChannel_Delay.mat'], 'avgTFRperChannel','-v7.3')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/'];
    save([savefolder band '_avgTFRperChannel_Fix.mat'], 'avgTFRperChannel','-v7.3')

elseif task == 'Resting_State'
     savefolder = [outpath band '/'];
    save([savefolder band '_avgTFRperChannel_RS.mat'], 'avgTFRperChannel','-v7.3')

end

%% plot avg TFR for all subjects

for c = 1:63
    channel = avgTFRperChannel{c};
    avgTFR = squeeze(mean(channel(:,:,:),3));
    % Plot average raw TFR
    figure;
    subplot('Position',[0.08 0.75 0.75 0.17])
    imagesc([tVec(1) tVec(end)],[fVec(1) fVec(end)],avgTFR)
    x_tick = get(gca,'xtick');
    set(gca,'xtick',x_tick);
    set(gca,'ticklength',[0.0075 0.025])
    set(gca,'xticklabel',[])
    set(gca,'ytick',union(fVec([1,end]),eventBand))
    ylabel('Hz')
    pos = get(gca,'position');
    colormap jet
    cb = colorbar;
    cb.Position = [pos(1)+pos(3)+0.01 pos(2) 0.008 pos(4)];
    cb.Label.String = 'Spectral power';
    hold on
    line(tVec',repmat(eventBand,length(tVec),1)','Color','k','LineStyle',':')
    hold off
    title({'\fontsize{16}Raw TFR',['\fontsize{12}Channel ',num2str(c),', Trial class 4']})


    randTrial_inds = [1:10]; %randperm(numel(trial_inds),numSampTrials); %Sample trial indices
    eventBand_inds = find(fVec>=eventBand(1) & fVec<=eventBand(2)); %Indices of freq vector within eventBand

    % Plot 10 randomly sampled TFR trials
    for trl_i=1:10
        % Raw TFR trial
        rTrial_sub(trl_i) = subplot('Position',[0.08 0.75-(0.065*trl_i) 0.75 0.05]);
        %clims = [0 mean(eventThr(eventBand_inds))*1.5]; %Standardize upper limit of spectrogram scaling using the average event threshold
        %imagesc([tVec(1) tVec(end)],eventBand,TFR(eventBand_inds(1):eventBand_inds(end),:,trial_inds(randTrial_inds(trl_i))),clims)
        imagesc([tVec(1) tVec(end)],eventBand,channel(eventBand_inds(1):eventBand_inds(end),:,trl_i))
        x_tick_labels = get(gca,'xticklabels');
        x_tick = get(gca,'xtick');
        set(gca,'xtick',x_tick);
        set(gca,'ticklength',[0.0075 0.025])
        set(gca,'xticklabel',[])
        set(gca,'ytick',eventBand)
        rTrial_pos = get(gca,'position');
        colormap jet
        cb = colorbar;
        cb.Position = [rTrial_pos(1)+rTrial_pos(3)+0.01 rTrial_pos(2) 0.008 rTrial_pos(4)];


    end
end

%% find the adults

ages = readtable(hera('Projects\7TBrainMech\scripts\eeg\Shane\Results\Power_Analysis\Spectral_events_analysis\agefile_20210204.csv'));
ages.Subject = char(ages.Subject);
idvaluesTable = table((idvalues));
idvaluesTable.Properties.VariableNames{1} = 'Subject';
visitIdx = find((ages.visitno) == 1);
visit1 = ages(visitIdx, :);
test = rmmissing(outerjoin(idvaluesTable, visit1,'Keys','Subject', 'MergeKeys', 1));

adultIdx = find(categorical(test.Adult) == 'Adult');
childIdx = find(categorical(test.Adult) == 'Child');

% Average all adult TRFs together for each channel

for i = 1:63
    trial = nan(nFreq, nSamples, length(adultIdx));
    channel = TFRs_allChannels{i};
    for t = 1:size(channel,3)
        for s = 1:length(adultIdx)
            sub = adultIdx(s);
            if size(channel{sub},3) < t
                continue;
                trial(:,:,sub)
            else
                subject = channel{sub};
                trial(:,:,sub) = subject(:,:,t);
            end
        end
        avgTrial(:,:,t) = mean(trial, 3);
    end
    avgTFRperChannel{i} = avgTrial;
end

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/'];
    save([savefolder band '_avgTFRperChannel_Delay_childOnly.mat'], 'avgTFRperChannel','-v7.3')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/'];
    save([savefolder band '_avgTFRperChannel_Fix_childOnly.mat'], 'avgTFRperChannel','-v7.3')

elseif task == 'Resting_State'
     savefolder = [outpath band '/'];
    save([savefolder band '_avgTFRperChannel_RS_childOnly.mat'], 'avgTFRperChannel','-v7.3')

end

% average across time and frequency to get one measure of gamma power for each channel

for c = 1:63
    channel = avgTFRperChannel{c};
    avgchannel(1,c) = mean(mean(squeeze(mean(channel, 1)),2));
end

dbPower = pow2db(avgchannel);

channelLocations = readtable(hera('/Projects/7TBrainMech/scripts/eeg/Shane/Results/Power_Analysis/Spectral_events_analysis/ourChannels_locations.csv'));
channelLocations.Properties.VariableNames{1} = 'labels';

plot_topography('all', avgchannel, true, channelLocations)

for i = [36, 37, 38, 4, 5, 6]
    figure;
    imagesc([tVec(1) tVec(end)],[fVec(1) fVec(end)],mean(avgTFRperChannel{i},3))
    x_tick = get(gca,'xtick');
    set(gca,'xtick',x_tick);
    set(gca,'ticklength',[0.0075 0.025])
    set(gca,'xticklabel',[])
    set(gca,'ytick',union(fVec([1,end]),eventBand))
    ylabel('Hz')
    pos = get(gca,'position');

    colormap jet
    cb = colorbar;
    set(gca,'ticklength',[0.0075 0.025])
    set(gca,'xticklabel',x_tick_labels)
    xlabel('s')
end


n=round(.15*numel(avgchannel)); %how many channels is 15% of 63
idx = find(avgchannel > 1.1);

selectedChannels = channelLocations.labels(idx);


% plot the avg TFR and individual trials for the adult population for each channel
for c = 1:63
    channel = avgTFRperChannel{c};
    avgTFR = squeeze(mean(channel(:,:,:),3));
    % Plot average raw TFR
    figure;
    subplot('Position',[0.08 0.75 0.75 0.17])
    imagesc([tVec(1) tVec(end)],[fVec(1) fVec(end)],avgTFR)
    x_tick = get(gca,'xtick');
    set(gca,'xtick',x_tick);
    set(gca,'ticklength',[0.0075 0.025])
    set(gca,'xticklabel',[])
    set(gca,'ytick',union(fVec([1,end]),eventBand))
    ylabel('Hz')
    pos = get(gca,'position');
    colormap jet
    cb = colorbar;
    cb.Position = [pos(1)+pos(3)+0.01 pos(2) 0.008 pos(4)];
    cb.Label.String = 'Spectral power';
    hold on
    line(tVec',repmat(eventBand,length(tVec),1)','Color','k','LineStyle',':')
    hold off
    title({'\fontsize{16}Raw TFR',['\fontsize{12}Channel ',num2str(c),', Trial class 4']})


    randTrial_inds = randi([1 size(channel,3)], 1,10); %randperm(numel(trial_inds),numSampTrials); %Sample trial indices
    eventBand_inds = find(fVec>=eventBand(1) & fVec<=eventBand(2)); %Indices of freq vector within eventBand

    % Plot 10 randomly sampled TFR trials
    for trl_i=1:10
        % Raw TFR trial
        rTrial_sub(trl_i) = subplot('Position',[0.08 0.75-(0.065*trl_i) 0.75 0.05]);
        %clims = [0 mean(eventThr(eventBand_inds))*1.5]; %Standardize upper limit of spectrogram scaling using the average event threshold
        %imagesc([tVec(1) tVec(end)],eventBand,TFR(eventBand_inds(1):eventBand_inds(end),:,trial_inds(randTrial_inds(trl_i))),clims)
        imagesc([tVec(1) tVec(end)],eventBand,channel(eventBand_inds(1):eventBand_inds(end),:,trl_i))
        x_tick_labels = get(gca,'xticklabels');
        x_tick = get(gca,'xtick');
        set(gca,'xtick',x_tick);
        set(gca,'ticklength',[0.0075 0.025])
        set(gca,'xticklabel',[])
        set(gca,'ytick',eventBand)
        rTrial_pos = get(gca,'position');
        colormap jet
        cb = colorbar;
        cb.Position = [rTrial_pos(1)+rTrial_pos(3)+0.01 rTrial_pos(2) 0.008 rTrial_pos(4)];


    end

    subplot(rTrial_sub(end))
    %cb = colorbar;
    %cb.Position = [rTrial_pos(1)+rTrial_pos(3)+0.01 rTrial_pos(2) 0.008 rTrial_pos(4)];
    set(gca,'ticklength',[0.0075 0.025])
    set(gca,'xticklabel',x_tick_labels)
    xlabel('s')
end

%% group by region and then run spectral event toolbox
%frontal
frontalChannels_indx = find(contains(data.elec.label, ["F3","F4", "F1", "F2", "F5", "F6", "F7", "F8", "Fz"]));

for i = 1:length(x)
    %     clear frontal_trialLevel


    for j = 1:length(x{1, i})
        subjectData_trialLeavgTrialvel = x{1, i}{1,j};
        frontalChannels_trialLevel = subjectData_trialLevel(frontalChannels_indx,:);
        avg_frontalChannels_trialLevel = mean(frontalChannels_trialLevel);

        frontal_trialLevel(j,:) = avg_frontalChannels_trialLevel;

    end

    newX{i} = frontal_trialLevel';

end

for i = 1:length(newX)
    classLabels{i} = 4+zeros(1,size(newX{i},2));
end

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  newX ,classLabels);

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_selectiveFrontal_Delay' seconds '.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_selectiveFrontal_Delay' seconds '.mat'], 'TFRs','-v7.3')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_selectiveFrontal_Fix.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_selectiveFrontal_Fix.mat'], 'TFRs','-v7.3')

elseif task == 'Resting_State'
     savefolder = [outpath band '/'];
    save([savefolder band '_specEV_selectiveFrontal_RS.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_selectiveFrontal_RS.mat'], 'TFRs','-v7.3')
end

% occipital
occipitalChannels_indx = find(contains(data.elec.label, 'O'));
parietalChannels_indx = find(contains(data.elec.label, 'P'));
dlpfcChannels_indx = find(contains(data.elec.label, ["F3","F4"]));

clear newX

for i = 1:length(x)
    %     clear frontal_trialLevel


    for j = 1:length(x{1, i})
        subjectData_trialLevel = x{1, i}{1,j};
        occipitalChannels_trialLevel = subjectData_trialLevel(occipitalChannels_indx,:);
        avg_occipitalChannels_trialLevel = mean(occipitalChannels_trialLevel);

        occipital_trialLevel(j,:) = avg_occipitalChannels_trialLevel;

    end

    newX{i} = occipital_trialLevel';

end

for i = 1:length(newX)
    classLabels{i} = 4+zeros(1,size(newX{i},2));
end

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  newX ,classLabels);

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_Occipital_Delay' seconds '.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_Occipital_Delay' seconds '.mat'], 'TFRs','-v7.3')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_Occipital_Fix.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_Occipital_Fix.mat'], 'TFRs','-v7.3')

elseif task == 'Resting_State'
     savefolder = [outpath band '/'];
    save([savefolder band '_specEV_Occipital_RS.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_Occipital_RS.mat'], 'TFRs','-v7.3')
end


%Parietal
parietalChannels_indx = find(contains(data.elec.label, 'P'));
dlpfcChannels_indx = find(contains(data.elec.label, ["F3","F4"]));

clear newX

for i = 1:length(x)
    %     clear frontal_trialLevel


    for j = 1:length(x{1, i})
        subjectData_trialLevel = x{1, i}{1,j};
        parietalChannels_trialLevel = subjectData_trialLevel(parietalChannels_indx,:);
        avg_parietalChannels_trialLevel = mean(parietalChannels_trialLevel);

        parietal_trialLevel(j,:) = avg_parietalChannels_trialLevel;

    end

    newX{i} = parietal_trialLevel';

end

for i = 1:length(newX)
    classLabels{i} = 4+zeros(1,size(newX{i},2));
end

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  newX ,classLabels);


if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_Parietal_Delay' seconds '.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_Parietal_Delay' seconds '.mat'], 'TFRs','-v7.3')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_Parietal_Fix.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_Parietal_Fix.mat'], 'TFRs','-v7.3')

elseif task == 'Resting_State'
     savefolder = [outpath band '/'];
    save([savefolder band '_specEV_Parietal_RS.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_Parietal_RS.mat'], 'TFRs','-v7.3')
end

%DLPFC
dlpfcChannels_indx = find(contains(data.elec.label, ["F3","F4"]));

clear newX

for i = 1:length(x)
    %     clear frontal_trialLevel


    for j = 1:length(x{1, i})
        subjectData_trialLevel = x{1, i}{1,j};
        dlpfcChannels_trialLevel = subjectData_trialLevel(dlpfcChannels_indx,:);
        avg_dlpfcChannels_trialLevel = mean(dlpfcChannels_trialLevel);

        dlpfc_trialLevel(j,:) = avg_dlpfcChannels_trialLevel;

    end

    newX{i} = dlpfc_trialLevel';

end

for i = 1:length(newX)
    classLabels{i} = 4+zeros(1,size(newX{i},2));
end

[specEv_struct, TFRs, X] = spectralevents(cfg, eventBand, fVec, Fs, findMethod, vis,  newX ,classLabels);

if task == 'MGS' && epoch == 'Delay'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_DLPFC_Delay' seconds '.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_DLPFC_Delay' seconds '.mat'], 'TFRs','-v7.3')

elseif task == 'MGS' && epoch == 'Fix'
    savefolder = [outpath band '/'];
    save([savefolder band '_specEV_DLPFC_Fix.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_DLPFC_Fix.mat'], 'TFRs','-v7.3')

elseif task == 'Resting_State'
     savefolder = [outpath band '/'];
    save([savefolder band '_specEV_DLPFC_RS.mat'], 'specEv_struct','-v7.3')
    save([savefolder band '_TFRs_DLPFC_RS.mat'], 'TFRs','-v7.3')
end

