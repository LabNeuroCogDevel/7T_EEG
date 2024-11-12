function [EEG] =  remark(task, taskdirectory, dryrun)

path_file = hera('Raw/EEG/7TBrainMech');
outputpath = [taskdirectory '/remarked/'];

%directory of EEG data
[path,folder] = fileparts(path_file);
d =[path,'/',folder,'/'];

namesOri = dir([d,'*/*.bdf']);

if task == "MGS"
    mgsIDX = find (cellfun (@any,regexpi ( {namesOri.name}.', 'mgs')));

    for idx = mgsIDX'

        currentName = namesOri(idx).name(1:end-4);
        d = [namesOri(idx).folder '/'];

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
        EEG = pop_biosig([d currentName '.bdf']);
        EEG.setname=[currentName 'Rem']; %name the EEGLAB set (this is not the set file itself)

        eeglab redraw

        [micromed_time,mark]=make_photodiodevector(EEG);

        mark = mark - min(mark(mark>0));
        mark(mark>65000) = 0;
        % isi (150+x) and iti (254) are different
        %    event inc in 50: (50-200: cue=50,img=100,isi=150,mgs=200)
        %    category inc in 10 (10->30: None,Outdoor,Indoor)
        %    side inc in 1 (1->4: Left -> Right)
        %        61 == cue:None,Left
        %        234 == mgs:Indoor,Right
        %1 254 = ITI
        %2 50<cue<100 [50+(c 10,20,30)+(s,1-4)]
        %3 100<img.dot<150 [100+(c 10,20,30)+(s,1-4)] +/-
        %4 150<delay<200 [150+(c 10,20,30)]
        %5 200<mgs<250 [200+(c 10,20,30)+(s,1-4)]

        ending = mod(mark',10);

        simple = nan(size(mark));
        simple(mark == 254)= 1;
        simple(mark>=50 & mark<100)= 2;
        simple(mark>=100 & mark<150 & ((ending == 1) + (ending == 2))')= -3;
        simple(mark>=100 & mark<150 & ((ending == 3) + (ending == 4))')= 3;

        simple(mark>=150 & mark<200)= 4;
        simple(mark>=200 & mark<250 & ((ending == 1) + (ending == 2))')= -5;
        simple(mark>=200 & mark<250 & ((ending == 3) + (ending == 4))')= 5;

        for i=unique(simple)
            mmark=find(simple==i);
            if ~isempty(mmark)
                for j = 1:length(mmark)
                    EEG = pop_editeventvals(EEG,'changefield',{mmark(j) 'type' i});
                end
            end
        end

        EEG = pop_saveset( EEG, 'filename',[currentName '_Rem.set'],'filepath',outputpath);

    end


elseif task == "anti"
    antiIDX = find (cellfun (@any,regexpi ( {namesOri.name}.', 'anti')));

    for idx = antiIDX'

        currentName = namesOri(idx).name(1:end-4);
        d = [namesOri(idx).folder '/'];

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
        EEG = pop_biosig([d currentName '.bdf']);
        EEG.setname=[currentName 'Rem']; %name the EEGLAB setclc (this is not the set file itself)
    
        eeglab redraw

        [micromed_time, mark]=make_photodiodevector(EEG);
        if mark > 500

        iti = mode(mark); 

        mark = mark - iti + 254;
        end
        % 101-105: anti cue 
        % 151-155: target (dot on, look away)
        % 254 = back to fixation

        simple = nan(size(mark));
        simple(mark == 254)= 1; % (New ITI)
        simple(mark>=100 & mark<110)= 2; % (new Anti cue - red fixation cross, prepatory)
        simple(mark>=150 & mark<= 155)= 3; % (new dot on, look away) 


        for i=unique(simple)
            mmark=find(simple==i);
            if ~isempty(mmark)
                for j = 1:length(mmark)
                    EEG = pop_editeventvals(EEG,'changefield',{mmark(j) 'type' i});
                end
            end
        end

        EEG = pop_saveset( EEG, 'filename',[currentName '_Rem.set'],'filepath',outputpath);

    end




elseif task == "Resting_State"

    restIDX = find (cellfun (@any,regexpi ( {namesOri.name}.', 'rest')));

    for idx = restIDX'

        currentName = namesOri(idx).name(1:end-4);
        d = [namesOri(idx).folder '/'];

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
        EEG = pop_biosig([d currentName '.bdf']);
        EEG.setname=[currentName 'Rem']; %name the EEGLAB set (this is not the set file itself)

        eeglab redraw

        [micromed_time,mark]=make_photodiodevector(EEG);

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


        EEG = pop_resample(EEG, 512);

        EEG = pop_saveset( EEG, 'filename',[currentName '_Rem.set'],'filepath',outputpath);

    end

elseif task == "SNR"

    AudSSIDX = find (cellfun (@any,regexpi ( {namesOri.name}.', 'ss')));

    for idx = AudSSIDX'

        currentName = namesOri(idx).name(1:end-4);
        d = [namesOri(idx).folder '/'];

        % skip if we've already done
        finalfile=fullfile(outputpath, [currentName '_Rem.set']);
        % if exist(finalfile,'file')
        %     fprintf('already have %s\n', finalfile)
        %     continue
        % end

        fprintf('making %s\n',finalfile);

        %% load EEG set
        EEG = pop_biosig([d currentName '.bdf']);
        if isempty(EEG.event)
            continue
        else
            EEG.setname=[currentName '_Rem']; %name the EEGLAB set (this is not the set file itself)

            eeglab redraw

            [micromed_time,mark]=make_photodiodevector(EEG); % micromed_time: the time the trigger goes off; mark: the trigger value

            mark = mark - min(mark(mark>0));

            %changes the triggers to be single digit numbers
            for i=unique(mark)
                mmark=find(mark==i);
                if ~isempty(mmark)
                    for j = 1:length(mmark)
                        %             EEG.event(mmark).type = cond{i+1};
                        EEG = pop_editeventvals(EEG,'changefield',{mmark(j) 'type' i});
                    end
                end
            end

            EEG = pop_saveset(EEG, 'filename',[currentName '_Rem.set'],'filepath',outputpath);
        end
    end

end

end

