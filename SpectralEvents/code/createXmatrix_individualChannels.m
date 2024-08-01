function [x, subject] = createXmatrix_individualChannels(setfiles, task, epoch)

if task == 'MGS' && epoch == 'Delay'
    for i = 1:numSubj
        inputfile = (setfiles{i});
        subject{i} = inputfile(95:108);

        if any(strcmp(alreadyRun, subject{i}))
            warning('%s already complete', subject{i})
            continue

        else

            clear events
            clear data
            clear hdr
            clear avgData
            clear classLabels
            clear allData
            clear ChannelTrialData

            hdr = ft_read_header(inputfile);
            data = ft_read_data(inputfile, 'header', hdr);
            events = ft_read_event(inputfile, 'header', hdr);

            cfg = [];
            cfg.dataset = inputfile;
            cfg.headerfile = inputfile;
            cfg.channel =   {'all', '-POz', '-EX3', '-EX4'};
            cfg.trialdef.eventtype = 'trigger';
            cfg.trialdef.eventvalue = '4';
            cfg.trialdef.prestim = 0;
            cfg.trialdef.poststim = 6;
            cfg.trialfun = 'ft_trialfun_general';
            cfg.trialdef.ntrials = 95;
            cfg.event = events;

            if ~ischar(cfg.event(1).value)

                newevents = cfg.event;
                for j = 1:length(newevents)
                    newevents(j).value = num2str(newevents(j).value);
                end
                cfg.event = newevents;

            end

            try

                [cfg]= ft_definetrial(cfg);

            catch
                warning('%s wont run through feildtrip', subject{i})
                wontRun{i} = subject{i};
                continue;
            end

            cfg.continuous = 'no';
            cfg.checkmaxfilter = 'true';
            [data] = ft_preprocessing(cfg);


            %redefine the trails to be between 3-4 seconds of the delay period
            cfg.trl = [];
            cfg = rmfield(cfg, 'trl');
            cfg.toilim = [3 4];
            data = ft_redefinetrial(cfg, data);

            for j = 1:length(data.trial) %trial

                trialData = (data.trial{1,j});
                ChannelTrialData{j} = trialData(:,:);

            end

            %     avgData = squeeze(mean(allData, 1));

            x{i} = ChannelTrialData;
            additionalSubject{i} = inputfile(95:108);

        end

    end

elseif task == 'MGS' && epoch == "Fix"
    for i = 1:numSubj
        inputfile = setfiles{i};
        subject = inputfile(95:108);

        if any(strcmp(alreadyRun, subject))
            warning('%s already complete', subject)
            continue

        else

            clear events
            clear data
            clear hdr
            clear avgData
            clear classLabels
            clear allData
            clear ChannelTrialData

            hdr = ft_read_header(inputfile);
            data = ft_read_data(inputfile, 'header', hdr);
            events = ft_read_event(inputfile, 'header', hdr);

            cfg = [];
            cfg.dataset = inputfile;
            cfg.headerfile = inputfile;
            cfg.channel =   {'all', '-POz', '-EX3', '-EX4'};
            cfg.trialdef.eventtype = 'trigger';
            cfg.trialdef.eventvalue = '2';
            cfg.trialdef.prestim = 0;
            cfg.trialdef.poststim = 1;
            cfg.trialfun = 'ft_trialfun_general';
            cfg.trialdef.ntrials = 95;
            cfg.event = events;

            if ~ischar(cfg.event(1).value)

                newevents = cfg.event;
                for j = 1:length(newevents)
                    newevents(j).value = num2str(newevents(j).value);
                end
                cfg.event = newevents;

            end

            try

                [cfg]= ft_definetrial(cfg);

            catch
                warning('%s wont run through feildtrip', subject)
                wontRun{i} = subject;
                continue;
            end


            [data] = ft_preprocessing(cfg);


            %redefine the trails to be between 3-4 seconds of the Fix period
            cfg.trl = [];
            cfg = rmfield(cfg, 'trl');
            cfg.toilim = [0 1];
            data = ft_redefinetrial(cfg, data);

            for j = 1:length(data.trial) %trial

                trialData = (data.trial{1,j});
                ChannelTrialData{j} = trialData(:,:);

            end

            %     avgData = squeeze(mean(allData, 1));

            x{i} = ChannelTrialData;
            additionalSubject{i} = inputfile(95:108);

        end
    end


elseif task == 'Resting_State'
    for i = 1:numSubj
        inputfile = setfiles{i};
        Subjects{:,i} = inputfile(93:106);


        if any(strcmp(alreadyRun, subject))
            warning('%s already complete', subject)
            continue

        else

            inputfile = setfiles{i};

            clear events
            clear data
            clear hdr
            clear avgData
            clear classLabels
            clear allData

            hdr = ft_read_header(inputfile);
            data = ft_read_data(inputfile, 'header', hdr);
            events = ft_read_event(inputfile, 'header', hdr);

            cfg = [];
            cfg.dataset = inputfile;
            cfg.headerfile = inputfile;
            cfg.channel =   {'all', '-POz'};
            cfg.trialdef.eventtype = 'trigger';
            cfg.trialdef.eventvalue = '16129';
            cfg.trialdef.prestim = 0;
            cfg.trialdef.poststim = 4;
            cfg.trialfun = 'ft_trialfun_general';
            cfg.trialdef.ntrials = 57;
            cfg.event = events;

            if ~ischar(cfg.event(1).value)

                newevents = cfg.event;
                for j = 1:length(newevents)
                    newevents(j).value = num2str(newevents(j).value);
                end
                cfg.event = newevents;

            end

            try

                [cfg]= ft_definetrial(cfg);

            catch
                warning('%s wont run through feildtrip', subject)
                wontRun{i} = subject;
                continue;
            end


            [data] = ft_preprocessing(cfg);


            %redefine the trails to be between 2 and 3 of the eyes open resting
            %state period
            cfg.trl = [];
            cfg = rmfield(cfg, 'trl');
            cfg.toilim = [2 3];
            data = ft_redefinetrial(cfg, data);

            for j = 1:length(data.trial) %trial

                trialData = (data.trial{1,j});
                ChannelTrialData{j} = trialData(:,:);

            end

            %     avgData = squeeze(mean(allData, 1));

            x{i} = ChannelTrialData;
        end

    end
end


