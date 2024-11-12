
for i= 1:length(didntRun)
    currentFile = didntRun{i};

    path_file = hera('Raw/EEG/7TBrainMech');
    outputpath = ['MGS' '/remarked/'];

    %directory of EEG data
    [path,folder] = fileparts(path_file);
    d =[path,'/',folder,'/'];

    namesOri = dir([d,currentFile(112:125),'*/*.bdf']);

    mgsIDX = find (cellfun (@any,regexpi ( {namesOri.name}.', 'mgs')));

    if numel(mgsIDX) > 1
        mgsIDX = mgsIDX(2);
        currentName = namesOri(mgsIDX).name(1:end-5);
    end

    if isempty(mgsIDX)
        continue;
    end

    currentName = namesOri(mgsIDX).name(1:end-4);
    d = [namesOri(mgsIDX).folder '/'];

    %% load raw EEG set
    rawEEG = pop_biosig([d currentName '.bdf']);

    eeglab redraw

   if class(rawEEG.event(1).type) == "char"
       for i = 1:length(rawEEG.event)
            rawEEG.event(i).type = str2double(rawEEG.event(i).type);
       end
   end

   for e = 1:length(rawEEG.event)
       if isnan(rawEEG.event(e).type)
           rawEEG.event(e).type = rawEEG.event(e).edftype;
       end
   end

    [micromed_time,mark]=make_photodiodevector(rawEEG);
    mark = mark - min(mark(mark>0))+61;
    mark(mark>65000) = 0;

    ending = mod(mark',10);

    simple = nan(size(mark));
    simple(mark == 254)= 1;
    simple(mark>=50 & mark<100)= 2;
    simple(mark>=100 & mark<150 & ((ending == 1) + (ending == 2))')= -3;
    simple(mark>=100 & mark<150 & ((ending == 3) + (ending == 4))')= 3;

    simple(mark>=150 & mark<200)= 4;
    simple(mark>=200 & mark<250 & ((ending == 1) + (ending == 2))')= -5;
    simple(mark>=200 & mark<250 & ((ending == 3) + (ending == 4))')= 5;

    % load preprocessed EEG
    EEG = pop_loadset(currentFile); % load in eeg file

    for e = 1:length(simple)
        EEG.event(e).type = simple(e);
        EEG.event(e).seconds = micromed_time(e);
    end

    EEG = pop_saveset(EEG, currentFile);
end


