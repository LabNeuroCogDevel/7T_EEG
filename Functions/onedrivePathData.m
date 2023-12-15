function savedir = onedrivePathData(pathEnding)
allpaths = {'/home/sdm63/onedrive/data/',...
    '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/1. Recording & Stimulation of iEEG'};
for path = allpaths
    path = path{1};
    if exist(path, 'dir')
        savedir = fullfile(path,pathEnding);
        return;
    end
end

error('You dont have directory mounted')
end