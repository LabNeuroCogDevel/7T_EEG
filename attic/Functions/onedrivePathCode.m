function savedir = onedrivePathCode(pathEnding)
allpaths = {'/home/sdm63/onedrive/code/',...
    '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/abelCode'};
for path = allpaths
    path = path{1};
    if exist(path, 'dir')
        savedir = fullfile(path,pathEnding);
        return;
    end
end

error('You dont have directory mounted')
end