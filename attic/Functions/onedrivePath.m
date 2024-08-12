function savedir = onedrivePath()
allpaths = {'/home/sdm63/onedrive', '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh'};
for path = allpaths
    path = path{1};
    if exist(path, 'dir')
        savedir = path;
        return;
    end
end

error('You dont have directory mounted')
end