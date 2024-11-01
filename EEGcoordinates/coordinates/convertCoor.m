
ourchannelECOG = ourChannels;
for i = 1:length(ourchannelECOG.Channel)
   ourchannelECOG.Channel(i).Type = 'ECOG';
end

% for i = 1:length(ALLEEG(1).chanlocs)
% eeglabcor(i,1) = ALLEEG(1).chanlocs(i).X; 
% eeglabcor(i,2) = ALLEEG(1).chanlocs(i).Y; 
% eeglabcor(i,3) = ALLEEG(1).chanlocs(i).Z; 
% 
% end

for i = 1:length(ourChannelsCortex.Channel)
    row = (ourChannelsCortex.Channel(i).Loc);
    locs2(i,1) = row(1,1);
    locs2(i,2) = row(2,1);
    locs2(i,3) = row(3,1);
end

for i = 1:length(ourChannelsCortex.Channel)
   names{i,1} = ourChannelsCortex.Channel(i).Name;
end


newCoor = cs_convert(mri, 'scs', 'mni', locs2)*1000;

newCoor = num2cell(newCoor);

for j = 1:length(newCoor)
    newCoor{j,4} = names{j,:};
end

writetable(cell2table(newCoor), '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/electrodeMNIcoordinatesCortex_20240202.csv')