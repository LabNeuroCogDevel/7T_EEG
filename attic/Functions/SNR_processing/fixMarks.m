
function [newMark] = fixMarks(mark)

u = unique(mark);
keep = [];

for i = 1:length(u)
    ind = find(mark == u(i));

   if length(ind) >120
        keep = [keep,u(i)];

   end
end

keepIdx = ismember(mark,keep);

newMark = mark;
newMark(~keepIdx)=0;

if length(unique(newMark)) ~= 4
    fprintf('%d number of unique triggers, %d number of kept triggers, %d length of mark array \n',length(u),length(keepIdx), length(mark));
    error('bad triggers');
end

newMark = newMark - min(newMark(newMark>0))+2; % change them to single digit

