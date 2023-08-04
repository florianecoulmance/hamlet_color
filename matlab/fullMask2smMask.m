function smMask = fullMask2smMask(fullMask)

fullMask = logical(fullMask);
s = size(fullMask);
fullMask = reshape(fullMask,[s(1)*s(2) 1]);
ind = find(fullMask);
smMask.size = s;
smMask.pts = ind;

