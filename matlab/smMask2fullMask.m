function fullMask = smMask2fullMask(smMask)

s = smMask.size;
fullMask = zeros(s(1)*s(2),1);
pts = smMask.pts;
fullMask(pts) = 1;
fullMask = reshape(fullMask,[s(1),s(2)]);
fullMask = logical(fullMask);


