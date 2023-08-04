function flipFlag = isFishUpsideDown(mask)

flipFlag = 0;
% Assumes the fish is already horizontal
s = size(mask);
midpoint = round(s(2)/2);
fatpoint = find(sum(mask,1) == max(sum(mask,1)),1);

if fatpoint < midpoint
    flipFlag = 1;
end
