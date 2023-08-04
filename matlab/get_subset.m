function [Y,thesespecies] = get_subset(allImages,keepInd,numTextons,field)

keepImages = allImages(keepInd);

all_N = reshape([keepImages.textons],[numTextons numel(keepImages)])';

D = pdist2(all_N,all_N,@distChiSq);
Y = mdscale(D,3);
thesespecies = {keepImages.(field)};