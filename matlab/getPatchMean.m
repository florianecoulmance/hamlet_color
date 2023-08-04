function RGB = getPatchMean(I,mask)
RedLayer = double(I(:,:,1));
GreenLayer = double(I(:,:,2));
BlueLayer = double(I(:,:,3));

assert(size(RedLayer,1)==size(mask,1),'the mask is not the same size as the image')
assert(size(RedLayer,2)==size(mask,2),'the mask is not the same size as the image')

R = mean(RedLayer(mask));
G = mean(GreenLayer(mask));
B = mean(BlueLayer(mask));

RGB = [R G B];
