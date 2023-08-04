function save_thumbnails(imageFolder,file_name,resizeFactor)

short_name = file_name(1:end-13);

load(fullfile(imageFolder,file_name),'fish_mask')
load(fullfile(imageFolder,[short_name,'.mat']),'Irgb')


% Rotate the fish horizontally
stats = regionprops(fish_mask,'Orientation');
if stats(1).Orientation>90
    theta = 90-stats(1).Orientation;
elseif stats(1).Orientation<45
    theta = -stats(1).Orientation;
else
    theta = stats(1).Orientation;
end
mask = imresize(imrotate(fish_mask,theta),resizeFactor);
I = imresize(imrotate(Irgb,theta),resizeFactor);

I = cropImageToMask(I,mask);
stats = regionprops(mask,'BoundingBox');
Icrop = imcrop(I,stats.BoundingBox);
mask = imcrop(mask,stats.BoundingBox);
imwrite(Icrop,fullfile(imageFolder,[short_name,'_thumb.jpg']))
imwrite(mask,fullfile(imageFolder,[short_name,'_thumbMask.jpg']))