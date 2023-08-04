function saveOrientedFish(imageFolder,file_name)

% The original name of the file
short_name = file_name(1:end-13);

% Right or left side?
leftFlag = isPhotoLeftSideOfFish(short_name);

% Load the fish mask
load(fullfile(imageFolder,file_name),'fish_mask')

% Load the color-corrected image
load(fullfile(imageFolder,[short_name,'.mat']),'Irgb')

% Rotate the fish horizontally
stats = regionprops(fish_mask,'Orientation');
if stats(1).Orientation>90
    theta = 90-stats(1).Orientation;
elseif stats(1).Orientation<45
    theta = -stats(1).Orientation;
elseif stats(1).Orientation>= 45 & stats(1).Orientation<90
    theta = 90+stats(1).Orientation;
else
    theta = stats(1).Orientation;
end
mask = imrotate(fish_mask,theta);
I = imrotate(Irgb,theta);

% Check again!
% Rotate the fish horizontally
stats = regionprops(mask,'Orientation');
if stats(1).Orientation>5
    theta = -stats(1).Orientation; 
    mask = imrotate(mask,theta);
    I = imrotate(I,theta);
end

% Crop to the smallest rectangle around
stats = regionprops(mask,'BoundingBox','Area');
if numel(stats) > 1
    keepInd = ([stats.Area] == max([stats.Area]));
    stats(~keepInd) = [];
end

I = imcrop(I,stats(1).BoundingBox);
mask = imcrop(mask,stats(1).BoundingBox);

flipFlag = isFishUpsideDown(mask);
if flipFlag & ~leftFlag
    mask = imrotate(mask,180);
    I = imrotate(I,180);
end
if ~flipFlag & leftFlag
    I = imrotate(I,180);
    mask = imrotate(mask,180);
end

% Save the cropped fish and its mask
save(fullfile(imageFolder,[short_name,'_oriented_cropped_fish.mat']),'I','mask');

% Save thumbnails
imwrite(imresize(cropImageToMask(I,mask),0.2),fullfile(imageFolder,[short_name,'_thumb.jpg']))
imwrite(imresize(mask,0.2),fullfile(imageFolder,[short_name,'_thumbMask.jpg']))
