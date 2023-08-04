function extract_color_textons(imageFolder,file_name,patch_size,resize_factor)

short_name = file_name(1:end-26);

load(fullfile(imageFolder,[short_name,'_oriented_cropped_fish.mat']),'I','mask')

if isPhotoLeftSideOfFish(short_name)
    I = flipud(I);
    mask = flipud(mask);
end

% Find patches for this image
this_image_patches = find_image_patches(imresize(I,resize_factor),imresize(mask,resize_factor),patch_size);

% Save patches so we don't have to extract them every time
save_image_patches(this_image_patches,fullfile(imageFolder,short_name));
