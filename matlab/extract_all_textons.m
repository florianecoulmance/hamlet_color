function extract_all_textons(imageFolder,file_name,patch_size,resize_factor)

short_name = file_name(1:end-26);

load(fullfile(imageFolder,[short_name,'_oriented_cropped_fish.mat']),'I','mask')

if isPhotoLeftSideOfFish(short_name)
    I = flipud(I);
    mask = flipud(mask);
end

% Find patches for this image
[color_textons,gray_textons,just_color] = find_image_patches(imresize(I,resize_factor),imresize(mask,resize_factor),patch_size);

% Save patches so we don't have to extract them every time
save_image_patches(color_textons,gray_textons,just_color,fullfile(imageFolder,short_name));
