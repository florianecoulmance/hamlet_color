function extract_only_color(imageFolder,file_name,patch_size,resize_factor)

short_name = file_name(1:end-26);

load(fullfile(imageFolder,[short_name,'_oriented_cropped_fish.mat']),'I','mask')

% Find patches for this image
this_image_patches = find_color_patches(imresize(I,resize_factor),imresize(mask,resize_factor),patch_size);

% Save patches so we don't have to extract them every time
save_image_patches(this_image_patches,fullfile(imageFolder,short_name));
