function save_image_patches_grayscale(patches,save_path)

save(fullfile([save_path,'_grayscale_textons.mat']),'patches')