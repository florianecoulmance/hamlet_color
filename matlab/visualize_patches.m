function texton_images = visualize_patches(this_image_patches,patch_size)

s = size(this_image_patches);
texton_images = cell(s(1),1);

for i = 1:s(1)
    this_texton = zeros(patch_size,patch_size,3);
    for j = 1:3
        this_texton(:,:,j) = reshape(this_image_patches(i,:,j),[patch_size,patch_size]);
    end
    texton_images{i} = this_texton;
end

