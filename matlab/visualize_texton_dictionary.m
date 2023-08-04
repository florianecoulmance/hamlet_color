function texton_images = visualize_texton_dictionary(texton_dictonary,num_textons,patch_size)


texton_dictonary = reshape(texton_dictonary,[num_textons, patch_size^2 3]);
texton_images = cell(num_textons,1);
for i = 1:num_textons
    this_texton = zeros(patch_size,patch_size,3);
    for j = 1:3
        this_texton(:,:,j) = reshape(texton_dictonary(i,:,j),[patch_size,patch_size]);
    end
    texton_images{i} = this_texton;
end

