function classify_image_patches_grayscale(folders,texton_dictionary_left,texton_dictionary_right)

num_textons = size(texton_dictionary_left,1);


for i = 1:numel(folders)
    this_folder = folders{i};
    files = dir(fullfile(this_folder,'*_grayscale_textons.mat'));
    for j = 1:numel(files)
        this_file = files(j).name;
        short_name = this_file(1:end-13);
        leftFlag = isPhotoLeftSideOfFish(short_name);
        load(fullfile(this_folder,files(j).name));
        s = size(patches);
        if leftFlag
            [idx,~] = kmeans(patches,num_textons,'start',texton_dictionary_left);
        else
            [idx,~] = kmeans(patches,num_textons,'start',texton_dictionary_right);
            
        end
        % Texton histogram
        [N,~] = histcounts(idx,'normalization','pdf');
        save(fullfile(this_folder,[files(j).name(1:end-4),'_grayscale_texton_histogram.mat']),'N')
        fprintf([num2str(j),' / ',num2str(numel(files)),'\n']);
    end
    i
end



