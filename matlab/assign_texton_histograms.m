function assign_texton_histograms(folders,texton_dictionary,texton_kind)

num_textons = size(texton_dictionary,1);


for i = 1:numel(folders)
    this_folder = folders{i};
    files = dir(fullfile(this_folder,'*_patches.mat'));
    for j = 1:numel(files)
        this_file = files(j).name;
        load(fullfile(this_folder,this_file),texton_kind);
        s = size(eval(texton_kind));
        if numel(s) == 2
            s(3) = 1;
        end
        [idx,~] = kmeans(reshape(eval(texton_kind),[s(1) s(2)*s(3)]),num_textons,'start',texton_dictionary);
        % Texton histogram
        [N,~] = histcounts(idx,'normalization','pdf');
        save(fullfile(this_folder,[files(j).name(1:end-numel('_patches.mat')),texton_kind,'_histogram.mat']),'N')
        fprintf([num2str(j),' / ',num2str(numel(files)),'\n']);
    end
    i
end



