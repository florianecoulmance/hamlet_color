function patch_list = aggregate_patch_list_grayscale(folders,side)


patch_list = [];
for i = 1:numel(folders)
    load(fullfile(folders{i},['all_textons_',side,'_grayscale.mat']),'texton_list');
    patch_list = [patch_list;texton_list];
    clear texton_list
    i
end