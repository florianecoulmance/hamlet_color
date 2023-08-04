function patch_list = aggregate_patch_list(folders,side)


patch_list = [];
for i = 1:numel(folders)
    load(fullfile(folders{i},['all_textons_',side,'.mat']),'texton_list');
    patch_list = [patch_list;texton_list];
    clear texton_list
    i
end