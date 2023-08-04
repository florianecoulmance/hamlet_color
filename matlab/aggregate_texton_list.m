function patch_list = aggregate_texton_list(folders,kind)


patch_list = [];
for i = 1:numel(folders)
    load(fullfile(folders{i},['all_textons_',kind,'.mat']),'texton_list');
    patch_list = [patch_list;texton_list];
    clear texton_list
    i
end