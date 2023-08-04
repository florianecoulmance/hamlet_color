function update_texton_list(folders,texton_kind)


for i = 1:numel(folders)
    this_folder = folders{i};
    files = dir(fullfile(this_folder,'*_patches.mat'));
    texton_list = [];
    for j = 1:numel(files)
        this_file = files(j).name;
        load(fullfile(this_folder,this_file),texton_kind);
        texton_list = [texton_list;eval(texton_kind)];
        
        %fprintf([num2str(j),' / ',num2str(numel(files)),'\n']);
    end
    % Save the textonlist
    save(fullfile(folders{i},['all_textons_',texton_kind,'.mat']),'texton_list');
    clc
    i
end