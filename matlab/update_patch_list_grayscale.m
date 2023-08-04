function texton_list = update_patch_list_grayscale(folders,side)


for i = 1:numel(folders)
    this_folder = folders{i};
    files = dir(fullfile(this_folder,'*_grayscale_textons.mat'));
    texton_list = [];
    for j = 1:numel(files)
        this_file = files(j).name;
        short_name = this_file(1:end-13);
        leftFlag = isPhotoLeftSideOfFish(short_name);
        
        if leftFlag
            y = 'l';
        else
            y = 'r';
        end
        % Is this the right side?
        if strcmp(side,y)
            load(fullfile(this_folder,files(j).name));
            texton_list = [texton_list;patches];
            clear textons
        end
        %fprintf([num2str(j),' / ',num2str(numel(files)),'\n']);
    end
    % Save the textonlist
    save(fullfile(folders{i},['all_textons_',side,'_grayscale.mat']),'texton_list');
    clc
    i
end