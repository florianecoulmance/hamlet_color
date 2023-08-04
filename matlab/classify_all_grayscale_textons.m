function allImages = classify_all_grayscale_textons(folders)

[num,txt,~] = xlsread('/Users/deryaakkaynak/Documents/GitHub/hamlets/image_metadata.xlsx');
imageNames = txt(2:end,13);
speciesNames = txt(2:end,9);
location = txt(2:end,10);
flash = txt(2:end,14);
bestImg = num(:,end);

allImages = struct();
n = 1;
for i = 1:numel(folders)
    this_folder = folders{i};
    files = dir(fullfile(this_folder,'*_grayscale_texton_histogram.mat'));
    for j = 1:numel(files)
        allImages(n).image = files(j).name(1:end-49);
        leftSideFlag = isPhotoLeftSideOfFish(allImages(n).image);
        ind = find(contains(imageNames,allImages(n).image));
        sp = speciesNames(ind);
        
        load(fullfile(this_folder,[files(j).name(1:end-4)]),'N')
        allImages(n).textons = N;
        allImages(n).species = sp{1};
        allImages(n).location = location(ind);
        allImages(n).folder = this_folder;
        allImages(n).leftSideFlag = leftSideFlag;
        if strcmp(flash(ind),'On')
            allImages(n).flash = 1;
        else
            allImages(n).flash = 0;
        end
        allImages(n).bestImg = bestImg(ind(1));
        n = n+1;
    end
    i
end