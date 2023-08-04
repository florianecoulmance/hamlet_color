function shrinkMasks

% Ask the user to point to a folder containing masks
folder = uigetdir();
savefolder = fullfile(folder,'small_files');
if ~exist(savefolder,'dir');mkdir(savefolder);end

% Get the list of mat files that contain masks
matFilesInFolder = dir(fullfile(folder,'*.mat'));
fishMasksInFolder = dir(fullfile(folder,'*fishMasks.mat'));

% Check that the number of fish masks equal the number of chart masks
%assert(2*numel(fishMasksInFolder)==numel(matFilesInFolder),'The number of files for fish masks is not the same as the number of files for color chart masks!');

% If everything is good, loop through all files in folder
for i = 1:numel(fishMasksInFolder)
    fishMaskFileName = fishMasksInFolder(i).name;
    shortName = fishMaskFileName(1:end-13);
    
%     % Load,shrink, and save fishMask
%     load(fullfile(folder,fishMaskFileName),'fish_mask');
%     small_fish_mask = fullMask2smMask(fish_mask);
%     save(fullfile(savefolder,[shortName,'_small_fishMask.mat']),'small_fish_mask');

    % Load,shrink, and save color chart masks
    load(fullfile(folder,[shortName,'.mat']),'masks','T','RGB');
    maskFields = fieldnames(masks);
  
    
    for j = 1:numel(maskFields)
        small_mask = fullMask2smMask(masks.(maskFields{j}).mask);
        masks.(maskFields{j}).small_mask = small_mask;
        masks.(maskFields{j}) = rmfield(masks.(maskFields{j}),'mask');
    end
    
    save(fullfile(savefolder,[shortName,'_small_chart_masks.mat']),'masks','T','RGB');
    %save(fullfile(savefolder,[shortName,'_Irgb.mat']),'Irgb');takes too
    %much space so we will not save Irgb. Derya will remake them.
end
    






