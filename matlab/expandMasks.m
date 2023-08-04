function expandMasks

% Ask the user to point to a folder containing masks
folder = uigetdir();
f = filesep;
dirs = split(folder,f);
savefolder = fullfile(f,dirs{2:end-1},'expanded_files');
if ~exist(savefolder,'dir');mkdir(savefolder);end

% Get the list of mat files that contain masks
matFilesInFolder = dir(fullfile(folder,'*small_chart_masks.mat'));

% If everything is good, loop through all files in folder
for i = 1:numel(matFilesInFolder)
    fileName = matFilesInFolder(i).name;
    shortName = fileName(1:end-22);
    
    % Load,shrink, and save color chart masks
    load(fullfile(folder,fileName),'masks','T','RGB');
    maskFields = fieldnames(masks);

    for j = 1:numel(maskFields)
        masks.(maskFields{j}).mask = smMask2fullMask(masks.(maskFields{j}).small_mask);
        masks.(maskFields{j}) = rmfield(masks.(maskFields{j}),'small_mask');
    end
    
    save(fullfile(savefolder,[shortName,'.mat']),'masks','T','RGB');
end
    






