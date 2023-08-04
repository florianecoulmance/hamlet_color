% Process data sent from Floriane

dataPath = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/IMAGE_PREPROCESS_DATA';

expandMasks;
%% Now, expanded masks should be deployed to relevant folders
% Notes: Derya renamed BOCAS_DEL_TORO in Floriane's data to "BOCAS",
% FLORIDA --> FLKEYS

expandedMasksFolder = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/expanded_files/';
fishMaskFolder = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/IMAGE_PREPROCESS_DATA';
originalsFolder = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/IMAGES_RAW_DATA';
sites = {'BELIZE','BOCAS','FLKEYS','PUERTORICO'};
destinationFolder = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/';

belizeData = dir(fullfile(originalsFolder,sites{1},'CR2','*.CR2'));
bocasData = dir(fullfile(originalsFolder,sites{2},'CR2','*.CR2'));
flkeysData = dir(fullfile(originalsFolder,sites{3},'CR2','*.CR2'));
prData = dir(fullfile(originalsFolder,sites{4},'CR2','*.CR2'));

expandedFiles = dir(fullfile(expandedMasksFolder,'*.mat'));
n = 1;
for i = 1:numel(expandedFiles)
    thisFile = expandedFiles(i).name;
    thisFileShort = thisFile(1:end-4);
    chartMasksFrom = fullfile(expandedMasksFolder,thisFile);
    fishMasksFrom = fullfile(fishMaskFolder,[thisFileShort,'fishMasks.mat']);
    try
        % Is it in Belize?
        if find(contains({belizeData.name},thisFileShort))
            ind = find(contains({belizeData.name},thisFileShort));
            chartMasksTo = fullfile(destinationFolder,'Belize','phenotypes','renamed');
            fishMasksTo = fullfile(destinationFolder,'Belize','phenotypes','renamed');
            
            % Bocas?
        elseif find(contains({bocasData.name},thisFileShort))
            ind = find(contains({bocasData.name},thisFileShort));
            chartMasksTo = fullfile(destinationFolder,'Bocas','phenotypes','renamed');
            fishMasksTo = fullfile(destinationFolder,'Bocas','phenotypes','renamed');
            
            %FL Keys?
        elseif find(contains({flkeysData.name},thisFileShort))
            ind = find(contains({flkeysData.name},thisFileShort));
            chartMasksTo = fullfile(destinationFolder,'FLKeys','renamed');
            fishMasksTo = fullfile(destinationFolder,'FLKeys','renamed');
            
            % Puerto Rico?
        elseif find(contains({prData.name},thisFileShort))
            ind = find(contains({prData.name},thisFileShort));
            chartMasksTo = fullfile(destinationFolder,'Puerto Rico','phenotypes','renamed');
            fishMasksTo = fullfile(destinationFolder,'Puerto Rico','phenotypes','renamed');
        end
        
        movefile(chartMasksFrom,chartMasksTo)
        movefile(fishMasksFrom,fishMasksTo)
    catch err
        errorList{n} = thisFile;
        n = n+1;
    end
    i
end

%% Now that everything is in its place, make and save Irgb for all new images
clear;close;clc

% image folder
folders{1} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Belize/phenotypes/renamed';
folders{2} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Bocas/phenotypes/renamed/';
folders{3} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/FLKeys/renamed/';
folders{4} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Puerto Rico/phenotypes/renamed/';

% color chart
load MacbethColorCheckerData.mat
% chart XYZ
XYZ = getChartXYZvalues(chart,colors);

for k = 2%1:numel(folders)
    imageFolder = folders{k};
    % file list
    files = dir(fullfile(imageFolder,'*.dng'));
    for i = 1%1:numel(files)
        file_name = files(i).name;
        short_name = files(i).name(1:end-4);
        % Open only if mask exists
        if exist(fullfile(imageFolder,[short_name,'.mat']),'file') == 2
            load(fullfile(imageFolder,[short_name,'.mat']),'masks');
            if exist('masks','var')
                save_Irgb(masks,colors,XYZ,file_name,imageFolder);
            end
        end
        clear masks
        fprintf(['Folder ',num2str(k),' , image ',num2str(i),' / ',num2str(numel(files)),' done. /n'])
    end
end

%% Thumbnail comparison

clear;close;clc

% image folder
folders{1} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Belize/phenotypes/renamed';
folders{2} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Bocas/phenotypes/renamed/';
folders{3} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/FLKeys/renamed/';
folders{4} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Puerto Rico/phenotypes/renamed/';

for k = 1%1:numel(folders)
    imageFolder = folders{k};
    % file list
    files = dir(fullfile(imageFolder,'*fishMasks.mat'));
    thumbs = dir(fullfile(imageFolder,'*thumb.jpg'));
    
    for i = 1:numel(files)
        files_new{i} = files(i).name(1:end-13);
    end
    for i = 1:numel(thumbs)
        thumbs_new{i} = thumbs(i).name(1:end-10);
    end
    diffs = find(~ismember(files_new,thumbs_new))
    
end

%%
%Florida

PL17_144gemflo-r11-s4-f4-c2-d1
