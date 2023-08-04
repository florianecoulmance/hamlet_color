% % Pattern analysis workflow for hamlets, with Oscar + Kosmas
% %
% % Derya Akkaynak 2018 | Questions to: derya.akkaynak@gmail.com
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Workflow Summary:
% %
% % STEP 1 - Mark the patches in the color chart, so colors can be
% % standardized. This saves color balanced image, and relevant parameters in
% % a folder.
% %
% % STEP 2 - This reads the color balanced images, and allows the user to mark the boundaries of the fish, so we can analyze the
% % colors/pattern on just its body.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% %% STEP 1: Make color chart masks
% clear;close;clc
%
% % image folder containing dngs, change this to the path that works for you.
% % Color chart masks and data are saved in the same folder.
% imageFolder = fullfile('.','testData');
%
% % file list of dngs
% files = dir([imageFolder,'/*.dng']);
%
% % color chart
% load MacbethColorCheckerData.mat
% % chart XYZ
% XYZ = getChartXYZvalues(chart,colors);
%
% for i = 1:numel(files)
%
%     % Do masks exist?
%     if ~exist(fullfile(imageFolder,[files(i).name(1:end-4),'.mat']),'file')
%
%         % read image
%         Ior = readDNGfile(fullfile(imageFolder,files(i).name));
%         Ior = double(Ior)./2^16;
%         s = size(Ior);
%         I = imrotate(Ior,-90);
%         % make masks
%         masks = makeChartMask(I*3,chart,colors,150);
%         % extract colors
%         RGB = getChartRGBvalues(I,masks,colors);
%
%         % transformation matrix
%         T = XYZ' * pinv(RGB)';
%
%         % save the corrected image for a visual check
%         Ixyz = reshape((T*reshape(Ior,[s(1)*s(2) 3])')',[s(1) s(2) 3]);
%         Irgb = XYZ2ProPhoto(Ixyz); % ProPhoto is a wide gamut RGB space that won't clip most colors.
%
%         % save the masks and RGB info
%         save(fullfile(imageFolder,[files(i).name(1:end-4),'.mat']),'masks','RGB','T','Irgb');
%         % write the transformed image, half size, for speed
%         imwrite(imresize(Irgb,0.5),fullfile(imageFolder,[files(i).name(1:end-4),'.jpg']))
%
%         % Finished mask for
%         fprintf(['Finished mask for image #',num2str(i),'.\n']);
%     else
%         fprintf(['Masks exist for image #',num2str(i),'.\n']);
%     end
% end
%
% %% STEP 2 : Make fish masks for images that have been converted
%
% clear;close;clc
%
% % Data folder. This is where
% imageFolder = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Belize/phenotypes/renamed';
%
% % file list
% files = dir([imageFolder,'/*.dng']);
%
%
% for i = 1:numel(files)
%     file_name = files(i).name(1:end-4);
%     % Do masks exist?
%     if exist(fullfile(imageFolder,[file_name,'.mat']),'file')
%         if ~exist(fullfile(imageFolder,[file_name,'fishMasks.mat']),'file')
%             load(fullfile(imageFolder,[file_name,'.mat']))
%             fish_mask = make_fish_mask(Irgb*2, fullfile(imageFolder,[file_name,'fishMasks.mat']));
%         end
%     end
%     i
% end
% %% Save a version of the IRGB, and its fish mask, already rotated horizontally
% clear;close;clc
% folders{1} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Belize/phenotypes/renamed';
% folders{2} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Bocas/phenotypes/renamed/';
% folders{3} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/FLKeys/renamed/';
% folders{4} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Puerto Rico/phenotypes/renamed/';
%
%
% for k = 1:numel(folders)
%     imageFolder = folders{k};
%     files = dir(fullfile(imageFolder,'*fishMasks.mat'));
%     for i = 1:numel(files)
%         file_name = files(i).name;
%         saveOrientedFish(imageFolder,file_name)
%         i
%     end
% end
%
% %% Contrast stretch thumbnails
% clear;close;clc
% folders{1} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Belize/phenotypes/renamed';
% folders{2} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Bocas/phenotypes/renamed/';
% folders{3} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/FLKeys/renamed/';
% folders{4} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Puerto Rico/phenotypes/renamed/';
%
%
% for k = 1:numel(folders)
%     imageFolder = folders{k};
%     files = dir(fullfile(imageFolder,'*fishMasks.mat'));
%     for i = 1:numel(files)
%         file_name = files(i).name;
%         saveContrastStretchedThumbs(imageFolder,file_name);
%         i
%     end
% end


%% STEP 3: Extract patches

clear;close;clc
folders{1} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Belize/phenotypes/renamed';
folders{2} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Bocas/phenotypes/renamed/';
folders{3} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/FLKeys/renamed/';
folders{4} = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Puerto Rico/phenotypes/renamed/';

% We will extract a texton of size patch_size x patch_size
patch_size = 5;
resize_factor = 0.1;

for k = 1:numel(folders)
    imageFolder = folders{k};
    % file list
    files = dir(fullfile(imageFolder,'*_oriented_cropped_fish.mat'));
    for i = 1:numel(files)
        file_name = files(i).name;
        extract_all_textons(imageFolder,file_name,patch_size,resize_factor)
        %extract_color_textons(imageFolder,file_name,patch_size,resize_factor)
        %extract_grayscale_textons(imageFolder,file_name,patch_size,resize_factor)
        i
    end
    clc;k
end

%% Gather all textons

rmpath(genpath('/Users/deryaakkaynak/Documents/Research/DownloadedCode/Matlab/vlfeat-0.9.21/'))

% Possible extensions are 'color_textons', gray_textons,just_color
texton_opts(1).class = 'color_textons'; texton_opts(1).numTextons = 50;
texton_opts(2).class= 'gray_textons';texton_opts(2).numTextons = 50;
texton_opts(3).class = 'just_color';texton_opts(3).numTextons = 10;

texton_dictionary = cell(numel(texton_opts),1);
allImages = cell(numel(texton_opts),1);
for i = 1:numel(texton_opts)
    update_texton_list(folders,texton_opts(i).class)
    texton_list = aggregate_texton_list(folders,texton_opts(i).class);
    s = size(texton_list);
    if numel(s) == 2; s(3) = 1; end
    texton_dictionary{i} = update_texton_dictionary(reshape(texton_list,[s(1) s(2)*s(3)]),texton_opts(i).numTextons);
    assign_texton_histograms(folders,texton_dictionary{i},texton_opts(i).class);
    allImages{i} = gather_textons_each_fish(folders,texton_opts(i).class);
end

save allHamletData allImages texton_dictionary texton_opts

%%
% Applies to all
allspecies = {color_textons.species}';
uniques = unique(allspecies);
colors = jet(numel(uniques));

% ALL DATA
for i = 1:3
    data = allImages{i};
    [Y,thesespecies] = get_subset(data,1:numel(data),texton_opts(i).numTextons,'species');
    title = texton_opts(i).class;
    title = [title,'all_data'];
    scatterplot_3D(Y,thesespecies,uniques,colors)
    set(gca,'view',[75.0704   38.7875])
    export_fig(fullfile('.','Figures',[title,'_3D.jpg']),'-transparent')
    close
    
    scatterplot_2D(Y,thesespecies,colors)
    export_fig(fullfile('.','Figures',[title,'_2D.jpg']),'-transparent')
    close
    
    
end

%% NO FLASH
for i = 1:3
    data = allImages{i};
    keepInd = [data.flash]==0;
    [Y,thesespecies] = get_subset(data,keepInd,texton_opts(i).numTextons,'species');
    title = texton_opts(i).class;
    title = [title,'no_flash'];
    scatterplot_3D(Y,thesespecies,uniques,colors)
    set(gca,'view',[75.0704   38.7875])
    export_fig(fullfile('.','Figures',[title,'3D.jpg']),'-transparent')
    close
    
    
    scatterplot_2D(Y,thesespecies,colors)
    export_fig(fullfile('.','Figures',[title,'_2D.jpg']),'-transparent')
    close
end

%% Thumbnail Scatters
data = allImages{1};
[thumbs_all,masks_all] = get_thumbs_masks(data);
keepInd = [data.flash]==0;

% NATURAL LIGHT
for i = 1:3
    data = allImages{i};
    
    [Y,thesespecies] = get_subset(data,keepInd,texton_opts(i).numTextons,'species');
    these_thumbs = thumbs_all(keepInd);
    these_masks = masks_all(keepInd);
    
    title = texton_opts(i).class;
    title = [title,'thumbs_no_flash'];
    
    thumbnailScatter_withThumbArray(Y(:,1),Y(:,2),these_thumbs,these_masks)
    export_fig(fullfile('.','Figures',[title,'Dim1_2_thumbs.jpg']),'-transparent')
    close
    
    thumbnailScatter_withThumbArray(Y(:,2),Y(:,3),these_thumbs,these_masks)
    export_fig(fullfile('.','Figures',[title,'Dim2_3_thumbs.jpg']),'-transparent')
    close
    
    thumbnailScatter_withThumbArray(Y(:,1),Y(:,3),these_thumbs,these_masks)
    export_fig(fullfile('.','Figures',[title,'Dim1_3_thumbs.jpg']),'-transparent')
    close
    
end

%% Location effects

data = allImages{1};
% ALL DATA
alllocations = {data.location}';
uniquelocations = unique(alllocations);
locspe = {'uni','nig','pue'};
colors = hsv(numel(uniquelocations));

keepInd = find(ismember({data.species},'uni'));
keepImages = data(keepInd);
all_N = reshape([keepImages.textons],[texton_opts(1).numTextons numel(keepImages)])';

D = pdist2(all_N,all_N,@distChiSq);
Y = mdscale(D,2);
thesespecies = {keepImages.location};
gscatter(Y(:,1),Y(:,2),theselocs',[],'.',30)


%% Location with 3 species

locspe = {'uni','nig','pue','ind'};
for i = 4%1:numel(locspe)
    for j = 1:numel(allImages)
        data = allImages{j};
        keepInd = find(ismember({data.species},locspe{i}));
        keepImages = data(keepInd);
        theselocs = {keepImages.location};
        all_N = reshape([keepImages.textons],[texton_opts(j).numTextons numel(keepImages)])';
        D = pdist2(all_N,all_N,@distChiSq);
        Y = mdscale(D,2);
        titleWord = [locspe{i},'_',texton_opts(j).class];
        location_scatterplot(Y,theselocs,locspe{i});
        
        export_fig(fullfile('.','Figures',[titleWord,'Dim1_Dim2.jpg']),'-transparent')
        close
    end
end

%% Thumbs with 3 species

data = allImages{1};
[thumbs_all,masks_all] = get_thumbs_masks(data);
keepInd = [data.flash]==0;

locspe = {'uni','nig','pue','ind'};
for i = 1:numel(locspe)
    for j = 1:numel(allImages)
        data = allImages{j};
        keepInd = find(ismember({data.species},locspe{i}));
        keepImages = data(keepInd);
        theselocs = {keepImages.location};
        all_N = reshape([keepImages.textons],[texton_opts(j).numTextons numel(keepImages)])';
        D = pdist2(all_N,all_N,@distChiSq);
        Y = mdscale(D,2);
        titleWord = [locspe{i},'_',texton_opts(j).class];
        keep_thumbs = thumbs_all(keepInd);
        keep_masks = masks_all(keepInd);
        thumbnailScatter_withThumbArray(Y(:,1),Y(:,2),keep_thumbs,keep_masks)
        
        export_fig(fullfile('.','Figures',[titleWord,'Dim1_Dim2_thumbs.jpg']),'-transparent')
        close
    end
end

%% Just look and ind and pue

% ALL DATA
markers = {'o','^'};
colors = jet(numel(uniques));
for i = 1
    data = allImages{i};
    [Y,thesespecies] = get_subset(data,1:numel(data),texton_opts(i).numTextons,'species');
    titleWord = texton_opts(i).class;
    
    for k = 1:numel(uniques)
        titleWord = [titleWord,'all_data_markers','_',uniques{k}];
        for j = 1:numel(data)
            if (strcmp(uniques{k},data(j).species))
                if data(j).flash==1
                    p(j) = plot3(Y(j,1),Y(j,2),Y(j,3),'marker',markers{1},'markerfacecolor','r',...
                        'markersize',10,'markeredgecolor','k');
                    hold on
                else
                    p(j) = plot3(Y(j,1),Y(j,2),Y(j,3),'marker',markers{2},'markerfacecolor','k',...
                        'markersize',10,'markeredgecolor','k');
                    hold on
                end
            end
        end
        axis square
        set(gca,'fontsize',40)
        set(gca,'fontname','Arial')
        xlabel('MDS Dim 1')
        ylabel('MDS Dim 2')
        zlabel('MDS Dim 3')
        set(gcf,'units','normalized')
        set(gcf,'outerposition',[0 0 1 1])
        set(gca,'xlim',[-0.27 0.29])
        set(gca,'xlim',[-0.15 0.17])
        set(gca,'xlim',[-0.14 0.16])
        title(uniques{k})
        export_fig(fullfile('.','Figures',[titleWord,'.jpg']),'-transparent')
        close
    end
    
end

%% Only pue, color textons

data = allImages{1};
ind_flash = find(ismember({data.species},'ind') & [data.flash]==1);
ind_noflash = find(ismember({data.species},'ind') & [data.flash]==0);
colors = jet(3);
[Y_flash,loc_flash] = get_subset(data,ind_flash,texton_opts(1).numTextons,'location');
[Y_noflash,loc_noflash] = get_subset(data,ind_noflash,texton_opts(1).numTextons,'location');

p1 = gscatter(Y_flash(:,1),Y_flash(:,2),loc_flash',[],'*',10);
hold on
p2 = gscatter(Y_noflash(:,1),Y_noflash(:,2),loc_noflash',[],'.',40);
legend([p1(1) p1(2) p2(1)],{'bel','por','bel'})
axis square
set(gca,'fontsize',40)
set(gca,'fontname','Arial')
set(gcf,'units','normalized')
set(gcf,'outerposition',[0 0 1 1])
title('ind')



