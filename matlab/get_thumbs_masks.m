function [thumbs,masks] = get_thumbs_masks(data)

thumbs = cell(numel(data),1);
masks = cell(numel(data),1);
for i = 1:numel(data)
    
    thumbs{i} = imread(fullfile(data(i).folder,[data(i).image,'_thumb_contrast.jpg']));
    masks{i}  = imread(fullfile(data(i).folder,[data(i).image,'_thumbMask.jpg']));
end
