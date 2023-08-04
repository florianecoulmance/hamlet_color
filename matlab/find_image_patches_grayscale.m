function out = find_image_patches_grayscale(I,mask,s)

if nargin<3
    s = 5; % we wukk always assume a 5 x 5 patch
end

assert(size(I,1) == size(mask,1),'Image and mask are not the same size');
assert(size(I,2) == size(mask,2),'Image and mask are not the same size');


mask_col = (im2col(mask,[s s],'distinct'))';
keep_ind = sum(mask_col,2) == s^2;


out = run_color_channel(rgb2gray(I),keep_ind,s);


end

function ch_col = run_color_channel(ch,keep_ind,s)

ch_col = (im2col(ch,[s s],'distinct'))';
ch_col(~keep_ind,:) = [];
% Center the patch; ie remove mean
ch_col = ch_col - mean(ch_col,2);
% normalize to have L1 norm = 1
ch_col = ch_col./(vecnorm(ch_col',1))';

% normalize
end

