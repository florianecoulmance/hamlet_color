function out = find_color_patches(I,mask,s)

if nargin<3
    s = 5; % we will always assume a 5 x 5 patch
end

assert(size(I,1) == size(mask,1),'Image and mask are not the same size');
assert(size(I,2) == size(mask,2),'Image and mask are not the same size');

% Rotating the fish is not necessary
mask_col = (im2col(mask,[s s],'distinct'))';
keep_ind = sum(mask_col,2) == s^2;

for j = 1:3
    out(:,:,j) = run_color_channel(I(:,:,j),keep_ind,s);
end

out = squeeze(out);

end

function ch_col = run_color_channel(ch,keep_ind,s)

ch_col = (im2col(ch,[s s],'distinct'))';
ch_col(~keep_ind,:) = [];
% get mean
ch_col = mean(ch_col,2);
end

