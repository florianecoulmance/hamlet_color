% Script to mark the body outline of the animal
function fish_mask = make_fish_mask(in_img, mask_save_path)

% Show the image so the user can mark the fish
% Note that we are excluding the fins since they aren't expressed
% consistently.
figure
imshow(in_img)
pix = impoly(gca);
wait(pix);
fish_mask = createMask(pix);

save(mask_save_path,'fish_mask');
close all

