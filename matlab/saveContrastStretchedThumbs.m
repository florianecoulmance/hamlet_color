function saveContrastStretchedThumbs(imageFolder,file_name)

% The original name of the file
short_name = file_name(1:end-13);

% Save thumbnails
I = imread(fullfile(imageFolder,[short_name,'_thumb.jpg']));
I = double(I);
I = uint8(255*(mat2gray(I)));
% Save thumbnails
imwrite(I,fullfile(imageFolder,[short_name,'_thumb_contrast.jpg']))