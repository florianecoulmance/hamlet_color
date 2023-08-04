function save_Irgb(masks,colors,XYZ,file_name,imageFolder)

short_name = file_name(1:end-4);

% read image
Ior = readDNGfile(fullfile(imageFolder,file_name));
Ior = mat2gray(Ior);
s = size(Ior);
I = imrotate(Ior,-90);

% extract colors
RGB = getChartRGBvalues(I,masks,colors);

% transformation matrix
T = XYZ' * pinv(RGB)';

% save the corrected image for a visual check
Ixyz = reshape((T*reshape(Ior,[s(1)*s(2) 3])')',[s(1) s(2) 3]);
Irgb = XYZ2ProPhoto(Ixyz); % ProPhoto is a wide gamut RGB space that won't clip most colors.

% save the masks and RGB info
save(fullfile(imageFolder,[short_name,'.mat']),'masks','RGB','T','Irgb');
% write the transformed image, half size, for speed
imwrite(imresize(Irgb,0.5),fullfile(imageFolder,[short_name,'.jpg']))

% Finished mask for
fprintf(['Finished mask for image ',short_name,'\n']);

