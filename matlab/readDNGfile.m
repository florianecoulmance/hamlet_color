function outImg = readDNGfile(path)
%
% OUTIMG = READDNGFILE(PATH)
% 
% PATH      : full or local path that points to the dng file
% OUTIMG    : the .dng file converted to a matrix, values scaled to be
% between 0 and 1
% 
% ************************************************************************
% If you use this code, please cite the following paper:
%
%
% <paper>
%
% ************************************************************************
% For questions, comments and bug reports, please use the interface at
% Matlab Central/ File Exchange. See paper above for details.
% ************************************************************************



warning off MATLAB:tifflib:TIFFReadDirectory:libraryWarning
t = Tiff(path,'r');
offsets = getTag(t,'SubIFD');
setSubDirectory(t,offsets(1));
outImg = read(t);
close(t);
outImg = (double(outImg));