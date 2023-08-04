% Count images that have not been processed yet
% (All images processed are from Belize, so only searching in this folder.)


folder = '/Volumes/DA_Sept18Exp/Hamlets from Kosmas/Belize/phenotypes/renamed';
rawFiles = dir(fullfile(folder,'*.CR2'));
filesList = [];
todo = [];
n = 1;
m = 1;
fid = fopen('filesToRun.txt','a+');
fopen(fid);
for i = 1:numel(rawFiles)
    if exist(fullfile(folder,[rawFiles(i).name(1:end-4),'.mat']))
        filesList{n} = [rawFiles(i).name(1:end-4),'.mat'];
        n = n + 1;
    else
        todo{m,1} = [rawFiles(i).name(1:end-4),'.CR2'];
        fprintf(fid,[todo{m,1},'\n']);
        m = m+1;
    end
    
end
fclose(fid);