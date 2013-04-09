function num = gennewimg(P, FolderPath, PosNameLst, OutputFolder)

% GENNEWIMG Generate new images according to neighbor samples
%
% num = gennewimg(P, FolderPath, PosNamelst, OutputFolder)

%
% Input: 
% P = Struct for neighbor indices and weights returned by GENVIRTUAL
% FolderPath = Folderpath of samples images
% PosNameLst = Cell array storing names of samples
% OutputFolder = Output folder path

%
% Output:
% num = Number of virtual images generated
%
% See also:
% GENVIRTUAL, SUBSAMPLE

%% Mix neighbour images
if ~ischar(FolderPath) || ~ischar(OutputFolder)
error('Must use a string as folder path.');
end

[k,M] = size(P.nidx);
tempImg = imread([FolderPath,'\',PosNameLst{1}]);
m = size(tempImg,1); n = size(tempImg,2);
clear tempImg;

delete([OutputFolder, '*.bmp']);% use functional form 
for i = 1:M
idx = P.nidx(:,i);
weight = P.weight(:,i);
virtualImg = zeros(m,n,'double');
% load img & add
for t = 1:k
% convert to double
neighImg = im2double(imread([FolderPath,'\', PosNameLst{idx(t)} ]));
virtualImg = virtualImg + neighImg.*weight(t);
end
virtualImg = im2uint8(virtualImg);
imwrite(virtualImg,[OutputFolder, '\', num2str(i), '.bmp']);
end

num = M;
fprintf(1,'\t%d virtual images generated\n',num);
