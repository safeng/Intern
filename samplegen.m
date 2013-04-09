function samplegen(NameListPath, InputFolderPath, numLimit)

% SAMPLEGEN performs sample resampling through mapping the manifold
% embedding onto 2D reduced space. In this 2D reduced space, it 
% removes sample points in overdense regions and adding virtual 
% samples in over sparse regions.
%
% samplegen(NameListPath, InputFolderPath, numLimits)
% 
% Input:
% NameListPath = path of list of sample file names
% InputFolderPath = path of sample data
% numLimit = upper limit of samples loaded. By default, there is no limit
%
% See also:
% PREPROCESSIMGDATA, ISOMAPII, SUBSAMPLE, GENVIRTUAL

clc; close all;
disp('Loading name list...');
fullList = importdata(NameListPath);
if nargin == 3
fullList = fullList(1:numLimit);
end
numTotalImg = length(fullList);

totalRound = ceil(numTotalImg/10000); % loops needed
OutputFolderPath = cell(totalRound,1);
disp('Preparing folders...');
% prepare output folders
for i = 1:totalRound
folderName = sprintf('virtual sample %d',i);
OutputFolderPath{i} = folderName;
if ~exist(folderName,'dir')
mkdir(folderName); % store newly generated samples
else
delete([folderName,'*.bmp']);
end
end
remainFolderPath = 'remaining samples';
if ~exist(remainFolderPath,'dir')
mkdir(remainFolderPath); % store remaining samples
else
delete([remainFolderPath,'*']);
end

curRound = 1;
global X;
options.dims = 2; % reduced to TWO dimension
options.display = 0; % display function in isomap
options.draw = 0; % for drawing function in subsample and resample
k = 12; % k nearest neighbors
disp('Start Resampling...');
for start=1:10000:numTotalImg % 10000 samples each round
close all;
if start + 10000 > numTotalImg
num = numTotalImg - start + 1;
else
num = 10000;
end

fprintf('\n*****************************\nROUND %d\n*****************************\n',curRound);
% preprocess face data and generate partial name list
[facedata, nameList] = preprocessimgdata(fullList, InputFolderPath, start, num);
numImg = size(facedata,2);

X = facedata;
% Perform isomap and get geodesic distances
[Y, ~, E, GeoD] = IsomapII('dfun', 'k', k, options);
clear facedata;
coords = Y.coords{1};
E = E(Y.index,Y.index); % may be more than one connected componet, E and namelist could contain dirty data
nameList = nameList(Y.index);
clear Y;
% Subsample over the sparse regions of embedding space
T = round(numImg/6);
disp('Subsampling...');
[dm,IDX] = subsample(GeoD, E, k, options.draw, coords);
GeoD = GeoD(IDX,IDX);
coords = coords(:,IDX);
E = E(IDX,IDX);
% Generate virtual samples
numGen = T; % number of virtual images to be generated
disp('Resampling...');
P = genvirtual(coords, E, GeoD, k, dm, numGen, options.draw);
% Store generated virtual images into disk
disp('Generating virtual samples...');
gennewimg(P, InputFolderPath, nameList, OutputFolderPath{curRound});
% Store remaining images into disk
disp('Saving remaining samples...')
for iremain = 1:length(IDX)
    ind = IDX(iremain);
    sample = imread([InputFolderPath,'\',nameList{ind}]);
    imwrite(sample,[remainFolderPath,'\',nameList{ind}]);
end
curRound = curRound + 1;
end
