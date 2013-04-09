function [dm,IDX] = subsample(GeoD, neighGraph, k, draw, coords)

% SUBSAMPLE Subsample in the overdense regions given the geodesic distance
% between data points
%
% IDX = subsample(GeoD, neighGraph k, draw, coords);
% 
% Input:
% GeoD = geodesic distances between data points, a full N-by-N matrix,
% where N is the number of data points
% k = number of neighbors for each data point
% draw = option about whether to draw the results. Default is 0
% coords = coordinates projected onto 2D space, 2-by-N matrix
%
% Output:
% IDX = Remaining indices of data points, m-d column vector
% 
% See also:
% ISOMAPII, GENVIRTUAL

%% Check inoput params & Prepare deletion threshold

if nargin <= 3
draw = 0;
end
%%% Test: delete points according to Euclidean distances in projected space
testEud=1;
if testEud
Eud=zeros(size(GeoD));
N = size(GeoD,1);
% derive Euclidean distances from 2D coordinates
for i=1:N
vecX = coords(1,i*ones(1,N));
vecY = coords(2,i*ones(1,N));
Eud(:,i) = sqrt((vecX-coords(1,:)).^2+(vecY-coords(2,:)).^2);
end
GeoD = Eud;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(GeoD,1);
if N~= size(GeoD,2)
error('Distance matrix must be square.');
end
if draw == 1
if nargin < 5
error('Too few arguments.');
end
if N~=size(coords,2) || 2~=size(coords,1)
error('Dimentionality disagrees');
end
end

% find threshold 
Thrh = median(GeoD(logical(neighGraph)));
dm = Thrh;
%Thrh = max(GeoD(:))*(T/n);

%% Search the over dense region and remove data points
IDX = true(N,1);
diagIdx = sub2ind(size(GeoD),1:N,1:N); % get indices along the main diagonal
GeoD(diagIdx) = intmax; % set the largest values to main diagonal
disp(' Sorting geodesic distances...');
[~, Idx] = sort(GeoD(:),1,'ascend');
ThD = GeoD<Thrh; % indication: smaller than threshold
%[I,J] = find(ThD);
[I,J] = ind2sub(size(GeoD),Idx);

disp(' Removing samples...');
count = 0;
for t = 1:length(I)
It = I(t); Jt = J(t); % 
if It~=Jt && It~=-1 && Jt~=-1
countI = sum(ThD(:,It));
countJ = sum(ThD(:,Jt));

    if countI < k || countJ < k % sparse
        continue;
    end

    if countI >countJ
        % delete Ith point
        IDX(It) = 0;
        % delete from index array
        ItReIdx = I == It | J == It;
        I(ItReIdx) = -1; % mark
        J(ItReIdx) = -1;
        % delete from indication array
        ThD(It,:) = 0;
        ThD(:,It) = 0;
    else
        % delete Jth point
        IDX(Jt) = 0;
        % delete from index array
        JtReIdx = I == Jt | J == Jt;
        I(JtReIdx) = -1;
        J(JtReIdx) = -1;
        % delete from indication array
        ThD(Jt,:) = 0;
        ThD(:,Jt) = 0;
    end
    count = count+1;
end
%if count == T % reach the limit
% break;
%end
end
fprintf('%d samples removed.\n',count);
indexFull = 1:N;
IDX = indexFull(IDX);
IDX = IDX';

%% Draw Result
if draw
figure; hold on;
plot(coords(1,:),coords(2,:),'bo');
plot(coords(1,IDX),coords(2,IDX),'ro');
hold off;
end
