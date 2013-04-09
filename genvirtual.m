function P = genvirtual(coords, neighgraph, geodist, k, dm, gennum, draw)

% GENVIRTUAL Generates virtual samples in sparse regions of embeddings
%
% P = genvirtual(coords, neighgraph, geodist, k, gennum, draw)
%
% Input:
% corrds = coordinates matrix in two-demensional embedding for N points,
% which should be 2-by-N
% neighgraph = N-by-N sparse graph indicating local neighborhood
% connection
% k = the number of neighours used to reconstruct the virtual data point
% geodist = the geodesic distances between all data points
% gennum = number of virtual samples generated. If omitted, the default
% is sqrt(2*N)
% draw = whether to draw the result. Default is 0
%
% Output:
% P = P.loc is a matrix containing M locations that shoulded be amended 
% with virtual data points. 
% P.nidx is a k-by-M matrix that contains the indices of neighbours. 
% P.weight is a k-by-M matrix consisting of weight of each neighbour.
%
% See also:
% ISOMAPII, KDTREESEARCHER, KNNSEARCH

%% Prepare necessary parameters for searching in 2-D embedding
if nargin <= 6
draw = 0;
end

[dim,N] = size(coords);
if dim~=2
error('Manifold embedding must be projected onto 2 dimension.')
end
if round(k)~=k
error('Number of neighbours must be integer.');
end
if N~=size(neighgraph,1) || N~=size(geodist,1)
error('Dimensionality diagrees.');
end
if N~=size(neighgraph,2) || N~=size(geodist,2)
error('Neighborhood graph and Geodesic distances must be square.')
end

% get the searching radius
%%%%%%%%%%%%%%%Use Euclidean distances in 2D to get radius%%%%%%%%%%%%%
%testEud = 1;
%if testEud
% Eud = zeros(size(geodist));
% N = size(geodist,1);
% derive Euclidean distances from 2D coordinates
% for i=1:N
% vecX = coords(1,i*ones(1,N));
% vecY = coords(2,i*ones(1,N));
% Eud(:,i) = sqrt((vecX-coords(1,:)).^2+(vecY-coords(2,:)).^2);
% end
% geodist = Eud;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dm = median(geodist(logical(neighgraph))); % the median of all geodesic distance between neighbor pairs
% construct kd tree with Euclidean distances between points
disp(' Building KD tree...');
kdtree = KDTreeSearcher(coords'); % observation in rows, variables in column

%% Search the 2D projected embedding to find over-sparse holes
mincoord = min(coords,[],2);
maxcoord = max(coords,[],2);
fprintf(1,'\tSearch sparse regions within range [<%0.2g,%0.2g>, <%0.2g,%0.2g>]\n',...
mincoord(1),mincoord(2),maxcoord(1),maxcoord(2));
P = struct('loc',[],'nidx',[],'weight',[]); % construct a struct
if nargin >= 6
numsearch = round(sqrt(gennum));
else
numsearch = round(sqrt(2*N));
end

count = 0;
%flag = 0;
for row = linspace(mincoord(1),maxcoord(1),numsearch*3)
for col = linspace(mincoord(2),maxcoord(2),numsearch*3)
x=round(row);y=round(col);
[idxsl,dists] = knnsearch(kdtree,[x y],'k',k);% for some residuals
% determin whether it is sparse region
numNeigh = sum(dists<=dm); % number of neighbors within the radius (max k)
if numNeigh < k && numNeigh >= floor(k/2) % make sure the sparsity, but not blank area
P.loc = cat(2,[x;y],P.loc);
P.nidx = cat(2,idxsl(:),P.nidx);
% weight is inverse proportional to distances
denorm = sum(1./dists);
weight = 1./dists./denorm;
P.weight = cat(2,weight(:),P.weight);
count = count+1;
%if count == gennum
% flag = 1; % want to break outer loop
% break;
%end
end
end
%if flag
% break;
%end
end

%% Release memory allocated & draw results
clear kdtree;
if draw
figure;
hold on;
plot(coords(1,:),coords(2,:),'ro');
plot(P.loc(1,:),P.loc(2,:),'bx');
% randomly draw some circles
randindex = ceil(1 + rand (1,10)*(size(P.loc,2)-1));
circle(P.loc(1,randindex),P.loc(2,randindex),dm); % draw circle
%for t = 1:size(P.loc,2) % draw line
%X = zeros(1,k*2);
%Y = zeros(1,k*2);
%idc = 1:2:length(X)-1;
%idr = 2:2:length(X);
%X(idc) = P.loc(1,t);
%Y(idc) = P.loc(2,t);
%X(idr) = coords(1,P.nidx(:,t));
%Y(idr) = coords(2,P.nidx(:,t));
%line(X,Y);
%end
hold off;
end

%% circle function
function circle(X,Y,R)
% CIRCLE Draws a circle at with radius R
N=length(X);
if N~=size(Y)
error('Dimensionality disagrees');
end

alpha=0:pi/50:2*pi;
x = repmat(X,length(alpha),1);
y = repmat(Y,length(alpha),1);
coordx = R.*cos(alpha);
coordy = R.*sin(alpha);
coordx = repmat(coordx',1,N) + x;
coordy = repmat(coordy',1,N) + y;

plot(coordx,coordy,'-');
