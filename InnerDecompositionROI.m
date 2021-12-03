
function InnerDecompositionROIChild()  
% inner decomposition
clc;
%load Child;
load eyu;
load 43;
load sheep;
% load sheepAuto;
% load PlaneAuto;
load plane
load moon;
load scurve;
%Point = [Point(:,2) Point(:,1)];

Im = imread('Training4.png');
Im = imread('0243.jpg');
Im = imread('43.png');
Im = imread('plane.jpg');
Im = imread('moon.jpg');
Im = imread('scurve.jpg');
%Im = imread('sheep.jpg');

figure(1);hold on;
imshow(Im); 
P = Point;
n = size(Point,1);
DT = delaunayTriangulation(P);
pointIndex = 1:1:n;
Point = [Point;Point(1,:)];
ConcaveSetInd = [];
xv = Point(:,1);yv = Point(:,2);
nn = size(DT.ConnectivityList,1);
for i=1:nn
    query = DT.ConnectivityList(i,:);
    Triangle = Point(query',:);
    center = mean(Triangle);

    xq = center(1);yq=center(2); 
    in = inpolygon(xq,yq,xv,yv);
    if (~in)
        hold on;
        %plot(center(:,2),center(:,1),'m*');
        ConcaveSetInd = [ConcaveSetInd i];
    end
end

FullSet = 1:1:nn;
ConvexSetInd = setdiff(FullSet,ConcaveSetInd);
ConcaveVertex = DT.ConnectivityList(ConcaveSetInd',:);
ConvexVertex = DT.ConnectivityList(ConvexSetInd',:);

%% To show inner triangles
%ToShowDesiredDelaunay(ConvexVertex,Point);
plot(Point(:,1),Point(:,2),'r-','linewidth',1);%ROI
%axis([col1-15 col+5 row1-15 row+15]);
%text(Point(1:end-1,1)-2,Point(1:end-1,2),num2str(pointIndex'),'color','b','fontsize',11);
%hold off;
IMM  = ones(size(Im))*255;
IMM  = uint8(IMM);
figure(2); hold on;
imshow(IMM);
hold on;
plot(Point(:,1),Point(:,2),'r-','linewidth',1);%ROI
VisualizationTriangles(Point,ConvexVertex);%show indexed centers
ToShowDesiredDelaunay(ConvexVertex,Point);


%load Tmp1117;
tic;
[ConvexSet,ConvexRegion]    = GetConvexSet(Point,ConvexVertex);
Time = toc;
[ConvexSet,ConvexRegion] = RemoveBlankCell(ConvexSet,ConvexRegion);
[OptimalConvexSet,OptimalConvexRegion] = ObtainFullRegion(ConvexVertex,ConvexSet,ConvexRegion);
ToShowFinalCover(OptimalConvexSet,Point,OptimalConvexRegion);
ConvexSetCoodinates = ObtainCoodinates(Point,OptimalConvexSet,OptimalConvexRegion); %closed hulls

fprintf('\n Points=%d, base regions= %d, maximum convex set=%d, ConvexSets=%d, Time=%f \n',size(Point,1)-1,size(ConvexVertex,1),length(ConvexRegion{1,1}),size(ConvexSetCoodinates,1),Time);
% ReIm is a logical matrix, those pixels with the values 1 (True) should be selected.
ReIm = ObtainSelectedPixel(Im, ConvexSetCoodinates);

Im2Pixel = RGB2Pixel(Im);
X = Im2Pixel(ReIm(:),:);


%% for visualization
[row,col,~] = size(Im);
Im2Pixel = RGB2Pixel(Im);
ReImage = ones(size(Im2Pixel))*255;
ReImage(ReIm(:),:) = Im2Pixel(ReIm(:),:);
ReImage = Pixel2RGB(ReImage,row,col);
figure(3);hold on;
imshow(uint8(ReImage));


end
function SelectedPixel = ObtainSelectedPixel(Im, ConvexSetCoodinates)
Im = double(Im);
ReIm = zeros(size(Im,1),size(Im,2));
ReIm = logical(ReIm);
% Im2Pixel = RGB2Pixel(Im);
% ReIm = zeros(size(Im2Pixel));
for i = 1:size(ConvexSetCoodinates,1)
    Region = ConvexSetCoodinates{i,1};%convex set
    minX = min(Region(:,1)); maxX = max(Region(:,1));
    minY = min(Region(:,2)); maxY = max(Region(:,2));
    Position = floor([minX maxX minY maxY]); %minimum rectangle region containing convex set
    OriIm = ReIm(Position(3):Position(4),Position(1):Position(2));
    OriIm2Pixel = OriIm(:);    
%     ConvexHull = convhull(Region);
%     ConvexHullPosition = Region(ConvexHull,:);
   [x,y] = meshgrid(Position(1):1:Position(2),Position(3):1:Position(4));
   x = x(:);   y=y(:);
   In = inpolygon(y,x,Region(:,2),Region(:,1));
   OriIm2Pixel(In,:) = true;
   OriIm = reshape(OriIm2Pixel, [size(OriIm,1) size(OriIm,2)]);
   ReIm(Position(3):Position(4),Position(1):Position(2),:) = ReIm(Position(3):Position(4),Position(1):Position(2),:) | OriIm;% logical 'or' operation
   
end
SelectedPixel = ReIm;


end
function Pixel2Im = Pixel2RGB(Im2Pixel,row,col)
R = reshape(Im2Pixel(:,1),[row col]);
G = reshape(Im2Pixel(:,2),[row col]);
B = reshape(Im2Pixel(:,3),[row col]);
Pixel2Im = cat(3,R,G,B);
Pixel2Im = uint8(Pixel2Im);

end

function Im2Pixel = RGB2Pixel(Im)
R = Im(:,:,1);G =Im(:,:,2); B = Im(:,:,3);
R = R(:); G=G(:); B = B(:);
[row,col,~] = size(Im);
[x,y] = meshgrid(1:1:row,1:1:col);
%Im2Pixel = [y(:) x(:)  R G B];
Im2Pixel = double([R G B]);
end
function ConvexSetCoodinates = ObtainCoodinates(Point,OptimalConvexSet,OptimalConvexRegion)
ConvexSetCoodinates = [];
n = size(OptimalConvexSet,1);
k = 0;
for i=1:n
    CurrentLayer = OptimalConvexSet{i,1};
    for ii = 1:size(CurrentLayer,1)
        CurrentRegion = CurrentLayer(ii,:);
        CurrentRegion = sort(CurrentRegion);
        TmpPoint = Point(CurrentRegion',:);
        TmpPoint = [TmpPoint;TmpPoint(1,:)]; %make it closed
        k = k+1;
        ConvexSetCoodinates{k,1}=TmpPoint;
        TmpPoint = [];
    end    
end

end
function ToShowFinalCover(OptimalConvexSet,Point,OptimalConvexRegion)
    for i = 1:size(OptimalConvexSet,1)-1
        V = OptimalConvexSet{i,1};
        for j =1: size(V,1)
            CColor = 'mcgymrgymcgymcgymcgymcgymcgymcgymcgymcgym';
            Vertex = V(j,:);
            TmpVertex = Point(Vertex',:);
            TmpVertex = [TmpVertex;TmpVertex(1,:)];
            hold on;
%             if (i<6)
%                 patch(TmpVertex(:,1),TmpVertex(:,2),CColor(i),'facealpha',0.3);
%             else
              Len = mod(j,length(CColor));
              if (Len==0) Len=1; end
              patch(TmpVertex(:,1),TmpVertex(:,2),CColor(Len),'facealpha',0.3);
%         end
        end
    end
end

function [ConvexSet,ConvexRegion] = ObtainFullRegion(ConvexVertex,ConvexSet,ConvexRegion)
TotalRegion = 1:1:size(ConvexVertex,1);
index = [];
for i = 1:size(ConvexRegion,1)
    Region = ConvexRegion{i,:};
    for j=1:size(Region,1)
        index = union(index,Region(j,:));
    end
end
ind = setdiff(TotalRegion,index');
RestVertex = ConvexVertex(ind',:);
RestRegion = ind';
ConvexSet{end+1,1} = RestVertex;
ConvexRegion{end+1,1} = RestRegion;

end


function [ConvexSet,ConvexRegion] = RemoveBlankCell(ConvexSet,ConvexRegion)

k = 0;
for i=1:size(ConvexRegion,1)
    Region = ConvexRegion{i,1};
    if (size(Region,1)>0)
        k = k+1;
        C{k,1} = ConvexSet{i,1};
        R{k,1} = ConvexRegion{i,1};
    end
end

ConvexSet = C;
ConvexRegion = R;
end

function ToShowDesiredDelaunay(ConvexVertex,Point)
n = size(ConvexVertex,1);% number of delaunay Triangulation
for i = 1:n
    ind = ConvexVertex(i,:);
    TmpSet = Point(ind',:);
    TmpSet = [TmpSet;TmpSet(1,:)];
    hold on
    plot(TmpSet(:,1),TmpSet(:,2),'b-.');
end
end

function [ConvexSet,ConvexRegion]  = GetConvexSet(Point,ConvexVertex)
ConvexSet = [];
n = size(ConvexVertex,1);
TmpConvexVertex = ConvexVertex;TmpConvexRegion = (1:1:n)';
RecursionFlag = 1;
ConvexVertex_bk = ConvexVertex;
LastConvexVertex = [];LastConvexRegion = [];
ConvexRegion = (1:1:n)';
k = 0;
%% Obtain the maximum convex set, may be not unique
while(RecursionFlag)
    %[TmpConvexVertex,TmpConvexRegion] = ConvexRecursion(Point,TmpConvexVertex,TmpConvexRegion);
    [TmpConvexVertex,TmpConvexRegion] = ConvexSetMerge(Point,ConvexVertex,TmpConvexVertex,ConvexRegion,TmpConvexRegion);
    for i = 1:size(TmpConvexRegion)
      TmpConvexRegion(i,:)=sort(TmpConvexRegion(i,:));    
    end
    [C,ia,~] = unique(TmpConvexRegion,'rows'); %
    TmpConvexRegion = C;
    TmpCV = TmpConvexVertex(ia,:);
    TmpConvexVertex = TmpCV;
    
    if isempty(TmpConvexRegion) 
        RecursionFlag = 0;
    end
    if (RecursionFlag)
        k = k+1;
        LastConvexVertex{k,1} = TmpConvexVertex;
        LastConvexRegion{k,1} = TmpConvexRegion;
    end
end

%% Backtracking k>0

ConvexRegion = [];
LastLayerVertex = LastConvexVertex{k,1};
LastLayerRegion = LastConvexRegion{k,1};
[MaxConvexVertex,MaxConvexRegion] = ObtainMaxConvex(Point,ConvexVertex,LastLayerVertex,LastLayerRegion);
ConvexSet{1,1} = MaxConvexVertex;
ConvexRegion{1,1} = MaxConvexRegion;

LastConvexVertex(end) = []; 
LastConvexRegion(end) = [];


[LastConvexVertex, LastConvexRegion] = RemoveRepeatSubset(ConvexSet,ConvexRegion,LastConvexVertex,LastConvexRegion);
k = size(LastConvexVertex,1);

layer = 1;
for i = k:-1:1
    CurrentLayerVertex = LastConvexVertex{i,1};
    CurrentLayerRegion = LastConvexRegion{i,1};
    if ( size(CurrentLayerRegion,1)<1 )
        LastConvexVertex(end) = []; 
        LastConvexRegion(end) = [];
        continue;
    end
    if (size(CurrentLayerRegion,1)==1) 
        layer = layer + 1;
        ConvexSet{layer,1} = CurrentLayerVertex;
        ConvexRegion{layer,1} = CurrentLayerRegion;
    end
    if (size(CurrentLayerRegion,1)>1)
        [MaxConvexVertex,MaxConvexRegion] = ObtainMaxConvex(Point,ConvexVertex,CurrentLayerVertex,CurrentLayerRegion);
        layer = layer + 1;
        ConvexSet{layer,1} = MaxConvexVertex;
        ConvexRegion{layer,1} = MaxConvexRegion;        
    end
    [LastConvexVertex, LastConvexRegion] = RemoveRepeatSubset(ConvexSet,ConvexRegion,LastConvexVertex,LastConvexRegion);
end


end
function [LastConvexVertex, LastConvexRegion] = RemoveRepeatSubset(ConvexSet,ConvexRegion,LastConvexVertex,LastConvexRegion)
    for i = 1:size(ConvexRegion,1) %get cell, same dimension
        C = ConvexSet{i,1};
        R = ConvexRegion{i,1};
        for ii = 1:size(C,1)
            CurrentConvexRegion = R(ii,:);
             for j=1:size(LastConvexRegion,1) %get cell, same dimension
                 RR = LastConvexRegion{j,1};
                 CC = LastConvexVertex{j,1};
                 index = []; %FullSet = 1:1:size(RR,1);
                 for jj= 1:size(CC,1)
                     CandidateRegion = RR(jj,:);
                     SubsetResult = intersect(CurrentConvexRegion,CandidateRegion);
                     if (length(SubsetResult)<1) %disjoint
                         index = union(index,jj);
                     end
                 end
                 RR = RR(index',:);  CC = CC(index',:);
                 LastConvexRegion{j,1} = RR;  LastConvexVertex{j,1} = CC;
             end
        end
    end
end

function [TmpConvexSet,TmpConvexRegion] = CheckDisjoint(TmpConvexSet,TmpConvexRegion,ConvexSet, ConvexRegion)

PreLayerRegion = TmpConvexRegion;
PreLayerVertex = TmpConvexSet;
index = [];
for j = 1:size(ConvexRegion,1)
        CurrentR = ConvexRegion{j,1};
        CurrentV = ConvexSet{j,:};
        if (size(CurrentR,1)<1), continue;end  %empty
       
       
        for jj = 1: size(CurrentR,1) 
        CurrentRegion = CurrentR(jj,:);
        CurrentVertex = CurrentV(jj,:);
        
       
            for ii = 1:size(PreLayerRegion,1)
                CurrentSubRegion = PreLayerRegion(ii,:);
                CurrentSubVertex = PreLayerVertex(ii,:);
                SubsetResult = intersect(CurrentSubRegion,CurrentRegion);
                if (length(SubsetResult)>0 ) % It is permitted for sharing the same regions between the selected regions and its previous layers.
                    index = union(index,ii);
                end

                SubsetResult = intersect(CurrentVertex,CurrentSubVertex);
                if (length(SubsetResult)>2)
                    index = union(index,ii);
                end
            end
       end
end
FullSet = 1:1:size(PreLayerRegion,1);
RestSet = setdiff(FullSet,index);
    if (length(RestSet)>0)
    TmpConvexSet = TmpConvexSet(RestSet',:);
    TmpConvexRegion = TmpConvexRegion(RestSet',:);
    else
        TmpConvexSet = [];
        TmpConvexRegion = [];
    end
end
function [TmpConvexSet,TmpConvexRegion] = ObtainMaxConvex(Point,ConvexVertex,LastLayerVertex,LastLayerRegion)
TmpConvexSet = [];TmpConvexRegion = [];
n = size(LastLayerRegion,1);
DisMatrix = zeros(n,n)-1; 
%DisMatrix = DisMatrix + eye(n)*1e10;
DisjointSet = [];JointSet = [];
DisJnt = 0; Joint = 0;
DisjointIndex = [];

for i = 1:n-1
    for j = i+1:n
        set1 = LastLayerVertex(i,:);
        set2 = LastLayerVertex(j,:);
        TmpSet = intersect(set1,set2);
        RegionVertex1 = Point(set1',:);  m1 = mean(RegionVertex1);
        RegionVertex2 = Point(set2',:);  m2 = mean(RegionVertex2);
        DisMatrix(i,j) = sqrt((m1-m2)*(m1-m2)'); %DisMatrix(j,i) = DisMatrix(i,j);
        if (length(TmpSet) <= 2) % neighor or disjoint, i.e., no same triangles
            DisJnt = DisJnt+1;
            DisjointSet(DisJnt,:) = [i j];
        else
            Joint = Joint+1;
            JointSet(Joint,:) = [i j];
        end
    end
end

if (isempty(DisjointSet)) % no disjoint maximum convex set, return the first one as the maximum convex set
    TmpConvexSet = LastLayerVertex(1,:);
    TmpConvexRegion = LastLayerRegion(1,:);
    return;
end

% if existing multiple maximum convex set,firstly put the longest two into temporary sets named TmpConvexSet and TmpConvexRegion
% tmp 
[x,y]  = find(DisMatrix == max(max(DisMatrix)));  %[3,5]
flag = ismember([x y],DisjointSet,'rows');
if (flag)
    TmpConvexSet = LastLayerVertex([x y]',:);
    TmpConvexRegion = LastLayerRegion([x y]',:);
    %DisMatrix(x,y) = -2;
else
    % only one maximum convex set, return
    TmpConvexSet = LastLayerVertex(x,:);
    TmpConvexRegion = LastLayerRegion(x,:);
    return;
end

ind = ismember(DisjointSet,[x y],'rows');
DisjointSet(ind,:) = [];   % remove [3,5] 

% restore distance matrix (Upper triangle) to full matrix
for i=1:n-1
    DisMatrix(i,i) = inf;
    for j=i+1:n
        DisMatrix(j,i) = DisMatrix(i,j);
    end
end
DisMatrix(5,5) = inf;


Current_Set = [x y];
Full_Set = DisjointSet(:);
Full_Set = (unique(Full_Set))';
Rest_Set = setdiff(Full_Set,Current_Set);
while(~isempty(Rest_Set))
   % obtain the most similar regions to selected ones from un-selected regions.
   row = Current_Set;
   col = Rest_Set;
   TmpDis = DisMatrix(row,col);
   [x,y] = find(TmpDis==min(min(TmpDis)));
   Selected = row(x);
   UnSelected = col(y);
   Rest_Set = setdiff(Rest_Set,UnSelected);% remove unselected one since it has been investiaged
   SelectedPair = [Selected UnSelected];
   SelectedPair = sort(SelectedPair);% ascending order
   
   % Does it appear in disjoint set?
   ind = ismember(DisjointSet,SelectedPair,'rows');
   if (sum(ind)>0) %Yes
       Current_Set = union(Current_Set,UnSelected);
      TmpConvexSet = [TmpConvexSet;LastLayerVertex(UnSelected,:);];
      TmpConvexRegion = [TmpConvexRegion;LastLayerRegion(UnSelected,:)];
   end
end
end
function [TmpConvexSet,TmpConvexRegion] = BackTrack(Point,ConvexVertex,LastLayerVertex,LastLayerRegion)
TmpConvexSet = [];TmpConvexRegion = [];
n = size(LastLayerRegion,1);
DisMatrix = zeros(n,n)-1; 
%DisMatrix = DisMatrix + eye(n)*1e10;
DisjointSet = [];JointSet = [];
DisJnt = 0; Joint = 0;
DisjointIndex = [];

for i = 1:n-1
    for j = i+1:n
        set1 = LastLayerVertex(i,:);
        set2 = LastLayerVertex(j,:);
        TmpSet = intersect(set1,set2);
        RegionVertex1 = Point(set1',:);  m1 = mean(RegionVertex1);
        RegionVertex2 = Point(set2',:);  m2 = mean(RegionVertex2);
        DisMatrix(i,j) = sqrt((m1-m2)*(m1-m2)'); %DisMatrix(j,i) = DisMatrix(i,j);
        if (length(TmpSet) <= 2) % neighor or disjoint, i.e., no same triangles
            DisJnt = DisJnt+1;
            DisjointSet(DisJnt,:) = [i j];
        else
            Joint = Joint+1;
            JointSet(Joint,:) = [i j];
        end
    end
end

if (isempty(DisjointSet)) % no disjoint maximum convex set, return the first one as the maximum convex set
    TmpConvexSet = LastLayerVertex(1,:);
    TmpConvexRegion = LastLayerRegion(1,:);
end

% if existing multiple maximum convex set,firstly put the longest two into temporary sets named TmpConvexSet and TmpConvexRegion
% tmp 
[x,y]  = find(DisMatrix == max(max(DisMatrix)));  %[3,5]
flag = ismember([x y],DisjointSet,'rows');
if (flag)
    TmpConvexSet = LastLayerVertex([x y]',:);
    TmpConvexRegion = LastLayerRegion([x y]',:);
    %DisMatrix(x,y) = -2;
else
    % only one maximum convex set, return
    TmpConvexSet = LastLayerVertex(x,:);
    TmpConvexRegion = LastLayerRegion(x,:);
    return;
end

ind = ismember(DisjointSet,[x y],'rows');
DisjointSet(ind,:) = [];   % remove [3,5] 

% restore distance matrix (Upper triangle) to full matrix
for i=1:n-1
    DisMatrix(i,i) = inf;
    for j=i+1:n
        DisMatrix(j,i) = DisMatrix(i,j);
    end
end
DisMatrix(5,5) = inf;


Current_Set = [x y];
Full_Set = DisjointSet(:);
Full_Set = (unique(Full_Set))';
Rest_Set = setdiff(Full_Set,Current_Set);
while(~isempty(Rest_Set))
   % obtain the most similar regions to selected ones from un-selected regions.
   row = Current_Set;
   col = Rest_Set;
   TmpDis = DisMatrix(row,col);
   [x,y] = find(TmpDis==min(min(TmpDis)));
   Selected = row(x);
   UnSelected = col(y);
   Rest_Set = setdiff(Rest_Set,UnSelected);% remove unselected one since it has been investiaged
   SelectedPair = [Selected UnSelected];
   SelectedPair = sort(SelectedPair);% ascending order
   
   % Does it appear in disjoint set?
   ind = ismember(DisjointSet,SelectedPair,'rows');
   if (sum(ind)>0) %Yes
       Current_Set = union(Current_Set,UnSelected);
      TmpConvexSet = [TmpConvexSet;LastLayerVertex(UnSelected,:);];
      TmpConvexRegion = [TmpConvexRegion;LastLayerRegion(UnSelected,:)];
   end
end


%% for visualization
%    ind = [4];
%     CColor = 'cgym';
%     for i=1:length(ind)
%         Vertex = LastLayerVertex(ind(i),:);
%         TmpVertex = Point(Vertex',:);
%         TmpVertex = [TmpVertex;TmpVertex(1,:)];
%         hold on;
%         patch(TmpVertex(:,1),TmpVertex(:,2),CColor(i),'facealpha',0.3);
%     end
end
function [TmpConvexVertex,TmpConvexRegion] = ConvexSetMerge(Point,ConvexVertex,LastConvexVertex,ConvexRegion,LastConvexRegion)
TmpConvexVertex =[];TmpConvexRegion=[];
k=0;
for i=1:size(ConvexVertex,1)
   for j = 1:size(LastConvexVertex)
       set1 = ConvexVertex(i,:);
       set2 = LastConvexVertex(j,:);
       TmpSet = intersect(set1,set2);
       if (length(TmpSet)==2)
            %NeighborMatrix(i,j)=1;
            set1 = union(set1,set2);  % union neigbor
            TmpSet = Point(set1',:);
             kk =  convhull(TmpSet(:,1),TmpSet(:,2));  % counterclockwise closed points
            if (length(kk)-1 == length(set1))
                k = k+1;
                ConvexVertex2Tri(k,:) = set1;
                Convex2Region(k,:) = [ConvexRegion(i,:)  LastConvexRegion(j,:)];
            end
       end
   end
end   
if k>0
    TmpConvexVertex = ConvexVertex2Tri;
    TmpConvexRegion = Convex2Region;
end
    
end
function [TmpConvexVertex,TmpConvexRegion] = ConvexRecursion(Point,ConvexVertex,ConvexRegion)
n = size(ConvexVertex,1);% number of delaunay Triangulation
TmpConvexVertex = [];TmpConvexRegion = [];
k = 0;
for i=1:n-1
    for j=i+1:n
        set1 = ConvexVertex(i,:);
        set2 = ConvexVertex(j,:);
        TmpSet = intersect(set1,set2);
        if (length(TmpSet)==2)
            %NeighborMatrix(i,j)=1;
            set1 = union(set1,set2);  % union neigbor
            TmpSet = Point(set1',:);
            
            kk =  convhull(TmpSet(:,1),TmpSet(:,2));  % counterclockwise closed points
            if (length(kk)-1 == length(set1))
                k = k+1;
                ConvexVertex2Tri(k,:) = set1;
                Convex2Region(k,:) = [ConvexRegion(i,:)  ConvexRegion(j,:)];
            end
        end
    end
end
    if k>0 
        TmpConvexVertex = ConvexVertex2Tri;
        TmpConvexRegion = Convex2Region;
    end
end

function VisualizationTriangles(Point,ConvexVertex)
n = size(ConvexVertex,1);
for i = 1:n
    TriangleCoordinate = Point(ConvexVertex(i,:)',:);
    center(i,:) = mean(TriangleCoordinate);
end
hold on;
%text(center(:,1),center(:,2),num2str((1:1:n)'));
end


