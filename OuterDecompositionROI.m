
function [OutCoverage,ReIm] = OuterDecompositionROI(Im,Point)  
% outer decomposition

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
       % hold on;
        %plot(center(:,2),center(:,1),'m*');
        ConcaveSetInd = [ConcaveSetInd i];
    end
end

FullSet = 1:1:nn;
ConvexSetInd = setdiff(FullSet,ConcaveSetInd);
ConcaveVertex = DT.ConnectivityList(ConcaveSetInd',:);
ConvexVertex = DT.ConnectivityList(ConvexSetInd',:);


%% To show outer triangles and hull

IMM = ones(size(Im))*255;
% figure(2);hold on;
% imshow(IMM);
% hold on;
dt = DT;
k = convexHull(dt);
% plot(dt.Points(k,1),dt.Points(k,2), 'g-.','linewidth',4); %hull
% plot(Point(:,1),Point(:,2),'r-','linewidth',2);%ROI
% %text(Point(1:end-1,1)-10,Point(1:end-1,2),num2str(pointIndex'),'color','m','fontsize',8);
% %axis([col1-15 col+5 row1-15 row+15]);
% ToShowDesiredDelaunay(ConcaveVertex,Point);
% VisualizationTriangles(Point,ConcaveVertex);%show indexed centers

HullVertex = k;

Num = 0;% the number of convex set

[ConcaveRegion,SingleRegion] = GetConcaveRegion(Point,HullVertex,ConcaveVertex);

for i = 1:length(SingleRegion)
    Region = ConcaveVertex(SingleRegion(i),:);
    TmpCoodinate = Point(Region',:);
    TmpCoodinate = [TmpCoodinate; TmpCoodinate(1,:)  ];
    Num = Num+1;
    ConvexSetCoodinates{Num,1} = TmpCoodinate;%
end
mm = size(ConcaveRegion,1);

maxlayer = 0;ConvexSets = 0;
for i = 2:mm
    TmpRegion = ConcaveRegion{i,1}; 
    if (length(TmpRegion)>maxlayer), maxlayer=length(TmpRegion);end
    TmpConcaveVertex = ConcaveVertex(TmpRegion',:);
    ConcavePoint = GetConcavePoint(Point,TmpConcaveVertex);
    % need triangulation on concave regions
    TmpConcaveVertex = ReTriangulation(ConcavePoint); 
%     ToShowDesiredDelaunay(TmpConcaveVertex,ConcavePoint);
%     VisualizationTriangles(ConcavePoint,TmpConcaveVertex);%show indexed centers
    
    [TmpConvexSet,TmpConvexRegion]    = GetConvexSet(ConcavePoint,TmpConcaveVertex);
    [TmpConvexSet,TmpConvexRegion] = RemoveBlankCell(TmpConvexSet,TmpConvexRegion);
    [OptimalConvexSet,OptimalConvexRegion] = ObtainFullRegion(TmpConcaveVertex,TmpConvexSet,TmpConvexRegion);%best result
    ConvexSets = ConvexSets+size(OptimalConvexSet,1);
    % write convex sets to ConvexSetCoodinates
    for ii = 1:size(OptimalConvexSet,1)
        Region = OptimalConvexSet{ii,1};
        for jj = 1:size(Region,1)
            TmpRegion = Region(jj,:);
            TmpPoint = ConcavePoint(TmpRegion',:);
            ConvexHull = convhull(TmpPoint);
            TmpPoint = TmpPoint(ConvexHull,:);
            Num = Num+1;
            ConvexSetCoodinates{Num,1} = TmpPoint;
        end
    end
%     ToShowFinalCover(OptimalConvexSet,ConcavePoint,OptimalConvexRegion);    
end
fprintf('\n base regions= %d, maximum convex set=%d, ConvexSets=%d, \n',size(ConcaveVertex,1),maxlayer,Num);
ReIm = ObtainSelectedPixel(Im, ConvexSetCoodinates);
OutCoverage = ReIm; % outer coverage

%% visualization for outer coverage
% [row,col,~] = size(Im);
% Im2Pixel = RGB2Pixel(Im);
% ReImage = zeros(size(Im2Pixel));
% ReImage(ReIm(:),:) = Im2Pixel(ReIm(:),:);
% ReImage = Pixel2RGB(ReImage,row,col);
% figure(2);hold on;
% imshow(uint8(ReImage));
% 
% % solve ROI hull
ReIm1 = zeros(size(ReIm));
ReIm1 = logical(ReIm1);
ConvexHull = convhull(Point);
ConvexPoint{1,1} = Point(ConvexHull',:);
ReIm1 = ObtainSelectedPixel(Im, ConvexPoint);%pixels in the hull

ReIm = xor(ReIm,ReIm1);  %XOR for hull and outer coverage
% 
% [row,col,~] = size(Im);
% Im2Pixel = RGB2Pixel(Im);
% ReImage = zeros(size(Im2Pixel));
% ReImage(ReIm(:),:) = Im2Pixel(ReIm(:),:);
% ReImage = Pixel2RGB(ReImage,row,col);
% figure(3);hold on;
% imshow(uint8(ReImage));


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

function TmpConcaveVertex = ReTriangulation(Point)
P = Point;
n = size(Point,1);
DT = delaunayTriangulation(P);

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
        %hold on;       
        ConcaveSetInd = [ConcaveSetInd i];
    end
end

FullSet = 1:1:nn;
ConvexSetInd = setdiff(FullSet,ConcaveSetInd);
ConcaveVertex = DT.ConnectivityList(ConcaveSetInd',:);
ConvexVertex = DT.ConnectivityList(ConvexSetInd',:);
TmpConcaveVertex = ConvexVertex;
end

function ConcavePoint = GetConcavePoint(Point,TmpConcaveVertex)
n = size(TmpConcaveVertex);
ind = [];
for i=1:n
    ind = union(TmpConcaveVertex(i,:),ind);
end
ConcavePoint = Point(ind',:);
end


function [ConcaveRegion,SingleRegion] = GetConcaveRegion(Point,HullVertex,ConcaveVertex)
ConcaveRegion = [];SingleRegion =[];
n = size(ConcaveVertex,1);
m = n;
NeighMatrix = zeros(n,n);%region neighboring

for i=1:n-1
    Region1 = ConcaveVertex(i,:);
    for j=i+1:n
        Region2 = ConcaveVertex(j,:);
        TmpSet = intersect(Region1,Region2);
        if (length(TmpSet)==2)
            NeighMatrix(i,j) = 1;%NeighMatrix(j,i)=1;
        end
    end
end

k = 0;
ConcaveRegion = [];

for i = 1:n 
    tmp = find(NeighMatrix(i,:)==1);
    if (~isempty(tmp))
        k = k+1;
        ConcaveRegion{k,1} = [i tmp];
    end
end

n = size(ConcaveRegion,1);
for i = 1:n-1
    Region1 = ConcaveRegion{i,1};
    if (length(Region1)<1),continue;end
    for j = i+1:n
        Region2 = ConcaveRegion{j,1};
        if (length(Region2)<1),continue;end
        TmpSet = intersect(Region1,Region2);
        if (~isempty(TmpSet))
            ConcaveRegion{i,1} = union(Region1,Region2);
            ConcaveRegion{j,1} = [];
            Region1 = ConcaveRegion{i,1};
        end
    end
end

for i = 1:n-1   %Check merge again
    Region1 = ConcaveRegion{i,1};
    if (length(Region1)<1),continue;end
    for j = i+1:n
        Region2 = ConcaveRegion{j,1};
        if (length(Region2)<1),continue;end
        TmpSet = intersect(Region1,Region2);
        if (~isempty(TmpSet))
            ConcaveRegion{i,1} = union(Region1,Region2);
            ConcaveRegion{j,1} = [];
            Region1 = ConcaveRegion{i,1};
        end
    end
end

n = size(ConcaveRegion,1);% remove empty cells
ind = [];
for i = 1:n
    if(length(ConcaveRegion{i,1})<1)
        ind = [ind i];
    end
end
ConcaveRegion(ind')=[];

% plus single concave regions
n = size(ConcaveRegion,1);
ind = [];
for i=1:n
    Region1 = ConcaveRegion{i,1};
    ind = union(ind,Region1);
end
FullSet = 1:1:m;
RestSet = setdiff(FullSet,ind);
SingleRegion = RestSet;
if (length(RestSet)>1)    
        SingleRegion = RestSet;    
end

end



function ToShowFinalCover(OptimalConvexSet,Point,OptimalConvexRegion)
    for i = 1:size(OptimalConvexSet,1)-1
        V = OptimalConvexSet{i,1};
        for j =1: size(V,1)
            CColor = 'gymcgymcgymcgym';
            Vertex = V(j,:);
            TmpVertex = Point(Vertex',:);
            TmpVertex = [TmpVertex;TmpVertex(1,:)];
            hold on;
%             if (j==2 && i==3)
%                patch(TmpVertex(:,1),TmpVertex(:,2),CColor(4),'facealpha',0.3);
%             else
                patch(TmpVertex(:,1),TmpVertex(:,2),CColor(j),'facealpha',0.3);
%             end
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
C = [];R = [];
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
if (size(LastConvexRegion,1)==1)
    Region = LastConvexRegion{1,1};
    if (size(Region,1)==1)
        ConvexSet = LastConvexVertex;
        ConvexRegion = LastConvexRegion;
        return;
    end
end
ConvexRegion = [];
if (k==0) 
    return;
end
LastLayerVertex = LastConvexVertex{k,1};
LastLayerRegion = LastConvexRegion{k,1};
[MaxConvexVertex,MaxConvexRegion] = ObtainMaxConvex(Point,ConvexVertex,LastLayerVertex,LastLayerRegion);
ConvexSet{1,1} = MaxConvexVertex;
ConvexRegion{1,1} = MaxConvexRegion;

LastConvexVertex(end) = []; 
LastConvexRegion(end) = [];


[LastConvexVertex, LastConvexRegion] = RemoveRepeatSubset(ConvexSet,ConvexRegion,LastConvexVertex,LastConvexRegion);
if isempty((LastConvexVertex))
    return;
end
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
if (isempty(LastConvexVertex))
    return;
end
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
if (n==1)
    TmpConvexSet = LastLayerVertex;
    TmpConvexRegion = LastLayerRegion;
    return;
end
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
text(center(:,1)-2,center(:,2),num2str((1:1:n)'));
%plot(center(:,1),center(:,2),'m.');
end


