function Comparison()
clc;
Ratio = 0.05;
      for num = 1:1       
        i = 1;
        filename = strcat(num2str(num),'.jpg');
        Im = imread(filename);
        MatFileName = strcat('Point',num2str(num),'.mat');
        load(MatFileName); % Point: for selected location coordinates ROI
        
        [f,~,~] = ToShowHulls(Point,Im);%for visualization;
        figure(1);hold on;set(f,'Visible','on');title('ROI');
       
        
        tic;
        [OutCoverage,ROI] = OuterDecompositionROI(Im,Point);        
        TrainSample = ObtainPixel(Im,ROI);
        Time(i,1) = toc;        
        MatFileName = strcat(num2str(num),'.mat');
        load(MatFileName);
        [ReIm,TP(i,1),FP(i,1)] = TrainingSVDD(TrainSample,Im,num,Ratio);
        figure(2);hold on;imshow(ReIm);title('ours');
        
        tic;
        ConvexROI = OutCoverage | ROI;
        ConvexSample = ObtainPixel(Im,ConvexROI);
        Time(i,2) = toc;
        [ConvexReIm,TP(i,2),FP(i,2)] = TrainingSVDD(ConvexSample,Im,num,Ratio); 
        figure(3);hold on;imshow(ConvexReIm);title('convexHull');
        
        tic;
        RectangleRegion = floor( ObtainRectangle(Point));        
        RetangleSample = Im(RectangleRegion(3):RectangleRegion(4),RectangleRegion(1):RectangleRegion(2),:);
        Time(i,3) = toc;
        RetangleSample = RGB2Pixel(RetangleSample);
        [RetangleReIm,TP(i,3),FP(i,3)] = TrainingSVDD(RetangleSample,Im,num,Ratio);       
        figure(4);hold on;imshow(RetangleReIm);title('RectangleHull');
        
        tic;
        [SegmentIm,TP(i,4),FP(i,4)] = ImageOtsuSegment(Im,TrainSample,num);
        Time(i,4) = toc;
        figure(5);hold on;imshow(SegmentIm);title('Otsu');
        
        tic;
        [MattingIm,TP(i,5),FP(i,5)] = ImageMatting(Im,TrainSample,ROI,OutCoverage,num);
        Time(i,5) = toc;
        figure(6);hold on;imshow(MattingIm);title('Matting');
        
        fprintf('%.3f %.3f %.3f %.3f %.3f\n',Time(i,1),Time(i,2),Time(i,3),Time(i,4),Time(i,5));%Time
        Time = TP;
        fprintf('%.3f %.3f %.3f %.3f %.3f\n',Time(i,1),Time(i,2),Time(i,3),Time(i,4),Time(i,5));%TP
        Time = FP;
        fprintf('%.3f %.3f %.3f %.3f %.3f\n',Time(i,1),Time(i,2),Time(i,3),Time(i,4),Time(i,5));%FP
        
      end
end

function [ReIm,TP,FP] = ImageMatting(Im,TrainSample,ROI,OutCoverage,num)


% MatFileName = strcat('Point',num2str(num),'.mat');
% load(MatFileName); % Point: for selected location coordinates ROI


ImPixel = RGB2Pixel(Im);
[row,col,~] = size(Im);
ReIm = zeros(row*col,1);
ReIm(ROI(:)) = 255;
ReIm(OutCoverage(:)) = 128;
ReIm = reshape(ReIm,[row col]);

image = Im;
trimap = uint8(ReIm);
ReIm = double(trimap);
a_cf = closedFormMatting(image, ReIm);

image = double(image);
ReIm = uint8(a_cf.*image);

ReR = logical(a_cf);ReR = ReR(:);
MatFileName = strcat(num2str(num),'.mat');
load(MatFileName);
TP = 100*sum(ReR(Label(:)))/sum(Label(:));
FP = 100*sum(ReR(~Label(:)))/sum(~Label(:));

end


function [R,Thre] = BinarySegment(R)
[counts,~] = imhist(R);
T = otsuthresh(counts);
% BW = imbinarize(R,T);
Thre = floor(T*255); 
R = double(R);
R(R<Thre) = floor((0+Thre)/2);
R(R>=Thre) = floor((255+Thre)/2);
R = uint8(R);
end

function [ReIm,TP,FP] = ImageOtsuSegment(Im,X,num)

        R = Im(:,:,1); G = Im(:,:,2); B = Im(:,:,3);
        [R,RThres] = BinarySegment(R);
        [G,GThres] = BinarySegment(G);
        [B,BThres] = BinarySegment(B);
        ReIm = cat(3,R,G,B);        
        R_Value = floor([RThres/2 (RThres+255)/2]);
        G_Value = floor([GThres/2 (GThres+255)/2]);
        B_Value = floor([BThres/2 (BThres+255)/2]);
        
        str = {'111','112','121','122','211','212','221','222'};
        SegPixelValue = [];
        for i = 1:size(str,2)
           TmpStr = str{1,i};
           Position = [str2double(TmpStr(1)) str2double(TmpStr(2)) str2double(TmpStr(3))];
           SegPixelValue = [SegPixelValue;R_Value(Position(1))  G_Value(Position(2)) B_Value(Position(3))];
        end
        

        M = mean(X);
        Dis = pdist2(SegPixelValue,M);
        [value,posi] = min(Dis);
        SelectRegion = str{1,posi};
        Selected = [SelectRegion(1)-'0'; SelectRegion(2)-'0'; SelectRegion(3)-'0'];
        ReR = GetLogicalMatrix(R,R_Value,Selected,1);
        ReG = GetLogicalMatrix(G,G_Value,Selected,2);
        ReB = GetLogicalMatrix(B,B_Value,Selected,3);
        
        ReR = ReR & ReG & ReB;
        ReR = ReR(:);
        Im2Pix = RGB2Pixel(Im);
        ReIm = zeros(size(Im2Pix));
        ReIm(ReR,:) = Im2Pix(ReR,:);
        ReIm = Pixel2RGB(ReIm,size(Im,1),size(Im,2));
        
        MatFileName = strcat(num2str(num),'.mat');
        load(MatFileName);
        TP = 100*sum(ReR(Label(:)))/sum(Label(:));
        FP = 100*sum(ReR(~Label(:)))/sum(~Label(:));
end


function ReR = GetLogicalMatrix(R,R_Value,Selected,i)
ReR = zeros(size(R));
ReR = logical(ReR);
ReR(R==R_Value(Selected(i))) = 1;
end

function [ReIm,TP,FP] = TrainingSVDD(TrainSample,Im,num,ratio)
 MatFileName = strcat(num2str(num),'.mat');
 load(MatFileName);
 
[row,col,~] = size(Im);
Y = ones(size(TrainSample,1),1);
Mdl = fitcsvm(TrainSample,Y,'KernelScale','auto','OutlierFraction',ratio); 
ImPixel = RGB2Pixel(Im);
[Predict,score] = predict(Mdl,ImPixel); 
Predict = (score>=0);
RePixel = zeros(size(ImPixel));
RePixel(score>=0,:) = ImPixel(score>=0,:);
ReIm = Pixel2RGB(RePixel,row, col);

TP = 100*sum(Predict(Label(:)) )/sum(Label(:));
FP = 100*sum(Predict(~Label(:)))/sum(~Label(:));


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

function [f,ConvexPoint,RectangleRegion] = ToShowHulls(Point,Im)
   Point = [Point;Point(1,:)];
   ConvexHull = convhull(Point);
   ConvexPoint = Point(ConvexHull',:);
   RectangleRegion = ObtainRectangle(Point); %RectangleRegion = [min_x max_x min_y max_y]
   RectangleHull = [RectangleRegion(1) RectangleRegion(3);                 
                 RectangleRegion(1) RectangleRegion(4);
                 RectangleRegion(2) RectangleRegion(4);
                 RectangleRegion(2) RectangleRegion(3)];
   RectanglePoint = [RectangleHull;RectangleHull(1,:)];
   f =  imshow(Im);   
   hold on;
   plot(Point(:,1),Point(:,2),'g-','linewidth',2);
   plot(ConvexPoint(:,1),ConvexPoint(:,2),'m-','linewidth',2);
   plot(RectanglePoint(:,1),RectanglePoint(:,2),'y-','linewidth',2); 
   set(f,'Visible','off');
end

function Sample = ObtainPixel(Im,ROI)
ImPixel = RGB2Pixel(Im);
Sample = ImPixel(ROI(:),:);
end

function RectangleRegion = ObtainRectangle(Point) %RectangleRegion = [min_x max_x min_y max_y]
minx = min(Point(:,1)); maxx = max(Point(:,1));
miny = min(Point(:,2)); maxy = max(Point(:,2));
RectangleRegion = [minx maxx miny maxy];
end

function Im2Pixel = RGB2Pixel(Im)
R = Im(:,:,1);G =Im(:,:,2); B = Im(:,:,3);
R = R(:); G=G(:); B = B(:);
[row,col,~] = size(Im);
[x,y] = meshgrid(1:1:row,1:1:col);
%Im2Pixel = [y(:) x(:)  R G B];
Im2Pixel = double([R G B]);
end

function Pixel2Im = Pixel2RGB(Im2Pixel,row,col)
R = reshape(Im2Pixel(:,1),[row col]);
G = reshape(Im2Pixel(:,2),[row col]);
B = reshape(Im2Pixel(:,3),[row col]);
Pixel2Im = cat(3,R,G,B);
Pixel2Im = uint8(Pixel2Im);
end