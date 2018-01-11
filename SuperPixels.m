%im = (imread('TP/im/scotland_house.png')); 
pathImA = 'TP/im/scotland_house.png';
pathImB = 'TP/im/scotland_plain.png';
nameA = strsplit(pathImA,{'/','.'});
nameA = nameA{3};
nameB = strsplit(pathImB,{'/','.'});
nameB = nameB{3};
im = (imread(pathImA));
imB = (imread(pathImB)); 
numRows = size(im,1);
numCols = size(im,2);
numRowsB = size(imB,1);
numColsB = size(imB,2);

nbSuperPixelsWantedA = round(numRows*numCols/500);
nbSuperPixelsWantedB = round(numRowsB*numColsB/500);

epsilon = 3;
global R;
R = 50;
%this step computes and return the super pixels
%L is the matrix which has the same size as the image, each value
%correspond to the superpixel the pixel belongs
%N is the number of superpixels actually computed
[L,N] = superpixels(im, nbSuperPixelsWantedA);
[LB,NB] = superpixels(imB, nbSuperPixelsWantedB);
%figure;
BW = boundarymask(L);
%imshow(imoverlay(im,BW,'cyan'),[])
%figure;
BWB = boundarymask(LB);
%imshow(imoverlay(imB,BWB,'cyan'),[])


outputImage = zeros(size(im),'like',im);
idx = label2idx(L);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = mean(im(redIdx));
    outputImage(greenIdx) = mean(im(greenIdx));
    outputImage(blueIdx) = mean(im(blueIdx));
end    

figure;

imshow(outputImage)

outputImageB = zeros(size(imB),'like',imB);
idxB = label2idx(LB);

meanB = zeros([NB,3]);
for labelVal = 1:NB
    redIdx = idxB{labelVal};
    greenIdx = idxB{labelVal}+numRowsB*numColsB;
    blueIdx = idxB{labelVal}+2*numRowsB*numColsB;
    mr = mean(imB(redIdx));
    mg = mean(imB(greenIdx));
    mb = mean(imB(blueIdx));
    meanB(labelVal,:)=[mr,mg,mb];
    outputImageB(redIdx) = mr;
    outputImageB(greenIdx) = mg;
    outputImageB(blueIdx) = mb;
end    

figure;

imshow(outputImageB)

centr = regionprops(L,'centroid');
g = adjacentRegionsGraph(L);

%nei = neighboors(67,g);
%theta = angleBetweenCentre(centr,83,67);

%argmin =minimumAngle(1,2.1,centr,neighboors(1,g));

SuperPatch = getSuperPatch(centr,R);

centrB = regionprops(LB,'centroid');
gB = adjacentRegionsGraph(LB);

%neiB = neighboors(13,gB);
%thetaB = angleBetweenCentre(centrB,13,1);

%argminB =minimumAngle(101,4.8,centrB,neighboors(101,gB));

SuperPatchB = getSuperPatch(centrB,R);



A = struct;
A.im = im;
A.L = L;
A.graph = neighboors(g);
A.centre = centr;
A.idx=idx;
A.SuperPatchs =SuperPatch;
B = struct;
B.im = imB;
B.L = LB;
B.graph = neighboors(gB);
B.centre = centrB;
B.idx=idxB;
B.SuperPatchs =SuperPatchB;
B.mean = meanB;

global distanceMatrix;
distanceMatrix = zeros([N,NB])-1;
global distanceSuperPixelMatrix;
distanceSuperPixelMatrix = zeros([N,NB])-1;
[matchA,matchB] = InitializeMatching(N,NB);
meanDistance(A,B,matchA)
plot(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA)

alpha =1.0;
for i = 1:12
    [matchA,matchB] = propagationStep(A,B,epsilon,matchA,matchB);
    meanDistance(A,B,matchA)
    [matchA,matchB] = randomSearchStep(A,B,alpha,epsilon,matchA,matchB,10);
    meanDistance(A,B,matchA)
    alpha=0.7*alpha;
    %plot(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA)
end
plot(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA)

%global Q;
Q = cell([1,N]);
%global A_bar;
A_bar = cell([1,N]);
delta_s = 10;
delta_c = 0.1;

for i=1:N
    A_i = idxToCoord(A.idx{i},numRows,numCols,A.im);
    covX_i = cov(A_i(:,1:2));
    covC_i = cov(A_i(:,3:5));
    Qi = blkdiag(delta_s^2*covX_i,delta_c^2*covC_i);
    Qi = inv(Qi);
    Q{i} = Qi;
    Abar = mean(A_i);
    A_bar{i} = Abar;
end



TransformedImage = zeros(size(A.im));
% for p = 1:(numRows*numCols)
%    x = ceil(p/numRows);
%    y = mod(p-1,numRows)+1;
%    TransformedImage(y,x,:) = A_t(p,A,B,matchA,Q,A_bar) ;
% end

TransformedImage1 = zeros([numRows,numCols]);
TransformedImage2 = zeros([numRows,numCols]);
TransformedImage3 = zeros([numRows,numCols]);

parfor p = 1:(numRows*numCols)
   a = A_t(p,A,B,matchA,Q,A_bar);
   TransformedImage1(p) = a(1);
   TransformedImage2(p) = a(2);
   TransformedImage3(p) = a(3);
end

TransformedImage(:,:,1)= TransformedImage1;
TransformedImage(:,:,2)= TransformedImage2;
TransformedImage(:,:,3)= TransformedImage3;

figure;imshow(TransformedImage/255,[]);
Final = regrain(double(im)/255,TransformedImage/255);
figure;imshow(Final,[]);

%imwrite(TransformedImage/255, strcat('Results/',nameB,'_on_',nameA,'_eps_',num2str(epsilon),'_R_',num2str(R),'.png'));
%imwrite(Final, strcat('Results/',nameB,'_on_',nameA,'_eps_',num2str(epsilon),'_R_',num2str(R),'_regrain.png'));

function plot(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA)
    for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    redIdxM = idxB{matchA(labelVal)};
    greenIdxM = idxB{matchA(labelVal)}+numRowsB*numColsB;
    blueIdxM = idxB{matchA(labelVal)}+2*numRowsB*numColsB;
    outputImage(redIdx) = mean(imB(redIdxM));
    outputImage(greenIdx) = mean(imB(greenIdxM));
    outputImage(blueIdx) = mean(imB(blueIdxM));
    end
    figure; imshow(outputImage);
end

function meanDist = meanDistance(A,B,matchA)
    [~,n] = size(matchA);
    dist = 0;
    for i = 1:n
       dist = dist + CheckDistance(i,matchA(i),A,B); 
    end
    meanDist = dist/n;
end

