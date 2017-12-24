im = (imread('TP/im/scotland_house.png')); 
imB = (imread('TP/im/scotland_plain.png')); 
nbSuperPixelsWanted = 300;
epsilon =5;
global R;
R = 1;
%this step computes and return the super pixels
%L is the matrix which has the same size as the image, each value
%correspond to the superpixel the pixel belongs
%N is the number of superpixels actually computed
[L,N] = superpixels(im, nbSuperPixelsWanted);
[LB,NB] = superpixels(imB, nbSuperPixelsWanted);
figure;
BW = boundarymask(L);
imshow(imoverlay(im,BW,'cyan'),[])
figure;
BWB = boundarymask(LB);
imshow(imoverlay(imB,BWB,'cyan'),[])


outputImage = zeros(size(im),'like',im);
idx = label2idx(L);
numRows = size(im,1);
numCols = size(im,2);
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
numRowsB = size(imB,1);
numColsB = size(imB,2);
for labelVal = 1:NB
    redIdx = idxB{labelVal};
    greenIdx = idxB{labelVal}+numRowsB*numColsB;
    blueIdx = idxB{labelVal}+2*numRowsB*numColsB;
    outputImageB(redIdx) = mean(imB(redIdx));
    outputImageB(greenIdx) = mean(imB(greenIdx));
    outputImageB(blueIdx) = mean(imB(blueIdx));
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
    alpha=0.8*alpha;
    plot(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA)
end


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