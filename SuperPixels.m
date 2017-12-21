im = (imread('TP/im/simpson512.png')); 
nbSuperPixelsWanted = 100;
radius = 70;
%this step computes and return the super pixels
%L is the matrix which has the same size as the image, each value
%correspond to the superpixel the pixel belongs
%N is the number of superpixels actually computed
[L,N] = superpixels(im, nbSuperPixelsWanted);
im(28,34,:) = 255;
im(55,79,:) = 255;
im(94,99,:) = 255;
figure;
BW = boundarymask(L);
imshow(imoverlay(im,BW,'cyan'),[])


outputImage = zeros(size(im),'like',im);
idx = label2idx(L);
numRows = size(im,1);
numCols = size(im,2);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;1
    outputImage(redIdx) = mean(im(redIdx));
    outputImage(greenIdx) = mean(im(greenIdx));
    outputImage(blueIdx) = mean(im(blueIdx));
end    

figure;

imshow(outputImage)

centr = regionprops(L,'centroid');
g = adjacentRegionsGraph(L);

nei = neighboors(13,g);
%theta = angleBetweenCentre(centr,13,2);

argmin =minimumAngle(13,2,centr,nei);

SuperPatch = getSuperPatch(centr,radius);

outputImage = im;
idx = label2idx(L);
numRows = size(im,1);
numCols = size(im,2);

listLabelVal = [33, 65]
for labelVal = listLabelVal
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = 255;
    outputImage(greenIdx) = 255;
    outputImage(blueIdx) = 255;
end
figure;
imshow(outputImage)

labelSuperPixelCentralA = 33;
labelSuperPixelCentralB = 65;
superPatchA = SuperPatch{labelSuperPixelCentralA};
superPatchB = SuperPatch{labelSuperPixelCentralB};
centerA = centr;
centerB = centr;
imageA=im;
imageB = im;
idxA = idx;
idxB = idx;
sigma21 = 0.01;
sigma22 = 0.01;
res = distanceSuperPatchL2(superPatchA, superPatchB, labelSuperPixelCentralA ,...
   labelSuperPixelCentralB, centerA, centerB, imageA, imageB, idxA, idxB);
res

labelSuperPixelCentralA = 3 ;
labelSuperPixelCentralB = 20;
superPatchA = SuperPatch{labelSuperPixelCentralA};
superPatchB = SuperPatch{labelSuperPixelCentralB};
centerA = centr;
centerB = centr;
imageA=im;
imageB = im;
idxA = idx;
idxB = idx;
sigma21 = 0.01;
sigma22 = 0.01;
res = distanceSuperPatchL2(superPatchA, superPatchB, labelSuperPixelCentralA ,...
   labelSuperPixelCentralB, centerA, centerB, imageA, imageB, idxA, idxB);
res

[matchA,matchB] = InitializeMatching(99,99);
[matchA,matchB] = propagationStep(imageA,imageB,g,g,centerA,centerB,SuperPatch,SuperPatch,idxA,idxB,matchA,matchB,1);
[matchA,matchB] = propagationStep(imageA,imageB,g,g,centerA,centerB,SuperPatch,SuperPatch,idxA,idxB,matchA,matchB,1);
[matchA,matchB] = propagationStep(imageA,imageB,g,g,centerA,centerB,SuperPatch,SuperPatch,idxA,idxB,matchA,matchB,1);
[matchA,matchB] = propagationStep(imageA,imageB,g,g,centerA,centerB,SuperPatch,SuperPatch,idxA,idxB,matchA,matchB,1);
[matchA,matchB] = propagationStep(imageA,imageB,g,g,centerA,centerB,SuperPatch,SuperPatch,idxA,idxB,matchA,matchB,1);
