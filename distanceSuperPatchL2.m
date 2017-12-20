function dist = distanceSuperPatchL2(superPatchA, superPatchB, image, idx)
%superPatchA is the list of the label of the superpixels in the superpatchA
    dist = 0;
    for color = 0:2
        color
        histA = getHist(superPatchA, image, idx,  color);  
        histB = getHist(superPatchB, image, idx,  color) ;
        diff2 = (histA - histB).^2;
        dist = dist + mean(diff2);
    end
end

 function histo = getHist(superPatchA, im, idx,  n)
    n
    [numRows,numCols, plop]  = size(im);
    numRows
    numCols
    
    redIdxA = idx{superPatchA};
    redIdxA + n*numRows*numCols;
    redValues = im(redIdxA);
    imhistoA = hist(redValues,[0:255]);  % or imhist(imgray(:),256); with Scilab 
    imhistoA =imhistoA/sum(imhistoA);
    histo = cumsum(imhistoA); 
 end
 