function dist = distanceSuperPatchL2(superPatchA, superPatchB, imageA, imagB, idxA, idxB)
end


function dist = distanceSuperPixelL2(superPatchA, superPatchB, imageA, imagB, idxA, idxB)
%superPatchA is the list of the label of the superpixels in the superpatchA
    dist = 0;
    for color = 0:2
        histA = getHist(superPixelA, imageA, idxA,  color);  
        histB = getHist(superPixelB, imageB, idxB,  color) ;
        diff2 = (histA - histB).^2;
        dist = dist + mean(diff2);
    end
end

 function histo = getHist(superPixel, im, idx,  n)
    [numRows,numCols, dim]  = size(im);
    Idx = idx{superPixel};
    Idx = Idx + n*numRows*numCols;
    Values = im(Idx);
    imhisto = hist(Values,[0:255]);  % or imhist(imgray(:),256); with Scilab 
    imhisto =imhisto/sum(imhisto);
    histo = cumsum(imhisto); 
 end
