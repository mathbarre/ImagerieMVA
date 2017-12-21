function dist = distanceSuperPatchL2(superPatchA, superPatchB, labelSuperPixelCentralA,labelSuperPixelCentralB, centerA, centerB, imageA, imageB, idxA, idxB)
    sigma21 = 0.01
    sigma22 = 0.01
    ci = centerA(labelSuperPixelCentralA).Centroid;
    cj = centerB(labelSuperPixelCentralB).Centroid;
    [~, n] = size(superPatchA);
    [~, m] = size(superPatchB);
    weigthMatrix = zeros(n,m);
    weightTimeDistBetweenPixelsMatrix = zeros(n,m);
    
    for iPrime = 1:n
        labelSuperPixelA = superPatchA(iPrime);
        for jPrime = 1:m
            labelSuperPixelB = superPatchB(jPrime);
            ciPrime = centerA(labelSuperPixelA ).Centroid;
            cjPrime = centerB(labelSuperPixelB).Centroid;
            weigthMatrix(iPrime,jPrime) = weightSuperPixels(ci, cj, ciPrime, cjPrime, sigma21, sigma22);
            distSuperPixel = distanceSuperPixelL2(labelSuperPixelA, labelSuperPixelB, imageA, imageB, idxA, idxB);
            weightTimeDistBetweenPixelsMatrix(iPrime,jPrime) =  weigthMatrix(iPrime,jPrime) *distSuperPixel ;
        end
    end
    
    dist = sum(sum(weightTimeDistBetweenPixelsMatrix));
    dist = dist / sum(sum(weigthMatrix));

end

function w = weightSuperPixels(ci, cj, ciPrime, cjPrime, sigma21, sigma22)
    xiPrimejPrime = (cjPrime-cj) - (ciPrime - ci);
    w = exp(- sum( xiPrimejPrime.^2)/sigma21);
    w = w  * weigthSuperPixel(ci, ciPrime, sigma22) * weigthSuperPixel(cj, cjPrime, sigma22);
end
    
     
 function w = weigthSuperPixel(ci, ciPrime, sigma22)
    w = exp(-sum((ci-ciPrime).^2) / sigma22);
 end
    


function dist = distanceSuperPixelL2(labelSuperPixelA, labelSuperPixelB, imageA, imageB, idxA, idxB)
%superPatchA is the list of the label of the superpixels in the superpatchA
    dist = 0;
    for color = 0:2;
        histA = getHist(labelSuperPixelA, imageA, idxA,  color);  
        histB = getHist(labelSuperPixelB, imageB, idxB,  color) ;
        diff2 = (histA - histB).^2;
        dist = dist + mean(diff2);
    end
end

 function histo = getHist(labelSuperPixel, im, idx,  n)
    [numRows,numCols, ~]  = size(im);
    Idx = idx{labelSuperPixel};
    Idx = Idx + n*numRows*numCols;
    Values = im(Idx);
    imhisto = hist(Values,[0:255]);  % or imhist(imgray(:),256); with Scilab 
    imhisto =imhisto/sum(imhisto);
    histo = cumsum(imhisto); 
 end

 
    
    