function dist = distanceSuperPatchL2(SuperPixelCentralA,SuperPixelCentralB,A,B)
    %Comput the distance between two superpatchs, currently compute
    % histogram L2 distance between central superpixels.
    %To compute the distance between global superpatchs, decomment the
    %following code

%     [h,w,~] = size(A.im);
%     [~,K] = size(A.SuperPatchs);
%     sigma21 = h*w/K/4;
%     global R;
%     sigma22 = 2*R^2;
%     ci = A.centre(SuperPixelCentralA).Centroid;
%     cj = B.centre(SuperPixelCentralB).Centroid;
%     [~, n] = size(A.SuperPatchs{SuperPixelCentralA});
%     [~, m] = size(B.SuperPatchs{SuperPixelCentralB});
%     weigthMatrix = zeros(n,m);
%     weightTimeDistBetweenPixelsMatrix = zeros(n,m);
%     for iPrime = 1:n
%         labelSuperPixelA = A.SuperPatchs{SuperPixelCentralA}(iPrime);
%         for jPrime = 1:m
%             labelSuperPixelB = B.SuperPatchs{SuperPixelCentralB}(jPrime);
%             ciPrime = A.centre(labelSuperPixelA).Centroid;
%             cjPrime = B.centre(labelSuperPixelB).Centroid;
%             weigthMatrix(iPrime,jPrime) = weightSuperPixels(ci, cj, ciPrime, cjPrime, sigma21, sigma22);
%             distSuperPixel =checkDistanceSuperpixel(labelSuperPixelA, labelSuperPixelB,A,B);
%             weightTimeDistBetweenPixelsMatrix(iPrime,jPrime) =  weigthMatrix(iPrime,jPrime) *distSuperPixel ;
%         end
%     end
%     
%     dist = sum(sum(weightTimeDistBetweenPixelsMatrix));
%     dist = dist / sum(sum(weigthMatrix));
    dist =checkDistanceSuperpixel(SuperPixelCentralA, SuperPixelCentralB,A,B);
end

function w = weightSuperPixels(ci, cj, ciPrime, cjPrime, sigma21, sigma22)
    xiPrimejPrime = (cjPrime-cj) - (ciPrime - ci);
    w = exp(- sum( xiPrimejPrime.^2)/sigma21);
    w = w  * weigthSuperPixel(ci, ciPrime, sigma22) * weigthSuperPixel(cj, cjPrime, sigma22);
end
    
     
 function w = weigthSuperPixel(ci, ciPrime, sigma22)
    w = exp(-sum((ci-ciPrime).^2) / sigma22);
 end
    


function dist = distanceSuperPixelL2(labelSuperPixelA, labelSuperPixelB,A,B)
%superPatchA is the list of the label of the superpixels in the superpatchA
    histA = A.hist{labelSuperPixelA};
    histB = B.hist{labelSuperPixelB};
    dist = sqrt(sum(sum((histA-histB).^2)));
end


 function Distance=checkDistanceSuperpixel(superpixelA, superpixelB,A,B)
    global distanceSuperPixelMatrix;
     if distanceSuperPixelMatrix(superpixelA,superpixelB) == -1 
        Distance = distanceSuperPixelL2(superpixelA,superpixelB,A,B);
        distanceSuperPixelMatrix(superpixelA,superpixelB) = Distance;
    else 
        Distance = distanceSuperPixelMatrix(superpixelA,superpixelB);
    end
 end
    
    