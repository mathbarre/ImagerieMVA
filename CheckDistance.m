function Distance=CheckDistance(superpixelA,superpixelB,A,B)
    %To avoid doing computation twice 
    global distanceMatrix;
    if distanceMatrix(superpixelA,superpixelB) == -1 
        Distance = distanceSuperPatchL2(superpixelA,superpixelB,A,B);
        distanceMatrix(superpixelA,superpixelB) = Distance;
    else 
        Distance = distanceMatrix(superpixelA,superpixelB);
    end
end