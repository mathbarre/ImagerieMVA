function [newMatchA,newMatchB]=propagationStep(imA,imB,graphA,graphB,centreA,centreB,SuperPatchsA,SuperPatchsB,idxA,idxB,matchA,matchB,epsilon)

[~,nbSuperPatchsA] = size(SuperPatchsA);
[~,nbSuperPatchsB] = size(SuperPatchsB);

newMatchA=matchA;
newMatchB=matchB;

for i = 1:nbSuperPatchsA
   neighboorsA = neighboors(i,graphA);
   for nei = neighboorsA
       if nei < i
          [newMatchA,newMatchB]= permutMatch(imA,imB,graphB,idxA,idxB,i,nei,centreA,centreB,epsilon,newMatchA,newMatchB,superPatchsA,SuperPatchsB);
       end
   end
end

for j = 1:nbSuperPatchsA
   neighboorsA = neighboors(i,graphA);
   for nei = neighboorsA
       if nei < i
          [newMatchA,newMatchB]= permutMatch(imA,imB,graphB,idxA,idxB,i,nei,centreA,centreB,epsilon,newMatchA,newMatchB,superPatchsA,SuperPatchsB);
       end
   end
end

end

%superPixelA = label du superPixel, neighboorA = label d'un superPixels
%voisin de A, matchA = array Superpixel de A --> superpixel de B,
%matchB = cell array , superpixel de B --> list superpixel de A
function [newMatchA,newMatchB]=permutMatch(imA,imB,graphB,idxA,idxB,superpixelA,neighboorA,centreA,centreB,epsilon,matchA,matchB,SuperPatchsA,SuperPatchsB)
    newMatchA = matchA;
    newMatchB = matchB;
    theta = angleBetweenCentre(centreA,superpixelA,neighboorA); % angle i_i'
    matchSuperPixelAInB = matchA(superpixelA); %le label du superPatch = label du superPixel match� � superPixelA 
    neighboorsMatchSuperPixelAInB = neighboors(matchSuperPixelAInB,graphB); % les voisins du superPixel precedent
    goodOrientationB = minimumAngle(matchSuperPixelAInB,theta,neighboorsMatchSuperPixelAInB); % le superPixel de B qui correspond � la meilleur orientation
    newDist = distanceSuperPatchL2(SuperPatchsA{superpixelA},SuperPatchsB{goodOrientationB});% nouvelle distance entre patchs centr�s en superpixelA et goodOrientationB
    DistA = distanceSuperPatchL2(SuperPatchsA{superpixelA},SuperPatchsA{matchSuperPixelAInB});
    if(newDist < DistA )
       [~,eps] =  size(matchB{goodOrientationB});
       if (eps < epsilon)%on met goodOrientation comme match pour superPixelsA et on enl�ve superpixelA de la liste des patch match� avec matchSuperPixelAInB
           newMatchA(superPixelA) = goodOrientationB;
           newMatchB{goodOrientationB} = [matchB{goodOrientationB},superPixelA];
           newMatchB{matchSuperPixelAInB} = matchB{matchSuperPixelAInB}((matchB{matchSuperPixelAInB}~=superPixelA));
       else
           cost = 100000;
           argmin = 0;
           for SwitchCandidate = matchB{goodOrientationB}
               switchCost = SwitchingCost(SuperPatchA,SuperPatchB,matchA,superPixelA,SwitchCandidate,newDist,DistA);
               if(switchCost < cost)
                  cost = switchCost;
                  argmin = SwitchCandidate;
               end
           end
           if(cost < 0)
              newMatchA(argmin) = match(superPixelA);
              newMatchA(superPixelA) = goodOrientationB;
              newMatchB{goodOrientationB} = [matchB{goodOrientationB},superPixelA];
              newMatchB{goodOrientationB} = matchB{goodOrientationB}((matchB{goodOrientationB}~=argmin));
              newMatchB{matchSuperPixelAInB} = [matchB{matchSuperPixelAInB},argmin];
              newMatchB{matchSuperPixelAInB} = matchB{matchSuperPixelAInB}((matchB{matchSuperPixelAInB}~=superPixelA));
           end
       end
    end
end

function Cost = SwitchingCost(SuperPatchA,SuperPatchB,matchA,superpixel1,superpixel2,dist1_2,dist1_1)
    Cost = dist1_2 -dist1_1...
        + distanceSuperPatchL2(SuperPatchA{superpixel2},SuperPatchB(matchA(superpixel1))) ...
        - distanceSuperPatchL2(SuperPatchA{superpixel2},SuperPatchB{matchA(superpixel2)});
end