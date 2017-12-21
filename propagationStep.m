function [newMatchA,newMatchB]=propagationStep(imA,imB,graphA,graphB,centreA,centreB,SuperPatchsA,SuperPatchsB,idxA,idxB,matchA,matchB,epsilon)

[~,nbSuperPatchsA] = size(SuperPatchsA);


newMatchA=matchA;
newMatchB=matchB;

for i = 1:nbSuperPatchsA
   neighboorsA = neighboors(i,graphA);
   for nei = neighboorsA
       if nei < i
          [newMatchA,newMatchB]= permutMatch(imA,imB,graphB,idxA,idxB,i,nei,centreA,centreB,epsilon,newMatchA,newMatchB,SuperPatchsA,SuperPatchsB);
       end
   end
end

%for j = fliplr(1:nbSuperPatchsA)
%   neighboorsA = neighboors(j,graphA);
%   for nei = neighboorsA
%       if nei > j
%          [newMatchA,newMatchB]= permutMatch(imA,imB,graphB,idxA,idxB,j,nei,centreA,centreB,epsilon,newMatchA,newMatchB,SuperPatchsA,SuperPatchsB);
%       end
%   end
%end

end

%superPixelA = label du superPixel, neighboorA = label d'un superPixels
%voisin de A, matchA = array Superpixel de A --> superpixel de B,
%matchB = cell array , superpixel de B --> list superpixel de A
function [newMatchA,newMatchB]=permutMatch(imA,imB,graphB,idxA,idxB,superpixelA,neighboorA,centreA,centreB,epsilon,matchA,matchB,SuperPatchsA,SuperPatchsB)
    newMatchA = matchA;
    newMatchB = matchB;
    theta = angleBetweenCentre(centreA,superpixelA,neighboorA); % angle i_i'
    matchSuperPixelAInB = matchA(superpixelA); %le label du superPatch = label du superPixel matché à superPixelA 
    neighboorsMatchSuperPixelAInB = neighboors(matchSuperPixelAInB,graphB); % les voisins du superPixel precedent
    goodOrientationB = minimumAngle(matchSuperPixelAInB,theta,centreB,neighboorsMatchSuperPixelAInB); % le superPixel de B qui correspond à la meilleur orientation
    newDist = distanceSuperPatchL2(SuperPatchsA{superpixelA},SuperPatchsB{goodOrientationB},superpixelA,goodOrientationB, centreA, centreB, imA, imB, idxA, idxB);% nouvelle distance entre patchs centrés en superpixelA et goodOrientationB
    DistA = distanceSuperPatchL2(SuperPatchsA{superpixelA},SuperPatchsA{matchSuperPixelAInB},superpixelA,matchSuperPixelAInB, centreA, centreB, imA, imB, idxA, idxB);
    if(newDist < DistA )
       [~,eps] =  size(matchB{goodOrientationB});
       if (eps < epsilon)%on met goodOrientation comme match pour superPixelsA et on enlève superpixelA de la liste des patch matché avec matchSuperPixelAInB
           newMatchA(superpixelA) = goodOrientationB;
           newMatchB{goodOrientationB} = [matchB{goodOrientationB},superpixelA];
           newMatchB{matchSuperPixelAInB} = matchB{matchSuperPixelAInB}((matchB{matchSuperPixelAInB}~=superpixelA));
       else
           cost = 100000;
           argmin = 0;
           for SwitchCandidate = matchB{goodOrientationB}
               switchCost = SwitchingCost(SuperPatchsA,SuperPatchsB,matchA,superpixelA,SwitchCandidate,newDist,DistA,centreA, centreB, imA, imB, idxA, idxB);
               if(switchCost < cost)
                  cost = switchCost;
                  argmin = SwitchCandidate;
               end
           end
           if(cost < 0)
              newMatchA(argmin) = matchA(superpixelA);
              newMatchA(superpixelA) = goodOrientationB;
              newMatchB{goodOrientationB} = [matchB{goodOrientationB},superpixelA];
              newMatchB{goodOrientationB} = matchB{goodOrientationB}((matchB{goodOrientationB}~=argmin));
              newMatchB{matchSuperPixelAInB} = [matchB{matchSuperPixelAInB},argmin];
              newMatchB{matchSuperPixelAInB} = matchB{matchSuperPixelAInB}((matchB{matchSuperPixelAInB}~=superpixelA));
           end
       end
    end
end

function Cost = SwitchingCost(SuperPatchA,SuperPatchB,matchA,superpixel1,superpixel2,dist1_2,dist1_1, centreA, centreB, imA, imB, idxA, idxB)
    Cost = dist1_2 -dist1_1...
        + distanceSuperPatchL2(SuperPatchA{superpixel2},SuperPatchB{matchA(superpixel1)},superpixel1,matchA(superpixel2), centreA, centreB, imA, imB, idxA, idxB) ...
        - distanceSuperPatchL2(SuperPatchA{superpixel2},SuperPatchB{matchA(superpixel2)},superpixel2,matchA(superpixel2), centreA, centreB, imA, imB, idxA, idxB);
end