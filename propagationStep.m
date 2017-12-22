function [newMatchA,newMatchB]=propagationStep(A,B,epsilon,matchA,matchB)

[~,nbSuperPatchsA] = size(A.SuperPatchs);


newMatchA=matchA;
newMatchB=matchB;

for i = 1:nbSuperPatchsA
   neighboorsA = neighboors(i,A.graph);
   for nei = neighboorsA
       if nei < i
          [newMatchA,newMatchB]= permutMatch(A,B,i,nei,epsilon,newMatchA,newMatchB);
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
function [newMatchA,newMatchB]=permutMatch(A,B,superpixelA,neighboorA,epsilon,matchA,matchB)
    newMatchA = matchA;
    newMatchB = matchB;
    theta = angleBetweenCentre(A.centre,superpixelA,neighboorA); % angle i_i'
    matchSuperPixelAInB = matchA(superpixelA); %le label du superPatch = label du superPixel matché à superPixelA 
    neighboorsMatchSuperPixelAInB = neighboors(matchSuperPixelAInB,B.graph); % les voisins du superPixel precedent
    goodOrientationB = minimumAngle(matchSuperPixelAInB,theta,B.centre,neighboorsMatchSuperPixelAInB); % le superPixel de B qui correspond à la meilleur orientation
    newDist = distanceSuperPatchL2(superpixelA,goodOrientationB,A,B);% nouvelle distance entre patchs centrés en superpixelA et goodOrientationB
    DistA = distanceSuperPatchL2(superpixelA,matchSuperPixelAInB,A,B);
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
               switchCost = SwitchingCost(matchA,superpixelA,SwitchCandidate,newDist,DistA,A,B);
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

function Cost = SwitchingCost(matchA,superpixel1,superpixel2,dist1_2,dist1_1,A,B)
    Cost = dist1_2 -dist1_1...
        + distanceSuperPatchL2(superpixel1,matchA(superpixel2),A,B) ...
        - distanceSuperPatchL2(superpixel2,matchA(superpixel2),A,B);
end