function [newMatchA,newMatchB]=propagationStep(A,B,epsilon,matchA,matchB)

[~,nbSuperPatchsA] = size(A.SuperPatchs);


newMatchA=matchA;
newMatchB=matchB;

for i = 1:nbSuperPatchsA
   neighboorsA = neighboors(i,A.graph);
   for nei = neighboorsA
       if nei < i
          [newMatchForA,oldMatchToSwitch]= permutMatch(A,B,i,nei,epsilon,newMatchA,newMatchB);
          if newMatchForA ~= -1
              if oldMatchToSwitch == -1
                  matchedSuperPixelAInB = newMatchA(i);
                  newMatchA(i) = newMatchForA;
                  newMatchB{newMatchForA} = [newMatchB{newMatchForA},i];
                  newMatchB{matchedSuperPixelAInB} = newMatchB{matchedSuperPixelAInB}((newMatchB{matchedSuperPixelAInB}~=i));
              else
                  matchSuperPixelAInB = newMatchA(i);
                  newMatchA(oldMatchToSwitch) = matchSuperPixelAInB;
                  newMatchA(i) = newMatchForA;
                  newMatchB{newMatchForA} = [newMatchB{newMatchForA},i];
                  newMatchB{newMatchForA} = newMatchB{newMatchForA}((newMatchB{newMatchForA}~=oldMatchToSwitch));
                  newMatchB{matchSuperPixelAInB} = [newMatchB{matchSuperPixelAInB},oldMatchToSwitch];
                  newMatchB{matchSuperPixelAInB} = newMatchB{matchSuperPixelAInB}((newMatchB{matchSuperPixelAInB}~=i));
              end
          end
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
function [newMatchForA,oldMatchToSwitch]=permutMatch(A,B,superpixelA,neighboorA,epsilon,matchA,matchB)
 
    newMatchForA = -1;
    oldMatchToSwitch = -1;
    
    theta = angleBetweenCentre(A.centre,superpixelA,neighboorA); % angle i_i'
    matchSuperPixelAInB = matchA(superpixelA); %le label du superPatch = label du superPixel matché à superPixelA 
    neighboorsMatchSuperPixelAInB = neighboors(matchSuperPixelAInB,B.graph); % les voisins du superPixel precedent
    
    goodOrientationB = minimumAngle(matchSuperPixelAInB,theta,B.centre,neighboorsMatchSuperPixelAInB); % le superPixel de B qui correspond à la meilleur orientation
    newDist = CheckDistance(superpixelA,goodOrientationB,A,B);% nouvelle distance entre patchs centrés en superpixelA et goodOrientationB
    DistA = CheckDistance(superpixelA,matchSuperPixelAInB,A,B);
    
    if(newDist < DistA )
       [~,eps] =  size(matchB{goodOrientationB});
       if (eps < epsilon)%on met goodOrientation comme match pour superPixelsA et on enlève superpixelA de la liste des patch matché avec matchSuperPixelAInB
            newMatchForA = goodOrientationB;
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
            newMatchForA = goodOrientationB;
            oldMatchToSwitch = argmin;
           end
       end
    end
end

function Cost = SwitchingCost(matchA,superpixel1,superpixel2,dist1_2,dist1_1,A,B)
    Cost = dist1_2 -dist1_1...
        + CheckDistance(superpixel1,matchA(superpixel2),A,B) ...
        - CheckDistance(superpixel2,matchA(superpixel2),A,B);
end

function Distance=CheckDistance(superpixelA,superpixelB,A,B)
    global distanceMatrix;
    if distanceMatrix(superpixelA,superpixelB) == -1 
        Distance = distanceSuperPatchL2(superpixelA,superpixelB,A,B);
        distanceMatrix(superpixelA,superpixelB) = Distance;
    else 
        Distance = distanceMatrix(superpixelA,superpixelB);
    end
end