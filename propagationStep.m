function [newMatchA,newMatchB]=propagationStep(A,B,epsilon,matchA,matchB)

[~,nbSuperPatchsA] = size(A.SuperPatchs);


newMatchA=matchA;
newMatchB=matchB;

for i = 1:nbSuperPatchsA
   neighboorsA = A.graph{i};
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
for j = fliplr(1:nbSuperPatchsA)
   neighboorsA = A.graph{j};
   for nei = neighboorsA
       if nei > j
          [newMatchForA,oldMatchToSwitch]= permutMatch(A,B,j,nei,epsilon,newMatchA,newMatchB);
          if newMatchForA ~= -1
              if oldMatchToSwitch == -1
                  matchedSuperPixelAInB = newMatchA(j);
                  newMatchA(j) = newMatchForA;
                  newMatchB{newMatchForA} = [newMatchB{newMatchForA},j];
                  newMatchB{matchedSuperPixelAInB} = newMatchB{matchedSuperPixelAInB}((newMatchB{matchedSuperPixelAInB}~=j));
              else
                  matchSuperPixelAInB = newMatchA(j);
                  newMatchA(oldMatchToSwitch) = matchSuperPixelAInB;
                  newMatchA(j) = newMatchForA;
                  newMatchB{newMatchForA} = [newMatchB{newMatchForA},j];
                  newMatchB{newMatchForA} = newMatchB{newMatchForA}((newMatchB{newMatchForA}~=oldMatchToSwitch));
                  newMatchB{matchSuperPixelAInB} = [newMatchB{matchSuperPixelAInB},oldMatchToSwitch];
                  newMatchB{matchSuperPixelAInB} = newMatchB{matchSuperPixelAInB}((newMatchB{matchSuperPixelAInB}~=j));
              end
          end
       end
   end
end


end

%superPixelA = label du superPixel, neighboorA = label d'un superPixels
%voisin de A, matchA = array Superpixel de A --> superpixel de B,
%matchB = cell array , superpixel de B --> list superpixel de A
function [newMatchForA,oldMatchToSwitch]=permutMatch(A,B,superpixelA,neighboorA,epsilon,matchA,matchB)
 
    newMatchForA = -1;
    oldMatchToSwitch = -1;
    
    theta = angleBetweenCentre(A.centre,superpixelA,neighboorA); % angle i_i'
    matchSuperPixelAInB = matchA(superpixelA); %le label du superPatch = label du superPixel matché à superPixelA 
    neighboorsMatchSuperPixelAInB = B.graph{matchSuperPixelAInB}; % les voisins du superPixel precedent
    
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
        + CheckDistance(superpixel2,matchA(superpixel1),A,B) ...
        - CheckDistance(superpixel2,matchA(superpixel2),A,B);
end

