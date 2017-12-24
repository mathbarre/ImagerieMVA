function [newMatchA,newMatchB] = randomSearchStep(A,B,scale,epsilon,matchA,matchB,nbOfRandomCandidate)
%radius is the search radius
[~,nbSuperPatchsA] = size(A.SuperPatchs);


newMatchA=matchA;
newMatchB=matchB;

for i = 1:nbSuperPatchsA
   [newMatchForA,oldMatchToSwitch]=search(i,A,B,scale,epsilon,newMatchA,newMatchB,nbOfRandomCandidate);
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

function [newMatchForA,oldMatchToSwitch]=search(superpixelA,A,B,scale,epsilon,matchA,matchB,nbOfRandomCandidate)
newMatchForA = -1;
oldMatchToSwitch = -1;
matchSuperPixelAInB = matchA(superpixelA);
[n,m,~] = size(B.im);
c = B.centre(matchSuperPixelAInB).Centroid;
%w = [min([(m-c(1)),(c(1)-1)]),min([(n-c(2)),(c(2)-1)])];
DistA = CheckDistance(superpixelA,matchSuperPixelAInB,A,B);
minDist = DistA;
finalCandidate = -1;
for i = 1:nbOfRandomCandidate
    u = 2*rand([1,2])-1;
    candidate =floor(c+scale*[m,n].*u);
    candidate = [mod(candidate(1),m)+1,mod(candidate(2),n)+1];
    superpixelCandidate = B.L(candidate(2),candidate(1));
    newDist = CheckDistance(superpixelA,superpixelCandidate,A,B);
    if(newDist < minDist )
      minDist = newDist;
      finalCandidate = superpixelCandidate;
    end
end

if(finalCandidate ~= -1)
    [~,eps] =  size(matchB{finalCandidate});
    if (eps < epsilon)
        newMatchForA = finalCandidate;
    else
       cost = 100000;
       argmin = 0;
       for SwitchCandidate = matchB{finalCandidate}
           switchCost = SwitchingCost(matchA,superpixelA,SwitchCandidate,minDist,DistA,A,B);
           if(switchCost < cost)
              cost = switchCost;
              argmin = SwitchCandidate;
           end
       end
       if(cost < 0)
        newMatchForA = finalCandidate;
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

