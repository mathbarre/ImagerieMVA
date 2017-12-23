%etant donne theta_i_i' on cherche parmis neighboors le super pixel k
%qui minimise abs(theta_i_i'+pi-theta_i'_k)
function argmin=minimumAngle(i,theta,centre,Neighboors)
min = 10;
argmin = i;
for k = Neighboors
    theta_k = angleBetweenCentre(centre,i,k);
    min_k =abs(mod(theta+pi,2*pi)-theta_k);
    if min_k <= min
       min = min_k;
       argmin=k;
    end
end
end



