%etant donne theta_i_i' on cherche parmis neighboors le super pixel k
%qui minimise abs(theta_i_i'+pi-theta_i'_k)
function argmin=minimumAngle(i,theta,centre,Neighboors)
min = 4.5;
argmin = i;
for k = Neighboors
    theta_k = angleBetweenCentre(centre,i,k);
    min_k =abs(theta+pi-theta_k);
    if min_k <= min
       min = min_k;
       argmin=k;
    end
end
end


