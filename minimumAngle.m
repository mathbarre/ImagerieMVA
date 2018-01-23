function argmin=minimumAngle(i,theta,centre,Neighboors)
%given theta_i_i' one looks among neighboors for the superpixel k
%which minimise abs(theta_i_i'+pi-theta_i'_k)
%slight changes from the formula in the paper
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



