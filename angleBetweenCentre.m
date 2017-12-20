%calcul le theta_i_i' du papier SuperPatch Match
function theta = angleBetweenCentre(centre,i,iprime)
p = centre(iprime).Centroid - centre(i).Centroid;
theta = angle2d([0,1],[p(2),p(1)]);
if theta < 0
    theta = 2*pi+theta;
end
end

function theta =angle2d(u,v)
theta =atan2(u(1)*v(2)-u(2)*v(1),dot(u,v));
end