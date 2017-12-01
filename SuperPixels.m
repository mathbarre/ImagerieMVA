im = (imread('TP/im/simpson512.png')); 
[L,N] = superpixels(im,100);
im(28,34,:) = 255;
im(55,79,:) = 255;
im(94,99,:) = 255;
figure;
BW = boundarymask(L);
imshow(imoverlay(im,BW,'cyan'),[])


outputImage = zeros(size(im),'like',im);
idx = label2idx(L);
numRows = size(im,1);
numCols = size(im,2);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = mean(im(redIdx));
    outputImage(greenIdx) = mean(im(greenIdx));
    outputImage(blueIdx) = mean(im(blueIdx));
end    

figure;

imshow(outputImage)

centr = regionprops(L,'centroid');
g = adjacentRegionsGraph(L);

nei = neighboors(13,g);
theta = angleBetweenCentre(centr,13,2);

argmin =minimumAngle(13,2,centr,nei);

function nei=neighboors(i,graph)
a = graph.Edges.EndNodes(graph.Edges.EndNodes(:,2)==i,1);
b = graph.Edges.EndNodes(graph.Edges.EndNodes(:,1)==i,2);
nei = transpose(cat(1,a,b));
end

function theta =angle2d(u,v)
theta =atan2(u(1)*v(2)-u(2)*v(1),dot(u,v));
end

%calcul le theta_i_i' du papier SuperPatch Match
function theta = angleBetweenCentre(centre,i,iprime)
p = centre(iprime).Centroid - centre(i).Centroid;
theta = angle2d([0,1],[p(2),p(1)]);
if theta < 0
    theta = 2*pi+theta;
end
end

%etant donné theta_i_i' on cherche parmis neighboors le super pixel k
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