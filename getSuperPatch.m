%very naive
%renvoie cell array, superPatch{i} donne la liste des label des superpixels
%contenu dans le superpatch de centre le superpixel-i
function superPatch=getSuperPatch(centre,r)
    [n,d]= size(centre);
    superPatch = {};
    for i = 1:n 
        l = [];
        for j = 1:n
           if(norm(centre(i).Centroid - centre(j).Centroid) <= r)
               l = [l j];
           end
        end
        superPatch{i}=l;
    end
end