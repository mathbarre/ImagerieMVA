function superPatch=getSuperPatch(centre,r)
    %very naive
    %return cell array, superPatch{i} list of labels of superpixels
    [n,~]= size(centre);
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