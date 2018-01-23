function [initA,initB] = InitializeMatching(nbSuperPixelsA,nbSuperPixelsB)
    %associate each superpixels in A with a different superpixel in B
    initA = randsample(1:nbSuperPixelsB,nbSuperPixelsA);
    initB = {};
    for i = 1:nbSuperPixelsB
       initB{i} = find(initA==i);
    end
end