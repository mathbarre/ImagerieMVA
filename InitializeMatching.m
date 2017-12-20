function [initA,initB] = InitializeMatching(nbSuperPixelsA,nbSuperPixelsB)
    initA = randsample(1:nbSuperPixelsB,nbSuperPixelsA);
    initB = {};
    for i = 1:nbSuperPixelsB
       initB{i} = find(initA==i);
    end
end