function resultIm=SCT(imA,imB,epsilon,lambda,R,alpha,nbIter,nbSuperpixelWantedA,nbSuperpixelWantedB,delta_s,delta_c,colorspace)
    prompt =sprintf('Preprocessing...\n');
    fprintf(prompt)
    A = create_object(imA,nbSuperpixelWantedA,'A',R,colorspace);
    B = create_object(imB,nbSuperpixelWantedB,'B',R,colorspace);
    [numRowsA,numColsA,~] = size(imA);
    global distanceMatrix;
    distanceMatrix = zeros([A.N,B.N])-1;
    global distanceSuperPixelMatrix;
    distanceSuperPixelMatrix = zeros([A.N,B.N])-1;
    [matchA,matchB] = InitializeMatching(A.N,B.N);
    prompt = sprintf('Initialisation du matching terminée\n');
    fprintf(prompt)
    
    for i = 1:nbIter
        [matchA,matchB] = propagationStep(A,B,epsilon,matchA,matchB);
        %meanDistance(A,B,matchA)
        [matchA,matchB] = randomSearchStep(A,B,alpha,epsilon,matchA,matchB,1);
        prompt =sprintf(strcat('Étape ',num2str(i)," du matching : Distance totale d'appareillage = ",num2str(meanDistance(A,B,matchA))," dans ",colorspace,'\n'));
        fprintf(prompt)

        %plot(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA)
        %save(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA,i)
    end
    
    
    %compute the matrix Q^-1 and the means A_bar in advance
    Q = cell([1,A.N]);

    A_bar = cell([1,A.N]);

    for i=1:A.N
        A_i = idxToCoord(A.idx{i},numRowsA,numColsA,A.im);
        covX_i = cov(A_i(:,1:2));
        covC_i = cov(A_i(:,3:5));
        Qi = blkdiag(delta_s^2*covX_i,delta_c^2*covC_i)+lambda*diag(ones([1,5]));
        Qi = inv(Qi);
        Q{i} = Qi;
        Abar = mean(A_i);
        A_bar{i} = Abar;
    end
    
    TransformedImage = zeros(size(A.im));


    TransformedImage1 = zeros([numRowsA,numColsA]);
    TransformedImage2 = zeros([numRowsA,numColsA]);
    TransformedImage3 = zeros([numRowsA,numColsA]);
    
    prompt = sprintf('Reconstruction de couleurs...\n');
    fprintf(prompt)

    parfor p = 1:(numRowsA*numColsA)
        a = A_t(p,A,B,matchA,Q,A_bar);
        TransformedImage1(p) = a(1);
        TransformedImage2(p) = a(2);
        TransformedImage3(p) = a(3);
    end

    TransformedImage(:,:,1)= TransformedImage1;
    TransformedImage(:,:,2)= TransformedImage2;
    TransformedImage(:,:,3)= TransformedImage3;
    
    resultIm = TransformedImage/255;
    disp("Fin");
end

function S=create_object(im,nbSuperPixelsWanted,which,R,colorspace)
    %preprocessing step which compute all the quantity we're going to need
    
    S=struct;
    [numRows,numCols,~] = size(im);
    [L,N] = superpixels(im, nbSuperPixelsWanted);
    idx = label2idx(L);
    if(which=='B')
       meanB = zeros([N,3]);
       for labelVal = 1:N
            redIdx = idx{labelVal};
            greenIdx = idx{labelVal}+numRows*numCols;
            blueIdx = idx{labelVal}+2*numRows*numCols;
            mr = mean(im(redIdx));
            mg = mean(im(greenIdx));
            mb = mean(im(blueIdx));
            meanB(labelVal,:)=[mr,mg,mb];
       end
       S.mean = meanB;
    end
    centr = regionprops(L,'centroid');
    g = adjacentRegionsGraph(L);
    SuperPatch = getSuperPatch(centr,R);
    S.im = im;
    S.L = L;
    S.graph = neighboors(g);
    S.centre = centr;
    S.idx=idx;
    S.N =N;
    S.SuperPatchs =SuperPatch;
    
    %usefull for colorspace change
    O = [1/sqrt(3), 1/sqrt(3), 1/sqrt(3) ; 1/sqrt(2) -1/sqrt(2) 0; 1/sqrt(6) 1/sqrt(6) -2/sqrt(6)];
    S.opp = reshape(reshape(double(im),numRows*numCols,3)*O,numRows,numCols,3);
    LMS = [0.3811 0.5783 0.0402 ; 0.1967 0.7244 0.0782 ; 0.0241 0.1288 0.8444];
    S.lms = log10(0.1+reshape(reshape(double(im),numRows*numCols,3)*LMS',numRows,numCols,3));
    S.lab = reshape(reshape(S.lms,numRows*numCols,3)*O,numRows,numCols,3);
    S.hsv = rgb2hsv(double(im));
    
    histS = cell(N);

    for i = 1:N 
        histo = zeros([256,3]);
        for c = 1:3
            histo(:,c) = getHist(i,S,c-1,colorspace);
        end
        histS{i} = histo;
    end
    S.hist = histS;
    


    
end

function histo = getHist(labelSuperPixel,A,n,colorspace)
    switch colorspace
        case 'rgb'
            A_im = A.im;
        case 'opp'
            A_im = A.opp;
        case 'lms'
            A_im = A.lms;
        case 'lab'
            A_im = A.lab;
        case 'hsv'
            A_im = A.hsv;
        otherwise
            A_im = A.im;
    end
    [numRows,numCols, ~]  = size(A.im);
    Idx = A.idx{labelSuperPixel};
    Idx = Idx + n*numRows*numCols;
    Values = A_im(Idx);
    imhisto = hist(Values,[0:255]);   
    imhisto =imhisto/sum(imhisto);
    histo = cumsum(imhisto); 
end


function save(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA,i)
    for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    redIdxM = idxB{matchA(labelVal)};
    greenIdxM = idxB{matchA(labelVal)}+numRowsB*numColsB;
    blueIdxM = idxB{matchA(labelVal)}+2*numRowsB*numColsB;
    outputImage(redIdx) = mean(imB(redIdxM));
    outputImage(greenIdx) = mean(imB(greenIdxM));
    outputImage(blueIdx) = mean(imB(blueIdxM));
    end
    imwrite(double(outputImage)/255, strcat('Results/',num2str(i),'step.png'));
end


function plot(N,idx,idxB,numRows,numCols,numRowsB,numColsB,outputImage,imB,matchA)
    %to get plot of averaged color
    for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    redIdxM = idxB{matchA(labelVal)};
    greenIdxM = idxB{matchA(labelVal)}+numRowsB*numColsB;
    blueIdxM = idxB{matchA(labelVal)}+2*numRowsB*numColsB;
    outputImage(redIdx) = mean(imB(redIdxM));
    outputImage(greenIdx) = mean(imB(greenIdxM));
    outputImage(blueIdx) = mean(imB(blueIdxM));
    end
    figure; imshow(outputImage);
end

function meanDist = meanDistance(A,B,matchA)
    %to see the global matching distance
    [~,n] = size(matchA);
    dist = 0;
    for i = 1:n
       dist = dist + CheckDistance(i,matchA(i),A,B); 
    end
    meanDist = dist/n;
end