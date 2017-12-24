function newColor = A_t(idx_pixel,A,B,matchA)
    [nRows,nCols,~] = size(A.im);
    [nRowsB,nColsB,~] = size(B.im);
    p = idxToCoord(idx_pixel,nRows,nCols,A.im);
    label_Ai = A.L(idx_pixel);
    Ai = idxToCoord(A.idx{label_Ai},nRows,nCols,A.im);
    Qi = inv_covQ(label_Ai,10,0.1,Ai);
    [~,sp]=size(A.SuperPatchs{label_Ai});
    w_ = zeros([1,sp]);
    colors = zeros([3,sp]);
    for i = 1:sp
       label_Aj =A.SuperPatchs{label_Ai}(i);
       Aj = idxToCoord(A.idx{label_Aj},nRows,nCols,A.im);
       a = p-a_bar(label_Aj,Aj);
       w_(i) = -a*Qi*transpose(a);
       label_Bj = matchA(label_Aj);
       redIdx = B.idx{label_Bj};
       greenIdx = B.idx{label_Bj}+nRowsB*nColsB;
       blueIdx = B.idx{label_Bj}+2*nRowsB*nColsB;
       colors(:,i) = [mean(B.im(redIdx)),mean(B.im(greenIdx)),mean(B.im(blueIdx))];
       
    end
    sigma_p = max(w_);
    w_ = exp(w_ - sigma_p);
    sum_color = sum(w_.*colors,2)/sum(w_);
%     for label_Aj =A.SuperPatchs{label_Ai}
%         w = exp(-w_
%         sum_w = sum_w + w;
%         label_Bj = matchA(label_Aj);
%         redIdx = B.idx{label_Bj};
%         greenIdx = B.idx{label_Bj}+nRowsB*nColsB;
%         blueIdx = B.idx{label_Bj}+2*nRowsB*nColsB;
%         sum_color = sum_color + w*[mean(B.im(redIdx)),mean(B.im(greenIdx)),mean(B.im(blueIdx))];
%     end
    newColor = sum_color;
    
end

function w = weight(p,Aj,label_Aj,Qi)
    a = p-a_bar(label_Aj,Aj);
    w = exp(-a*Qi*transpose(a));
end

function Qi = inv_covQ(i,delta_s,delta_c,A_i)
    global Q;
    if isempty(Q{i})
       covX_i = cov(A_i(:,1:2));
       covC_i = cov(A_i(:,3:5));
       Qi = blkdiag(delta_s^2*covX_i,delta_c^2*covC_i);
       Qi = inv(Qi);
       Q{i} = Qi;
       
    else
        Qi = Q{i};
    end
end

function Abar = a_bar(i,A_i)
global A_bar;
if isempty(A_bar{i})
   Abar = mean(A_i);
   A_bar{i} = Abar;
else
   Abar = A_bar{i};
end
end

function coord=idxToCoord(id,nRows,nCols,im)%x abscisses , y ordonnee
    x = ceil(id/nRows);
    y = mod(id-1,nRows)+1;
    red = double(im(id));
    green = double(im(id+nRows*nCols));
    blue = double(im(id+2*nRows*nCols));
    coord=[x/nCols,y/nRows,red/255,green/255,blue/255];
end