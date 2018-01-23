function newColor = A_t(idx_pixel,A,B,matchA,Q,A_bar)
    %Compute the new color A_t(p)
    [nRows,nCols,~] = size(A.im);
    [nRowsB,nColsB,~] = size(B.im);
    
    p = idxToCoord(idx_pixel,nRows,nCols,A.im);
    label_Ai = A.L(idx_pixel);
    
    Qinv_i = Q{label_Ai};
    
    %For the reconstruction, we use only superpixels at a certain distance
    %from the current superpixel
    [~,sp]=size(A.SuperPatchs{label_Ai});
    w_ = zeros([1,sp]);
    colors = zeros([3,sp]);
    for j = 1:sp
       label_Aj =A.SuperPatchs{label_Ai}(j);
       a = p-A_bar{label_Aj};
       w_(j) = -a*Qinv_i*transpose(a);
       label_Bj = matchA(label_Aj);
       colors(:,j) = B.mean(label_Bj,:);
    end
    sigma_p = max(w_);
    w_new = w_-sigma_p;
    exp_w_ = exp(w_new);
    sum_color = sum(exp_w_.*colors,2)/sum(exp_w_);

    newColor = sum_color;
    
end



