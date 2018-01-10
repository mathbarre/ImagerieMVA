function newColor = A_t(idx_pixel,A,B,matchA,Q,A_bar)
    [nRows,nCols,~] = size(A.im);
    [nRowsB,nColsB,~] = size(B.im);
    p = idxToCoord(idx_pixel,nRows,nCols,A.im);
    label_Ai = A.L(idx_pixel);
    %Ai = idxToCoord(A.idx{label_Ai},nRows,nCols,A.im);
    %Qi = inv_covQ(label_Ai,10,0.1,Ai);
    %global A_bar;
    %global Q;
    Qi = Q{label_Ai};
    [~,sp]=size(A.SuperPatchs{label_Ai});
    w_ = zeros([1,sp]);
    colors = zeros([3,sp]);
    for j = 1:sp
       label_Aj =A.SuperPatchs{label_Ai}(j);
       %Aj = idxToCoord(A.idx{label_Aj},nRows,nCols,A.im);
       %a = p-a_bar(label_Aj,Aj);
       a = p-A_bar{label_Aj};
       w_(j) = -a*Qi*transpose(a);
       label_Bj = matchA(label_Aj);
       redIdx = B.idx{label_Bj};
       greenIdx = B.idx{label_Bj}+nRowsB*nColsB;
       blueIdx = B.idx{label_Bj}+2*nRowsB*nColsB;
       colors(:,j) = [mean(B.im(redIdx)),mean(B.im(greenIdx)),mean(B.im(blueIdx))];
       
    end
    sigma_p = max(w_);
    w_ = exp(w_ - sigma_p);
    sum_color = sum(w_.*colors,2)/sum(w_);

    newColor = sum_color;
    
end



% function Qi = inv_covQ(i,delta_s,delta_c,A_i)
%     global Q;
%     if isempty(Q{i})
%        covX_i = cov(A_i(:,1:2));
%        covC_i = cov(A_i(:,3:5));
%        Qi = blkdiag(delta_s^2*covX_i,delta_c^2*covC_i);
%        Qi = inv(Qi);
%        Q{i} = Qi;
%        
%     else
%         Qi = Q{i};
%     end
% end


% function Abar = a_bar(i,A_i)
% global A_bar;
% if isempty(A_bar{i})
%    Abar = mean(A_i);
%    A_bar{i} = Abar;
% else
%    Abar = A_bar{i};
% end
% end

