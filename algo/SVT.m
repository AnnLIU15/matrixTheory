function SVT_recon = SVT(mask_image,mask,tau,step,maxIter,tol)
%
% This code implements the SVT algorithm
%
% More information about SVT can be found in the paper:
%    Cai, Jian-Feng and Candès, Emmanuel J and Shen, Zuowei, 
%    "A singular value thresholding algorithm for matrix completion", 
%    SIAM Journal on optimization, 2010.
%
%
% Inputs:
%    mask_image:  sampled image
%    mask:  sampled set
%    tau:   threshold parameter
%    step:  step size
%    increment: 5 
%    maxIter:  maximum allowable iterations
%    tol:   tolerance of convergence criterion
%
% Outputs:
%    SVT_recon:  recovered image, obtained by SVT
%
% Author: Hao Liang 
% Last modified by: 21/09/13
%

% Initialization
X = mask_image;
Y = zeros('like', mask_image);
M_mask_inv = 1 / norm(mask_image,'fro');
% r_k = 0;
% l = 5;
for i = 1:maxIter
    % s_k = r_k + 1;
    % Singular value decomposition to update X
    [U,S,V] = svd(Y); 

    % XTemp = X; 
    X = U * wthresh(S,'s',tau) *V';
    
    % Update Y
    Y = Y+step*(mask_image-X);
    Y = mask.*Y;
    
    % Stopping criteria
    TOLL = norm(mask.*X - mask_image,'fro')*M_mask_inv;
    if TOLL<tol
        break;
    end
    
end

SVT_recon = X;

end