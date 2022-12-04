function SVP_recon = SVP(mask_image,mask,step,k,maxIter,tol)
    %
    % This code implements the SVP algorithm
    %
    % More information about SVP can be found in the paper:
    %    Meka, Raghu and Jain, Prateek and Dhillon, Inderjit S, 
    %    "Guaranteed rank minimization via singular value projection", 
    %    arXiv preprint arXiv:0909.5457, 2009.
    %
    %
    % Inputs:
    %    mask_image:  sampled image
    %    mask:  sampled set
    %    step:  step size
    %    k:     the maximum allowable rank
    %    maxIter:  the maximum allowable iterations
    %    tol:   tolerance of convergence criterion
    %
    % Outputs:
    %    SVP_recon: recovered image, obtained by SVP
    %
    % Author: Hao Liang 
    % Last modified by: 21/09/13
    %
    
    % Initialization
    [nx,ny] = size(mask_image);
    X = zeros(nx,ny); 
    M_fro_inv=1/norm(mask_image,'fro');
    for iter = 1:maxIter  
        
        % Update Y
        Y = X-step*(mask.*X-mask_image);
        % Singular value decomposition 
        [U,S0,V] = svds(Y,k); 
        % Update X
        Xtemp = X; 
        X = U*S0*V';
        % Stopping criteria
        TOLL = norm(X-Xtemp,'fro') * M_fro_inv;
        if TOLL < tol
           break;
        end 
        
    end
 
    SVP_recon = X;
end
    