function SVP_recon = SVP(M_masked,mask,step,k,max_iter,tol)
    %--------------------------------------------------------------------------
    % use optspace way to reconstruct the iamge
    %--------------------------------------------------------------------------
    %     main part of optspace
    % 
    %     Inputs:
    %         M_masked             --- original image with masked
    %         mask                 --- index matrix of known elements
    %         r                    --- low rank value
    %         tau                  --- update-step
    %         max_iter             --- max iteration
    %         tol                  --- tolerance of convergence criterion
    %     Outputs: 
    %         M_approx             --- return M_recon
    %--------------------------------------------------------------------------
    
    % Initialization
    [nx,ny] = size(M_masked);
    X = zeros(nx,ny); 
    M_fro_inv=1/norm(M_masked,'fro');
    for iter = 1:max_iter  
        % Update Y
        Y = X-step*(mask.*X-M_masked);
        % Singular value decomposition 
        [U,S0,V] = svds(Y,k); 
        % Update X
        Xtemp = X; 
        X = U*S0*V';
        % Stopping criteria
        epsilon_X = norm(X-Xtemp,'fro') * M_fro_inv;
        if epsilon_X < tol
           break;
        end 
        
    end
 
    SVP_recon = X;
end
    