function SVT_recon = SVT(M_masked,mask,tau,step,max_iter,tol)
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
    X = M_masked;
    Y = zeros(size(M_masked));
    M_mask_inv = 1 / norm(M_masked,'fro');
    % r_k = 0;
    % l = 5;
    for i = 1:max_iter
        % s_k = r_k + 1;
        % Singular value decomposition to update X
        [U,S,V] = svd(Y); 
        % XTemp = X; 
        X = U * wthresh(S,'s',tau) *V';
        % Update Y
        Y = Y+step*(M_masked-X);
        Y = mask.*Y;
        % Stopping criteria
        epsilon_X = norm(mask.*X - M_masked,'fro')*M_mask_inv;
        if epsilon_X<tol
            break;
        end
    end
    SVT_recon = X;
end