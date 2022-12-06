function M_approx= optspace(M_masked, mask, r, tau, max_iter, tol)
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
    
    %% Projection
    [m, n] = size(M_masked);
    det_E = sum(mask,'all');
    
    %% trimming
    
    avg_row = 2 * det_E/ m;
    avg_col = 2 * det_E/ n;
    
    M_row = sum(mask, 2);
    M_col = sum(mask, 1);
    
    index_row = M_row > avg_row;
    index_col = M_col > avg_col;
    
    M_trim = M_masked;
    M_trim(index_row,:) = 0;
    M_trim(:, index_col) = 0;
    %% Project \tilde{M^mask} to T_r(\tilde{M^mask})
    
    [X, S, Y] = svds(M_trim, r);
    epsilon_inv = m*n/det_E;
    X = sqrt(m) * X;
    Y = sqrt(n) * Y;
    S = S * epsilon_inv;
    
    %%
    M_E_fro_inv = 1 / norm(M_masked, 'fro');
    
    for k = 1: max_iter
        S = optimize_F_S(X, Y, M_masked, mask);
        tmp_grad = (mask.* (X * S * Y' - M_masked));
        grad_X = tmp_grad * Y * S';
        grad_Y = tmp_grad' * X * S;
        t = tau;
        w_k = 0.5 * (norm(grad_X, 'fro')^2 + norm(grad_Y, 'fro')^2);
        t_w_k = tau * w_k;
        for i = 1: 50
            if(costFunc(X - t*grad_X, Y - t*grad_Y,S,M_masked, mask,0) - costFunc(X, Y,S,M_masked, mask,0) <  t_w_k)
                break;
            end
            t_w_k = t_w_k/2;
        end
        
        X = X - t*grad_X;
        Y = Y - t*grad_Y;
        if (norm(mask.*(M_masked - X*S*Y'), 'fro') * M_E_fro_inv < tol)
            break;
        end
    end
    
    M_approx = X*S*Y';
end


function S = optimize_F_S(X, Y, M_masked, mask)
    % X, Y: the input value of variable in the function F
    % M_E: the observed matrix
    % mask: the mask matrix
    % S: the matrix S that minimize the function F given value X, Y
    % solve A*S_ij = b problem

    [~, k] = size(X);
    b = X'* M_masked * Y;
    b = b(:); 
    A = zeros(k*k, k*k);
    for i = 1: k
        for j = 1:k
            index = (j-1) * k + i;
            coef = X' * (X(:,i) * Y(:, j)' .* mask) * Y;
            A(:, index) = coef(:);
        end
    end
    S = A\b;

    S = reshape(S, k, k);
end
%%
function F = costFunc(X, Y, S, M_masked, mask, lambda)
    %% Loss function
    XSY = X*S*Y';
    % loss = 0.5 * norm(M_masked - mask .* XSY,'fro')^2;
    loss = 0.5 * sum((M_masked - mask .* XSY).^2,'all'); 
    % reg  = lambda * 0.5 * norm((~mask).*XSY,'fro')^2;
    % reg = lambda * 0.5 * sum(((~mask).*XSY).^2,'all');
    % F = loss + reg;
    F = loss;
end