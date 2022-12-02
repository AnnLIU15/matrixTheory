function [X_rec, iter] = admmap_mat(M, X, mask, R, para)
%--------------------------------------------------------------------------
%     Inputs:
%         M                    --- original matrix
%         X                    --- matrix after add noise and mask
%         mask                 --- index matrix of known elements
%         R                    --- matrix rank
%         para                 --- struct of parameters  
%
%     Outputs: 
%         X_rec           --- recovered matrix
%--------------------------------------------------------------------------
max_iter = para.outer_iter;   % maximum number of outer iteration
tol      = para.outer_tol;    % tolerance of outer iteration
iter = 0;
M_fro = norm(M, 'fro');
for i = 1 : max_iter
    fprintf('iter %d\n', i);
    last_X = X;
    [U, ~, V] = svd(X);
    A = U(:, 1:R)'; B = V(:, 1:R)';
    [X, c] = admmapAXB(A, B, X, M, mask, para);
    
    delta = norm(X - last_X, 'fro') / M_fro;
    
    iter = iter + c;
    
    if delta < tol
        % fprintf('converged at iter: %d\n', i);
        break ;
    end
end
X_rec = X;
end