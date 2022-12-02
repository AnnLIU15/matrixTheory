function [X_recon, obj_iter] = admmap(M_masked, mask, A_l, B_l, args)
    beta_k = args.beta;
    beta_max = args.beta_max;
    rho0 = args.rho0;
    kappa = args.kappa;
    [m,n] = size(M_masked);
    
    AB = A_l'*B_l;
    M_fro_inv = 1/norm(M_masked, 'fro');
    X_k = M_masked;
    W_k = X_k;
    Y_k = Aop(M_masked);
    C = zeros('like', Y_k);
    C(m+1:2*m,1+n:2*n) = M_masked;
    C_fro_inv = 1/norm(C,'fro');
    obj_val = zeros(args.inner_iter,1);
    for k = 1:args.inner_iter
        %% step 1 singular value shrinkage operator
        [U, sigma, V] = svd(W_k - 1/ beta_k * Y_k(1:m,1:n));
        X_kp1 = U * wthresh(sigma,'s',1/ beta_k) * V';
        %% step 2
        W_kp1 = X_kp1 + 1/ beta_k * (AB + Y_k(1:m,1:n)) + ...
                0.5 * 1/ beta_k * (mask .* ...
                  (beta_k*(M_masked - X_kp1) - ...
                  (AB + Y_k(1:m,1:n) + Y_k(m+1:2*m,n+1:2*n))));
        % Fix values at observed entries
        
        %% step 3
        Y_kp1 = Y_k + beta_k *(Aop(X_kp1) + Bop(W_kp1,mask) - C);
        if beta_k * C_fro_inv * max(...
            norm(X_kp1 - X_k,'fro'),norm(W_kp1 - W_k,'fro') ...
            )  < kappa
            rho = rho0;
        else
            rho = 1;
        end
        beta_kp1 = min(beta_max, rho * beta_k);
        epsilon_X = norm(X_kp1 - X_k,'fro') * M_fro_inv;
        if epsilon_X <= args.tol
            break
        end
        X_k = X_kp1;
        W_k = W_kp1;
        Y_k = Y_kp1;
        beta_k = beta_kp1;
        % eq(37) page 6 col 2
        obj_val(k) = sum(svd(X_k)) - trace(A_l*W_k*B_l') + ...
            sum(vec(Y_k .* (Aop(X_k) + Bop(W_k,mask) - C))) + ...
            beta_k * 0.5 * norm(Aop(X_k) + Bop(W_k,mask) - C,'fro')^2;
            % lambda = 0.06
        
    end
    X_recon = X_kp1;
    obj_iter.k = k;
    obj_iter.obj_val = obj_val(1:k);
end

function AX = Aop(X)
    [m, n] = size(X);
    AX = zeros(2*m, 2*n);
    AX(1:m, 1:n) = X;
end

function BX = Bop(X, mask)
    [m, n] = size(X);
    BX = zeros(2*m, 2*n);
    BX(1:m, 1:n) = -X;
    BX(m+1:2*m, n+1:2*n) = X .* mask;
end