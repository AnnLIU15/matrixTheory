function [X_recon, obj_iter] = admm(M_masked, mask, A_l, B_l, args)
    beta = args.beta;
    beta_inv = 1 / args.beta;
    lambda = args.lambda;
    missing_val = ~mask;
    AB = A_l'*B_l;
    M_fro_inv = 1/norm(M_masked, 'fro');
    X_k = M_masked;
    W_k = X_k;
    Y_k = X_k;
    obj_val = zeros(args.inner_iter,1);
    for k = 1:args.inner_iter
        %% step 1 singular value shrinkage operator
        % https://www.cnblogs.com/kailugaji/p/14613210.html
        % wthresh(sigma,'s',1/beta)  equal to max(sigma - 1/beta,0)
        [U, sigma, V] = svd(W_k - beta_inv * Y_k);
        X_kp1 = U * wthresh(sigma,'s',beta_inv) * V';
        %% step 2
        W_kp1 = X_kp1 + beta_inv * (AB + Y_k);
        % Fix values at observed entries
        W_kp1 = W_kp1.* missing_val + M_masked;
        %% step 3
        Y_kp1 = Y_k + beta *(X_kp1 - W_kp1);
        epsilon_X = norm(X_kp1 - X_k,'fro') * M_fro_inv;
        if epsilon_X <= args.tol
            break
        end
        X_k = X_kp1;
        W_k = W_kp1;
        Y_k = Y_kp1;
        % eq(22) page 5 col 1
        obj_val(k) = sum(svd(X_k)) - trace(A_l*W_k*B_l') + ...
            lambda/2 * norm(X_k - W_k,'fro') ^2 + ...
            trace(Y_k'*(X_k-W_k));
            % lambda = 0.06
        
    end
    X_recon = X_kp1;
    obj_iter.k = k;
    obj_iter.obj_val = obj_val(1:k);
end