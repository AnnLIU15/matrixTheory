function [X_recon, obj_iter] = apgl(M_masked, mask, A_l, B_l, args)
    t_k = 1;
    lambda = args.lambda;
    M_fro_inv = 1/norm(M_masked, 'fro');
    AB = A_l'*B_l;
    X_k = M_masked;
    Y_k = X_k;
    obj_val = zeros(args.inner_iter,1);
    for k = 1:args.inner_iter
        %% step 1 singular value shrinkage operator
        % https://www.cnblogs.com/kailugaji/p/14613210.html
        % wthresh(sigma,'s',1/beta)  equal to max(sigma - 1/beta,0)
        [U, sigma, V] = svd(Y_k + t_k*(AB - lambda*(Y_k.*mask - M_masked)));
        X_kp1 = U * wthresh(sigma,'s',t_k) * V';
        %% step 2
        t_kp1 = 0.5 * (1 + sqrt(1 + 4*(t_k^2)));
        Y_kp1 = X_kp1 + (t_k - 1)/t_kp1 *(X_kp1 - X_k);
        epsilon_X = norm(X_kp1 - X_k,'fro') * M_fro_inv;
        if epsilon_X <= args.tol
            break
        end
        X_k = X_kp1;
        Y_k = Y_kp1;
        t_k = t_kp1;
        % eq(31) page 5 col 2
        obj_val(k) = sum(svd(X_k)) - trace(A_l*X_k*B_l') + ...
            lambda/2 * norm(X_k.*mask - M_masked,'fro') ^2;
            % lambda = 0.06
    end
    X_recon = X_kp1;
    obj_iter.k = k;
    obj_iter.obj_val = obj_val(1:k);
end