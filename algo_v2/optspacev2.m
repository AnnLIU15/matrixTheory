function M_approx= optspacev2(M_masked, mask, r,tau, max_iter, tol)
    % OPTSPACE 算法
    % 输入：
    %     M_masked: 待恢复的矩阵the observed matrix
    %     r: 希望恢复出来的矩阵的秩
    % 输出：
    %     M_approx: 恢复矩阵
    
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
    
    [U, S, V] = svds(M_trim, r);
    
    U = sqrt(m) * U;
    V = sqrt(n) * V;
    
    %S_0 = S(1:r, 1:r)/epsilon;  %运算中不需要
    %%
    M_E_fro_inv = 1 / norm(M_masked, 'fro');
    
    for k = 1: max_iter
        S = optimize_F_S(U, V, M_masked, mask);
        tmp_grad = (mask.* (U * S * V' - M_masked));
        grad_X = tmp_grad * V * S';
        grad_Y = tmp_grad' * U * S;
        
        t = tau;
        
        for i = 1: 50
            if(cost_F(U - t*grad_X, V - t*grad_Y,S,M_masked, mask) - cost_F(U, V,S,M_masked, mask) < ...
                    t * 0.5 * (norm(grad_X, 'fro')^2 + norm(grad_Y, 'fro')^2))
                break;
            end
            t = t/2;
        end
        
        U = U - t*grad_X;
        V = V - t*grad_Y;
        
        epsilon = norm(mask.*(M_masked - U*S*V'), 'fro') * M_E_fro_inv;
        % fprintf('%dth iteration: %e \n', k, epsilon);
        
        if (epsilon < tol)
            % fprintf('convergent');
            break;
        end
    end
    
    M_approx = U*S*V';
    end
    %%
    function S = optimize_F_S(U, V, M_masked, mask)
    % U, V: the input value of variable in the function F
    % M_E: the observed matrix
    % mask: the mask matrix
    % S：the matrix S that minimize the function F given value U, V
    % 想着化成A*S_ij = b的形式求解
    
    [~, k] = size(U);
    b = U'* M_masked * V;
    b = b(:); %拉成列向量
    
    A = zeros(k*k, k*k);
    for i = 1: k
        for j = 1:k
            index = (j-1) * k + i;
            coef = U' * (U(:,i) * V(:, j)' .* mask) * V;
            A(:, index) = coef(:);
        end
    end
    S = A\b;
    
    S = reshape(S, k, k);
    end
    %%
    function F = cost_F(U, V, S, M_masked, mask)
    %计算损失函数
    F = sum( sum( ( (U*S*V' - M_masked).*mask ).^2 ) )/2 ;
    end