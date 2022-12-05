function M_approx= optspacev2(M_masked, mask, r,tau, max_iter, tol)
    % OPTSPACE 绠楁硶
    % 杈撳叆锛?
    %     M_masked: 寰呮仮澶嶇殑鐭╅樀the observed matrix
    %     r: 甯屾湜鎭㈠鍑烘潵鐨勭煩闃电殑绉?
    % 杈撳嚭锛?
    %     M_approx: 鎭㈠鐭╅樀
    
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
    
    %S_0 = S(1:r, 1:r)/epsilon;  %杩愮畻涓笉闇?瑕?
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
    % S锛歵he matrix S that minimize the function F given value U, V
    % 鎯崇潃鍖栨垚A*S_ij = b鐨勫舰寮忔眰瑙?
    
    [~, k] = size(U);
    b = U'* M_masked * V;
    b = b(:); %鎷夋垚鍒楀悜閲?
    
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
    %璁＄畻鎹熷け鍑芥暟
    F = sum( sum( ( (U*S*V' - M_masked).*mask ).^2 ) )/2 ;
    end