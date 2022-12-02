function X_l_recon_rank_list = tnnr_recon(M_masked, mask, chosen_algo, args)
    %--------------------------------------------------------------------------
    % 
    %--------------------------------------------------------------------------
    %     main part of TNNR algorithm (via ADMM/APGL/ADMMAP)
    % 
    %     Inputs:
    %         M_masked             --- original image with masked
    %         mask                 --- index matrix of known elements
    %         args                 --- struct of parameters
    %         chosen_algo          --- current algo (type str)
    %     Outputs: 
    %         tnnr_algo_recon      --- result of algorithm
    %         X_rec                --- recovered image under the best rank
    %--------------------------------------------------------------------------
    min_R = args.min_R;        % minimum rank of chosen image
    max_R = args.max_R;        % maximum rank of chosen image
    rank_length = max_R - min_R + 1;
    [m,n,dim] = size(M_masked);
    M_fro = norm(M_masked, 'fro');

    if chosen_algo == "admm"
        algo = @(X, mask, A_l, B_l, args)admm(X, mask, A_l, B_l, args);
    elseif chosen_algo == "apgl"
        algo = @(X, mask, A_l, B_l, args)apgl(X, mask, A_l, B_l, args);
    elseif chosen_algo == "admmap"
        algo = @(X, mask, A_l, B_l, args)admmap(X, mask, A_l, B_l, args);
    else
        error('unknow algo %s',chosen_algo);
    end
    if dim == 1
        X_l_recon_rank_list = zeros(rank_length, m, n);
    else
        X_l_recon_rank_list = zeros(rank_length, m, n, dim);
    end
    for rank_idx = min_R : max_R
        X_l = M_masked;
        for outer_idx = 1 : args.outer_iter
            %% STEP 1 given X_l
            [U_l, ~, V_l] = svd(X_l);
            A_l = U_l(:,1:rank_idx)';
            B_l = V_l(:,1:rank_idx)';
            %% choose a function solve step 2
            [X_lp1, ~] = algo(M_masked,mask,A_l,B_l,args);
            early_stop_factor = norm(X_lp1 - X_l,'fro') / M_fro;
            display([num2str(outer_idx),'\t',num2str(early_stop_factor)]);
            if early_stop_factor <= args.outer_tol
                break;
            end
            X_l = X_lp1;
        end
        X_l = clip(X_l,0,255);
        X_l_recon_rank_list(rank_idx,:,:,1) = X_l;
        figure;
        imshow(X_l/255,[])
    end
end