%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% config loading
config = ReadYaml('../configs/fig3.yaml');
m = config.base.m; n =config.base.n;
sigma0 = config.base.sigma0; rank_length = config.base.rank_length;
tol = config.base.tol; max_iter = config.base.max_iter;
if ~exist(config.base.save_dir, 'dir'), mkdir(config.base.save_dir); end
obs_p_list = 0.6:0.1:0.9;
%% initial the output mat
SVT_erec     = zeros(4,rank_length);
SVP_erec     = zeros(4,rank_length);
Opt_erec     = zeros(4,rank_length);
admm_erec    = zeros(4,rank_length);
apgl_erec    = zeros(4,rank_length);
admmap_erec  = zeros(4,rank_length);
%% rank_level loop
for rank_level = 1:rank_length
    M = randn(m, rank_level) * randn(rank_level, n);
    Z = randn(m, n);
    B = M + sigma0 * Z;
    %% obs_percent loop
    for idx_obs_p = 1:4
        obs_p = obs_p_list(idx_obs_p);
        mask = zeros(m, n);
        mask(randperm(m*n, round(obs_p*m*n))) = 1;
        missing = ~mask;
        %% masked noise synethetic data
        X_masked = B .* mask;
        % SVT
        tau = sqrt(m*n);
        step = 1.2 * obs_p;
        SVT_recon = SVT(X_masked, mask, tau, step, max_iter, tol);

        % SVP
        step = 1/obs_p/sqrt(max_iter);
        SVP_recon = SVP(X_masked,mask,step,rank_level,max_iter,tol);

        % Optspace
        Opt_recon = OPTSPACE(X_masked, rank_level);

        % TNNR-admm
        admm_recon = admm_mat(M, X_masked, mask, rank_level, config);

        % TNNR-apgl
        [apgl_recon, ~] = apgl_mat(M, X_masked, mask, rank_level, config);

        % TNNR-admmap
        [admmap_recon, ~] = admmap_mat(M, X_masked, mask, rank_level, config);

        SVT_erec(idx_obs_p,rank_level)    = norm(vec((M-SVT_recon).*missing), 'fro');
        SVP_erec(idx_obs_p,rank_level)    = norm(vec((M-SVP_recon).*missing), 'fro');
        Opt_erec(idx_obs_p,rank_level)    = norm(vec((M-Opt_recon).*missing), 'fro');
        admm_erec(idx_obs_p,rank_level)   = norm(vec((M-admm_recon).*missing), 'fro');
        apgl_erec(idx_obs_p,rank_level)   = norm(vec((M-apgl_recon).*missing), 'fro');
        admmap_erec(idx_obs_p,rank_level) = norm(vec((M-admmap_recon).*missing), 'fro');
    end
end
figure;
for idx_obs_p = 1:4
    subplot(2,2,idx_obs_p)
    %% plot and save result
    obs_p = obs_p_list(idx_obs_p) * 100;
    idx_x = 1:rank_length;
    plot(idx_x, SVT_erec(idx_obs_p,:), '-sc', ...
         idx_x, SVP_erec(idx_obs_p,:), '-om', ...
         idx_x, Opt_erec(idx_obs_p,:), '-dg', ...
         idx_x, admm_erec(idx_obs_p,:), '-^k',...
         idx_x, apgl_erec(idx_obs_p,:), '->b',...
         idx_x, admmap_erec(idx_obs_p,:), '-+r');
    xlim([1, rank_length]);
    set(gcf,'color','none'); % background color -> none
    set(gca,'color','none'); % axis color -> none
    xlabel('Rank');
    ylabel('Total reconstruction error');
    legend('SVT', 'SVP', 'Optspace', 'admm', 'apgl', 'admmap','Location','northwest');
    legend('boxoff');
    title([num2str(obs_p),'% observed'])
end
save_path = [config.base.save_dir,'/',config.base.fig_name];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
close all;