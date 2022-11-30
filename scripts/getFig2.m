%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% config loading
config = ReadYaml('../configs/fig2.yaml');
m = config.base.m; n =config.base.n;
r0 = config.base.r0; max_iter = config.base.max_iter;
tol = config.base.tol;
if ~exist(config.base.save_dir, 'dir'), mkdir(config.base.save_dir); end
obs_p_list = 0.6:0.1:0.9;
%% Generating the matrix (only need once)
M = randn(m, r0) * randn(r0, n);
Z = randn(m, n);
%% initial the output mat
SVT_erec     = zeros(4,10);
SVP_erec     = zeros(4,10);
Opt_erec     = zeros(4,10);
admm_erec    = zeros(4,10);
apgl_erec    = zeros(4,10);
admmap_erec  = zeros(4,10);
%% sigma_level loop
for idx_sigma_level = 1:10
    cur_sigma = idx_sigma_level / 10;
    B = M + cur_sigma * Z;
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
        SVP_recon = SVP(X_masked,mask,step,r0,max_iter,tol);

        % Optspace
        Opt_recon = OPTSPACE(X_masked, r0);

        % TNNR-admm
        admm_recon = admm_mat(M, X_masked, mask, r0, config);

        % TNNR-apgl
        [apgl_recon, ~] = apgl_mat(M, X_masked, mask, r0, config);

        % TNNR-admmap
        [admmap_recon, ~] = admmap_mat(M, X_masked, mask, r0, config);

        SVT_erec(idx_obs_p,idx_sigma_level)    = norm(vec((M-SVT_recon).*missing), 'fro');
        SVP_erec(idx_obs_p,idx_sigma_level)    = norm(vec((M-SVP_recon).*missing), 'fro');
        Opt_erec(idx_obs_p,idx_sigma_level)    = norm(vec((M-Opt_recon).*missing), 'fro');
        admm_erec(idx_obs_p,idx_sigma_level)   = norm(vec((M-admm_recon).*missing), 'fro');
        apgl_erec(idx_obs_p,idx_sigma_level)   = norm(vec((M-apgl_recon).*missing), 'fro');
        admmap_erec(idx_obs_p,idx_sigma_level) = norm(vec((M-admmap_recon).*missing), 'fro');
    end
end
figure;
for idx_obs_p = 1:4
    subplot(2,2,idx_obs_p)
    %% plot and save result
    obs_p = obs_p_list(idx_obs_p) * 100;
    idx_x = 0.1:0.1:1;
    plot(idx_x, SVT_erec(idx_obs_p,:), '-sc', ...
         idx_x, SVP_erec(idx_obs_p,:), '-om', ...
         idx_x, Opt_erec(idx_obs_p,:), '-dg', ...
         idx_x, admm_erec(idx_obs_p,:), '-^k',...
         idx_x, apgl_erec(idx_obs_p,:), '->b',...
         idx_x, admmap_erec(idx_obs_p,:), '-+r');
    xlim([0.1 1]);
    set(gcf,'color','none'); % background color -> none
    set(gca,'color','none'); % axis color -> none
    xlabel('Noise level');
    ylabel('Total reconstruction error');
    legend('SVT', 'SVP', 'Optspace', 'admm', 'apgl', 'admmap','Location','northwest');
    legend('boxoff');
    title([num2str(obs_p),'% observed'])
end
print(gcf, '-dpdf',[config.base.save_dir,'/fig2.pdf'])
print(gcf,[config.base.save_dir,'/fig2.jpeg'] ,'-djpeg','-r300')
close all;