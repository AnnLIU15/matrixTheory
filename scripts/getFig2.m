%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% args loading
args = ReadYaml('../configs/fig2.yaml');
m = args.base.m; n =args.base.n;
r0 = args.base.r0; max_iter = args.base.max_iter;
tol = args.base.tol; noise_length= args.base.noise_length;
obs_p_list = 0.6:0.1:0.9;
%% Generating the matrix (only need once)
M = randn(m, r0) * randn(r0, n);
Z = randn(m, n);
%% initial the output mat
SVT_erec     = zeros(4,noise_length);
SVP_erec     = zeros(4,noise_length);
Opt_erec     = zeros(4,noise_length);
admm_erec    = zeros(4,noise_length);
apgl_erec    = zeros(4,noise_length);
admmap_erec  = zeros(4,noise_length);
% set time's format
datetime.setDefaultFormats('default',"yyyy/MM/dd HH:mm:ss:SSS")

%% sigma_level loop
for idx_sigma_level = 1:noise_length
    cur_sigma = idx_sigma_level / 10;
    fprintf("cur noise level: %.2f\n",cur_sigma);
    B = M + cur_sigma * Z;
    %% obs_percent loop
    for idx_obs_p = 1:4
        %% get current mask and masked noise synethetic data
        obs_p = obs_p_list(idx_obs_p);
        fprintf("\tcur obs level: %.2f\n",obs_p);
        mask = zeros(m, n);
        mask(randperm(m*n, round(obs_p*m*n))) = 1;
        missing = ~mask;
        M_masked = B .* mask;

        %% SVT
        cur_time = datetime('now');
        fprintf("\t\tStart SVT: %s ->",datestr(cur_time));
        tau = (m*n)^args.base.SVT_power;
        step = 1.2 * obs_p;
        SVT_recon = SVT(M_masked, mask, tau, step, max_iter, tol);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));

        %% SVP
        cur_time = datetime('now');
        fprintf("\t\tStart SVP: %s ->",datestr(cur_time));
        step = 1/(1+args.base.delta_2k); %\delta_{2k}<1/3
        SVP_recon = SVP(M_masked,mask,step,r0,max_iter,tol);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));

        %% Optspace
        cur_time = datetime('now');
        fprintf("\t\tStart Optspace: %s ->",datestr(cur_time));
        Opt_recon = optspace(M_masked, mask, r0, args.base.opt_tau, max_iter, tol);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));

        %% TNNR-admm
        cur_time = datetime('now');
        fprintf("\t\tStart TNNR-ADMM: %s ->",datestr(cur_time));
        [admm_recon, ~] = tnnr_recon(M_masked, mask, 'admm', args);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        admm_recon = reshape(admm_recon,m,n);

        %% TNNR-apgl
        cur_time = datetime('now');
        fprintf("\t\tStart TNNR-AGPL: %s ->",datestr(cur_time));
        [apgl_recon, ~] = tnnr_recon(M_masked, mask, 'apgl', args);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        apgl_recon = reshape(apgl_recon,m,n);

        %% TNNR-admmap
        cur_time = datetime('now');
        fprintf("\t\tStart TNNR-ADMMAP: %s ->",datestr(cur_time));
        [admmap_recon, ~] = tnnr_recon(M_masked, mask, 'admmap', args);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        admmap_recon = reshape(admmap_recon,m,n);
        %% Calculate the F-norm
        SVT_erec(idx_obs_p,idx_sigma_level)    = norm(vec((M-SVT_recon).*missing), 'fro');
        SVP_erec(idx_obs_p,idx_sigma_level)    = norm(vec((M-SVP_recon).*missing), 'fro');
        Opt_erec(idx_obs_p,idx_sigma_level)    = norm(vec((M-Opt_recon).*missing), 'fro');
        admm_erec(idx_obs_p,idx_sigma_level)   = norm(vec((M-admm_recon).*missing), 'fro');
        apgl_erec(idx_obs_p,idx_sigma_level)   = norm(vec((M-apgl_recon).*missing), 'fro');
        admmap_erec(idx_obs_p,idx_sigma_level) = norm(vec((M-admmap_recon).*missing), 'fro');
    end
end
%% plot the result
figure;
for idx_obs_p = 1:4
    subplot(2,2,idx_obs_p)
    % plot subfig with different noise level
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

%% saving (pay attention to the relative path)
if ~exist(args.base.save_dir, 'dir'), mkdir(args.base.save_dir); end
save_path = [args.base.save_dir,'/',args.base.fig_name];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
close all;