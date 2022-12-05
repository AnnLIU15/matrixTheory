%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% args loading
args = ReadYaml('../configs/fig3v2.yaml');
m = args.base.m; n =args.base.n;
sigma0 = args.base.sigma0; rank_length = args.base.rank_length;
max_iter = args.base.max_iter;
tol = args.base.tol; 
obs_p_list = 0.6:0.1:0.9;

%% initial the output mat
SVT_erec     = zeros(4,rank_length);
SVP_erec     = zeros(4,rank_length);
Opt_erec     = zeros(4,rank_length);
admm_erec    = zeros(4,rank_length);
apgl_erec    = zeros(4,rank_length);
admmap_erec  = zeros(4,rank_length);
datetime.setDefaultFormats('default',"yyyy/MM/dd HH:mm:ss:SSS")

%% rank_level loop
for rank_level = 1:rank_length
    M = randn(m, rank_level) * randn(rank_level, n);
    Z = randn(m, n);
    B = M + sigma0 * Z;
    fprintf("cur rank level: %d\n",rank_level);
    args.min_R = rank_level;
    args.max_R = rank_level;
    %% obs_percent loop
    for idx_obs_p = 1:4
        obs_p = obs_p_list(idx_obs_p);
        fprintf("\tcur obs level: %.2f\n",obs_p);
        mask = zeros(m, n);
        mask(randperm(m*n, round(obs_p*m*n))) = 1;
        missing = ~mask;
        %% masked noise synethetic data
        M_masked = B .* mask;
        % SVT
        cur_time = datetime('now');
        fprintf("\t\tStart SVT: %s ->",datestr(cur_time));
        tau = (m*n)^args.base.SVT_power;
        step = 1.2 * obs_p;
        SVT_recon = SVT(M_masked, mask, tau, step, max_iter, tol);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        % SVP
        cur_time = datetime('now');
        fprintf("\t\tStart SVP: %s ->",datestr(cur_time));
        step = 1/obs_p/sqrt(max_iter);
        SVP_recon = SVP(M_masked,mask,step,rank_level,max_iter,tol);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        % Optspace
        cur_time = datetime('now');
        opt_tau = 1e-2;
        fprintf("\t\tStart Optspace: %s ->",datestr(cur_time));
        Opt_recon = optspacev2(M_masked, mask, rank_level, opt_tau, max_iter, tol);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        % TNNR-admm
        cur_time = datetime('now');
        fprintf("\t\tStart TNNR-ADMM: %s ->",datestr(cur_time));
        [admm_recon, ~] = tnnr_recon(M_masked, mask, 'admm', args);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        admm_recon = reshape(admm_recon,m,n);
        
        % TNNR-apgl
        cur_time = datetime('now');
        fprintf("\t\tStart TNNR-AGPL: %s ->",datestr(cur_time));
        [apgl_recon, ~] = tnnr_recon(M_masked, mask, 'apgl', args);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        apgl_recon = reshape(apgl_recon,m,n);
        
        % TNNR-admmap
        cur_time = datetime('now');
        fprintf("\t\tStart TNNR-ADMMAP: %s ->",datestr(cur_time));
        [admmap_recon, ~] = tnnr_recon(M_masked, mask, 'admmap', args);
        fprintf(" %.3fs\n",getMSecDiff(cur_time));
        admmap_recon = reshape(admmap_recon,m,n);
        
        SVT_erec(idx_obs_p,rank_level)    = norm(vec((M-SVT_recon).*missing), 'fro');
        SVP_erec(idx_obs_p,rank_level)    = norm(vec((M-SVP_recon).*missing), 'fro');
        Opt_erec(idx_obs_p,rank_level)    = norm(vec((M-Opt_recon).*missing), 'fro');
        admm_erec(idx_obs_p,rank_level)   = norm(vec((M-admm_recon).*missing), 'fro');
        apgl_erec(idx_obs_p,rank_level)   = norm(vec((M-apgl_recon).*missing), 'fro');
        admmap_erec(idx_obs_p,rank_level) = norm(vec((M-admmap_recon).*missing), 'fro');
    end
end
%% plot the result
figure;
idx_x = 1:rank_length;
for idx_obs_p = 1:4
    subplot(2,2,idx_obs_p)
    % plot subfig with different noise level
    obs_p = obs_p_list(idx_obs_p) * 100;
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
if ~exist(args.base.save_dir, 'dir'), mkdir(args.base.save_dir); end

% saving
save_path = [args.base.save_dir,'/',args.base.fig_name];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
close all;