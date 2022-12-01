%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% config loading
config = ReadYaml('../configs/fig5.yaml');
tol = config.base.tol; max_iter = config.base.max_iter;
obs_p_list = 0.4:0.1:0.9;
obs_p_length = length(obs_p_list);
%% initial the output mat
SVT_psnr     = zeros(obs_p_length);
SVP_psnr     = zeros(obs_p_length);
Opt_psnr     = zeros(obs_p_length);
admm_psnr    = zeros(obs_p_length);
apgl_psnr    = zeros(obs_p_length);
admmap_psnr  = zeros(obs_p_length);
rgb_fig = imread(config.base.img_path);
gray_fig = rgb2gray(rgb_fig);
[m,n] = size(gray_fig);
gray_fig_d = double(gray_fig);
SVT_rec_fig    = zeros(m,n);
SVP_rec_fig    = zeros(m,n);
Opt_rec_fig    = zeros(m,n);
admm_rec_fig   = zeros(m,n);
apgl_rec_fig   = zeros(m,n);
admmap_rec_fig = zeros(m,n);
%% obs_percent loop
for idx_obs_p = 1:obs_p_length
    obs_p = obs_p_list(idx_obs_p);
    mask = zeros(m, n);
    mask(randperm(m*n, round(obs_p*m*n))) = 1;

    missing = ~mask;
    %% masked noise synethetic data
    X_masked = gray_fig_d .* mask;
    % SVT
    tau = sqrt(m*n);
    step = 1.2 * obs_p;
    SVT_recon = SVT(X_masked, mask, tau, step, max_iter, tol);
    [~, SVT_psnr(idx_obs_p)] = PSNR(gray_fig_d, SVT_recon, mask);
    for rank_level = 1:20
        % SVP
        step = 1/obs_p/sqrt(max_iter);
        SVP_recon = SVP(X_masked,mask,step,rank_level,max_iter,tol);
        % Optspace
        Opt_recon = OPTSPACE(X_masked, rank_level);
        
        [~, tmp_psnr] = PSNR(gray_fig_d, SVP_recon, mask);
        SVP_psnr(idx_obs_p)    = max(tmp_psnr,SVP_psnr(idx_obs_p));
        % if SVP_psnr(idx_obs_p) == max(SVP_psnr)
        %     SVP_rec_fig = SVP_recon;
        % end
        [~, tmp_psnr] = PSNR(gray_fig_d, Opt_recon, mask);
        Opt_psnr(idx_obs_p)    = max(tmp_psnr,Opt_psnr(idx_obs_p));
        % if Opt_psnr(idx_obs_p) == max(Opt_psnr)
        %     Opt_rec_fig = Opt_recon;
        % end
    end
    % TNNR-admm
    [admm_recon,~] = admm_pic(gray_fig_d, mask, config);
    % TNNR-apgl
    [apgl_recon, ~] = apgl_pic(gray_fig_d, mask, config);
    % TNNR-admmap
    [admmap_recon, ~] = admmap_pic(gray_fig_d, mask, config);
    admm_psnr  (idx_obs_p) = admm_recon.psnr_best;
    apgl_psnr  (idx_obs_p) = apgl_recon.psnr_best;
    admmap_psnr(idx_obs_p) = admmap_recon.psnr_best;
    % if SVT_psnr(idx_obs_p) == max(SVT_psnr)
    %     SVT_rec_fig = SVT_recon;
    % end
    % if admm_recon.psnr_best == max(admm_psnr  (idx_obs_p))
    %     admm_rec_fig = admm_recon.X_best_rec;
    % end
    % if apgl_recon.psnr_best == max(apgl_psnr  (idx_obs_p))
    %     apgl_rec_fig = apgl_recon.X_best_rec;
    % end
    
    % if admmap_recon.psnr_best == max(admmap_psnr  (idx_obs_p))
    %     admmap_rec_fig = admmap_recon.X_best_rec;
    % end
    
end

figure;

%% plot and save result
obs_p = obs_p_list(idx_obs_p) * 100;
plot(obs_p_list, SVT_psnr, '-sc', ...
     obs_p_list, SVP_psnr, '-om', ...
     obs_p_list, Opt_psnr, '-dg')
plot(obs_p_list, SVT_psnr, '-sc', ...
     obs_p_list, SVP_psnr, '-om', ...
     obs_p_list, Opt_psnr, '-dg', ...
     obs_p_list, admm_psnr, '-^k',...
     obs_p_list, apgl_psnr, '->b',...
     obs_p_list, admmap_psnr, '-+r');
xlim([min(obs_p_list),max(obs_p_list)]);
set(gcf,'color','none'); % background color -> none
set(gca,'color','none'); % axis color -> none
xlabel('Observed Ratio');
ylabel('PSNR');
legend('SVT', 'SVP', 'Optspace','admm','apgl','admmap','Location','northwest');
legend('boxoff');
%% saving the results
if ~exist(config.base.save_dir, 'dir'), mkdir(config.base.save_dir); end
save_path = [config.base.save_dir,'/',config.base.fig_name];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')




close all;