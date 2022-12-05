%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% config loading
args = ReadYaml('../configs/fig5.yaml');
tol = args.base.tol; max_iter = args.base.max_iter;
obs_p_list = 0.4:0.1:0.9;
obs_p_length = length(obs_p_list);
%% initial the output mat
SVT_psnr     = zeros(obs_p_length,1);
SVP_psnr     = zeros(obs_p_length,1);
Opt_psnr     = zeros(obs_p_length,1);
admm_psnr    = zeros(obs_p_length,1);
apgl_psnr    = zeros(obs_p_length,1);
admmap_psnr  = zeros(obs_p_length,1);

rgb_fig = imread(args.base.img_path);
gray_fig = rgb2gray(rgb_fig);
[m,n,dim] = size(gray_fig);
% SVT_rec_fig    = zeros(m,n);
% SVP_rec_fig    = zeros(m,n);
% Opt_rec_fig    = zeros(m,n);
% admm_rec_fig   = zeros(m,n);
% apgl_rec_fig   = zeros(m,n);
% admmap_rec_fig = zeros(m,n);

%% obs_percent loop
gray_fig_d = double(gray_fig);
args.base.SVT_power = 0.85;
for idx_obs_p = 1:obs_p_length
% for idx_obs_p = 3:3 
    obs_p = obs_p_list(idx_obs_p);
    fprintf("obs_p : %f \n", obs_p);
    mask = zeros(m, n);
    mask(randperm(m*n, round(obs_p*m*n))) = 1;

    missing = ~mask;
    % masked noise synethetic data
    X_masked = gray_fig_d .* mask;

%     % SVT
    tau = (m*n)^args.base.SVT_power;
    step = 1.2 * obs_p;
    SVT_recon = SVT(X_masked, mask, tau, step, max_iter, tol);
    
    [~, SVT_psnr(idx_obs_p)] = PSNR(gray_fig_d, SVT_recon, missing);
%     
%     
%     for rank_level = args.min_R:args.max_R
%         fprintf("rank : %d \n", rank_level);
%         % SVP
%         step = 1/obs_p/sqrt(max_iter);
%         SVP_recon = SVP(X_masked, mask, step, rank_level, max_iter, tol);
%         % Optspace
%         tau = 1;
%         Opt_recon = optspacev2(X_masked, mask, rank_level, tau, max_iter, tol);
%         
%         [~, tmp_psnr] = PSNR(gray_fig_d, SVP_recon, missing);
%         if tmp_psnr > SVP_psnr(idx_obs_p)
%             r_svp(idx_obs_p) = rank_level;
%         end
%         SVP_psnr(idx_obs_p)    = max(tmp_psnr,SVP_psnr(idx_obs_p));
%         % if SVP_psnr(idx_obs_p) == max(SVP_psnr)
%         %     SVP_rec_fig = SVP_recon;
%         % end
%         [~, tmp_psnr] = PSNR(gray_fig_d, Opt_recon, missing);
%         if tmp_psnr > Opt_psnr(idx_obs_p)
%             r_opt(idx_obs_p) = rank_level;
%         end
%         Opt_psnr(idx_obs_p)    = max(tmp_psnr,Opt_psnr(idx_obs_p));
%         % if Opt_psnr(idx_obs_p) == max(Opt_psnr)
%         %     Opt_rec_fig = Opt_recon;
%         % end
%     end
%     fprintf("end SVP OPT \n");

%     chosen_algo = "admm";
%     [X_admm_recon_rank_list, ] = tnnr_recon(X_masked, mask, chosen_algo, args);
%     chosen_algo = "apgl";
%     [X_apgl_recon_rank_list, ] = tnnr_recon(X_masked, mask, chosen_algo, args);
%     chosen_algo = "admmap";
%     [X_admmap_recon_rank_list, ] = tnnr_recon(X_masked, mask, chosen_algo, args);

%     for rank_level = args.min_R:args.max_R
%         admm_recon = reshape( X_admm_recon_rank_list(rank_level - args.min_R + 1, :, :, :), m,n,dim);
%         [~, tmp_psnr] = PSNR(gray_fig_d, admm_recon, missing);
%         if tmp_psnr > admm_psnr(idx_obs_p)
%             r_admm(idx_obs_p) = rank_level;
%         end
%         admm_psnr(idx_obs_p) = max(tmp_psnr,admm_psnr(idx_obs_p));
%         
%         apgl_recon = reshape( X_apgl_recon_rank_list(rank_level - args.min_R + 1, :, :, :), m,n,dim);
%         [~, tmp_psnr] = PSNR(gray_fig_d, apgl_recon, missing);
%         if tmp_psnr > apgl_psnr(idx_obs_p)
%             r_apgl(idx_obs_p) = rank_level;
%         end
%         apgl_psnr(idx_obs_p) = max(tmp_psnr,apgl_psnr(idx_obs_p));
%         
%         admmap_recon = reshape( X_admmap_recon_rank_list(rank_level - args.min_R + 1, :, :, :), m,n,dim);
%         [~, tmp_psnr] = PSNR(gray_fig_d, admmap_recon, missing);
%         if tmp_psnr > admmap_psnr(idx_obs_p)
%             r_admmap(idx_obs_p) = rank_level;
%         end
%         admmap_psnr(idx_obs_p) = max(tmp_psnr,admmap_psnr(idx_obs_p));
%     end

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


% plot and save result
figure;
obs_p = obs_p_list(idx_obs_p) * 100;
% plot(obs_p_list, SVT_psnr, '-sc', ...
%      obs_p_list, SVP_psnr, '-om', ...
%      obs_p_list, Opt_psnr, '-dg')
plot(obs_p_list, SVT_psnr, '-sc', ...
     obs_p_list, SVP_psnr, '-om', ...
     obs_p_list, Opt_psnr, '-dg', ...
     obs_p_list, admm_psnr, '-^k',...
     obs_p_list, apgl_psnr, '->b',...
     obs_p_list, admmap_psnr, '-+r');
xlim([min(obs_p_list),max(obs_p_list)]);
% set(gcf,'color','none'); % background color -> none
% set(gca,'color','none'); % axis color -> none
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