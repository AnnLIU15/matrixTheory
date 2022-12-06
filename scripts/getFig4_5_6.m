%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% args loading
args = ReadYaml('../configs/fig4_5_6.yaml');
tol = args.base.tol; max_iter = args.base.max_iter;
obs_p_fig = args.obs_p*10-3;
obs_p_list = 0.4:0.1:0.9;
obs_p_length = length(obs_p_list);
%% initial the output mat
SVT_psnr     = zeros(obs_p_length,1);
SVP_psnr     = zeros(obs_p_length,1);
Opt_psnr     = zeros(obs_p_length,1);
admm_psnr    = zeros(obs_p_length,1);
apgl_psnr    = zeros(obs_p_length,1);
admmap_psnr  = zeros(obs_p_length,1);

admm_iter    = zeros(obs_p_length,1);
apgl_iter    = zeros(obs_p_length,1);
admmap_iter  = zeros(obs_p_length,1);

rgb_fig = imread(args.base.img_path);
gray_fig = rgb2gray(rgb_fig);
gray_fig_d = double(gray_fig);
[m,n,dim] = size(gray_fig);

masked_fig     = zeros(obs_p_length,m,n);
SVT_rec_fig    = zeros(obs_p_length,m,n);
SVP_rec_fig    = zeros(obs_p_length,m,n);
Opt_rec_fig    = zeros(obs_p_length,m,n);
admm_rec_fig   = zeros(obs_p_length,m,n);
apgl_rec_fig   = zeros(obs_p_length,m,n);
admmap_rec_fig = zeros(obs_p_length,m,n);
%%
args.beta = 1;
args.lambda = 5e-2;
for idx_obs_p = 1:obs_p_length
% for idx_obs_p = 3:3
    obs_p = obs_p_list(idx_obs_p);
    
    fprintf("obs_p : %f \n", obs_p);
    
    mask = zeros(m, n);
    mask(randperm(m*n, round(obs_p*m*n))) = 1;
    missing = ~mask;
    %% masked noise synethetic data
    X_masked = gray_fig_d .* mask;
    masked_fig(idx_obs_p,:,:) = X_masked;
    
    % SVT
    tau = (m*n)^args.base.SVT_power;
    step = 1.2 * obs_p;
    SVT_recon = SVT(X_masked, mask, tau, step, max_iter, tol);
    [~, SVT_psnr(idx_obs_p)] = PSNR(gray_fig_d, SVT_recon, missing, 0);
    SVT_rec_fig(idx_obs_p,:,:) = SVT_recon;
    
    for rank_level = args.min_R:args.max_R
        fprintf("rank : %d \n", rank_level);
        % SVP
        step = 1/(1+args.base.delta_2k);
        SVP_recon = SVP(X_masked, mask, step, rank_level, max_iter, tol);
        % Optspace
        tau = 1;
        Opt_recon = optspace(X_masked, mask, rank_level, args.base.opt_tau, max_iter, tol);
        
        [~, tmp_psnr] = PSNR(gray_fig_d, SVP_recon, missing, 0);
        SVP_psnr(idx_obs_p)    = max(tmp_psnr,SVP_psnr(idx_obs_p));
        if SVP_psnr(idx_obs_p) == max(SVP_psnr)
            SVP_rec_fig(idx_obs_p,:,:) = SVP_recon;
        end
        [~, tmp_psnr] = PSNR(gray_fig_d, Opt_recon, missing, 0);
        Opt_psnr(idx_obs_p)    = max(tmp_psnr,Opt_psnr(idx_obs_p));
        if Opt_psnr(idx_obs_p) == max(Opt_psnr)
            Opt_rec_fig(idx_obs_p,:,:) = Opt_recon;
        end
    end
    fprintf("end SVP OPT \n");

    chosen_algo = "admm";
    [X_admm_recon_rank_list, admm_info] = tnnr_recon(X_masked, mask, chosen_algo, args);
    chosen_algo = "apgl";
    [X_apgl_recon_rank_list, apgl_info] = tnnr_recon(X_masked, mask, chosen_algo, args);
    chosen_algo = "admmap";
    [X_admmap_recon_rank_list, admmap_info] = tnnr_recon(X_masked, mask, chosen_algo, args);

    for rank_level = args.min_R:args.max_R
        admm_recon = reshape( X_admm_recon_rank_list(rank_level - args.min_R + 1, :, :, :), m,n,dim);
        [~, tmp_psnr] = PSNR(gray_fig_d, admm_recon, missing, 0);
        admm_psnr(idx_obs_p) = max(tmp_psnr,admm_psnr(idx_obs_p));
        if admm_psnr(idx_obs_p) == tmp_psnr
            admm_rec_fig(idx_obs_p,:,:) = admm_recon;
            admm_iter(idx_obs_p) =  sum(admm_info.iter_list(rank_level,:))/dim;
        end
        %%
        apgl_recon = reshape( X_apgl_recon_rank_list(rank_level - args.min_R + 1, :, :, :), m,n,dim);
        [~, tmp_psnr] = PSNR(gray_fig_d, apgl_recon, missing, 0);
        apgl_psnr(idx_obs_p) = max(tmp_psnr,apgl_psnr(idx_obs_p));
        if apgl_psnr(idx_obs_p) == tmp_psnr
            apgl_rec_fig(idx_obs_p,:,:) = apgl_recon;
            apgl_iter(idx_obs_p) =  sum(apgl_info.iter_list(rank_level,:))/dim;
        end
        
        admmap_recon = reshape( X_admmap_recon_rank_list(rank_level - args.min_R + 1, :, :, :), m,n,dim);
        [~, tmp_psnr] = PSNR(gray_fig_d, admmap_recon, missing, 0);
        admmap_psnr(idx_obs_p) = max(tmp_psnr,admmap_psnr(idx_obs_p));
        if admmap_psnr(idx_obs_p) == tmp_psnr
            admmap_rec_fig(idx_obs_p,:,:) = admmap_recon;
            admmap_iter(idx_obs_p) =  sum(admmap_info.iter_list(rank_level,:))/dim;
        end        
    end
    
end

%% plot Fig4
figure;
subplot(331)
imshow(gray_fig, [], 'border', 'tight')
xlabel('(a) Original image')
set(gcf,'color','none'); set(gca,'color','none');
subplot(332)
imshow(reshape(masked_fig(obs_p_fig,:,:),m,n), [], 'border', 'tight')
xlabel('(b) Masked image')
set(gcf,'color','none'); set(gca,'color','none');
subplot(333);
imshow(clip(reshape(SVT_rec_fig(obs_p_fig,:,:),m,n),0,255) ./ 255, 'border', 'tight');    % show the recovered image
% imshow(clip(reshape(SVT_rec_fig(obs_p_fig,:,:),m,n),0,1) , [], 'border', 'tight');    % show the recovered image
xlabel('(c) SVT');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(334);
imshow(clip(reshape(SVP_rec_fig(obs_p_fig,:,:),m,n),0,255) ./ 255, 'border', 'tight');    % show the recovered image
% imshow(clip(reshape(SVP_rec_fig(obs_p_fig,:,:),m,n),0,1) , [], 'border', 'tight');    % show the recovered image
xlabel('(d) SVP');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(335);
imshow(clip(reshape(Opt_rec_fig(obs_p_fig,:,:),m,n),0,255) ./ 255, 'border', 'tight');    % show the recovered image

% imshow(clip(reshape(Opt_rec_fig(obs_p_fig,:,:),m,n),0,1) , [], 'border', 'tight');    % show the recovered image
xlabel('(e) Opt');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(336);
imshow(clip(reshape(admm_rec_fig(obs_p_fig,:,:),m,n),0,255) ./ 255, 'border', 'tight');    % show the recovered image
% imshow(clip(reshape(admm_rec_fig(obs_p_fig,:,:),m,n),0,1) , [], 'border', 'tight');    % show the recovered image
xlabel('(f) tnnr-admm');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(337);
imshow(clip(reshape(apgl_rec_fig(obs_p_fig,:,:),m,n),0,255) ./ 255, 'border', 'tight');    % show the recovered image
% imshow(clip(reshape(apgl_rec_fig(obs_p_fig,:,:),m,n),0,1) , 'border', 'tight');    % show the recovered image
xlabel('(g) tnnr-apgl');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(338);
imshow(clip(reshape(admmap_rec_fig(obs_p_fig,:,:),m,n),0,255) ./ 255, 'border', 'tight');    % show the recovered image
% imshow(clip(reshape(admmap_rec_fig(obs_p_fig,:,:),m,n),0,1) , 'border', 'tight');    % show the recovered image
xlabel('(h) tnnr-admmap');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
%% save Fig4
if ~exist(args.base.save_dir4, 'dir'), mkdir(args.base.save_dir4); end
save_path = [args.base.save_dir4,'/',args.base.fig_name4];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
close all
%% Plot fig5
figure;
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
if ~exist(args.base.save_dir5, 'dir'), mkdir(args.base.save_dir5); end
save_path = [args.base.save_dir5,'/',args.base.fig_name5];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
close all
%% Plot fig6
figure;
obs_p = obs_p_list(idx_obs_p) * 100;
plot(obs_p_list, admm_iter, '-^k',...
     obs_p_list, apgl_iter, '->b',...
     obs_p_list, admmap_iter, '-+r');
xlim([min(obs_p_list),max(obs_p_list)]);
% set(gcf,'color','none'); % background color -> none
% set(gca,'color','none'); % axis color -> none
xlabel('Observed Ratio');
ylabel('Total Iterations');
legend('admm','apgl','admmap','Location','northwest');
legend('boxoff');
%% saving the results
if ~exist(args.base.save_dir6, 'dir'), mkdir(args.base.save_dir6); end
save_path = [args.base.save_dir6,'/',args.base.fig_name6];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
% close all;