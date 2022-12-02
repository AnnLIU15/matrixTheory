%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% config loading
config = ReadYaml('../configs/fig4.yaml');
rgb_fig = imread(config.base.img_path);
args.min_R = 1;
args.max_R = 10;
args.outer_iter = 100;
args.beta = 1e-3;
args.beta_max = 1e10;
args.kappa = 1e-3;
args.lambda = 0.001;
args.inner_iter = 200;
args.outer_tol = 1e-3;
args.tol = 1e-4;
args.rho0 = 1.1;
gray_fig = rgb2gray(rgb_fig);
[m,n] = size(gray_fig);
mask = zeros(m,n);
mask(randperm(m*n, round(config.base.obs_p*m*n))) = 1;
gray_fig_d = double(gray_fig);
masked_fig = gray_fig_d .* mask;
chosen_algo = "admmap";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
%%
for idx = 1: args.max_R - args.min_R + 1
    figure;
    subplot(131)
    imshow(gray_fig, [], 'border', 'tight')
    xlabel('(a) Original image')
    set(gcf,'color','none'); set(gca,'color','none');
    subplot(132)
    imshow(masked_fig, [], 'border', 'tight')
    xlabel('(b) Masked image')
    set(gcf,'color','none'); set(gca,'color','none');
    subplot(133);
    X_best_rec = clip(X_l_recon_rank_list(idx,:,:,:), 0, 255);
    X_best_rec = round(X_best_rec);
    if size(X_best_rec,1) == 1
        X_best_rec = reshape(X_best_rec,m,n);
    end
    imshow(X_best_rec ./ 255, 'border', 'tight');    % show the recovered image
    xlabel('(c) recovered image');  % below the img
    set(gcf,'color','none'); set(gca,'color','none');
end
% %% saving the results
% if ~exist(config.base.save_dir, 'dir'), mkdir(config.base.save_dir); end
% save_path = [config.base.save_dir,'/',config.base.fig_name];
% print(gcf, '-dpdf',[save_path,'.pdf'])
% print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
% close all;