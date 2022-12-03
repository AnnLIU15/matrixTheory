%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% args loading
args = ReadYaml('../configs/fig4v2.yaml');
rgb_fig = imread(args.base.img_path);

gray_fig = rgb2gray(rgb_fig);
[m,n] = size(gray_fig);
mask = zeros(m,n);
mask(randperm(m*n, round(args.base.obs_p*m*n))) = 1;
gray_fig_d = double(gray_fig);
masked_fig = gray_fig_d .* mask;
chosen_algo = "admm";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
if size(X_l_recon_rank_list,1) == 1
    admm_recon = reshape(X_l_recon_rank_list,m,n);
end
chosen_algo = "apgl";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
if size(X_l_recon_rank_list,1) == 1
    apgl_recon = reshape(X_l_recon_rank_list,m,n);
end
chosen_algo = "admmap";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
if size(X_l_recon_rank_list,1) == 1
    admmap_recon = reshape(X_l_recon_rank_list,m,n);
end
%%


max_iter = args.base.max_iter; tol = args.base.tol;r0 = args.min_R;
tau = (m*n)^args.base.SVT_power;
step = 1.2 * args.base.obs_p;
SVT_recon = SVT(masked_fig, mask, tau, step, max_iter, tol);
step = 1/args.base.obs_p/sqrt(max_iter);
SVP_recon = SVP(masked_fig,mask,step,r0,max_iter,tol);
Opt_recon = OPTSPACE(masked_fig, r0);
%% plot
figure;
subplot(331)
imshow(gray_fig, [], 'border', 'tight')
xlabel('(a) Original image')
set(gcf,'color','none'); set(gca,'color','none');
subplot(332)
imshow(masked_fig, [], 'border', 'tight')
xlabel('(b) Masked image')
set(gcf,'color','none'); set(gca,'color','none');
subplot(333);
imshow(clip(SVT_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
xlabel('(c) SVT');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(334);
imshow(clip(SVP_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
xlabel('(d) SVP');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(335);
imshow(clip(Opt_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
xlabel('(e) Opt');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(336);
imshow(clip(admm_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
xlabel('(f) tnnr-admm');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(337);
imshow(clip(apgl_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
xlabel('(g) tnnr-apgl');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
subplot(338);
imshow(clip(admmap_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
xlabel('(h) tnnr-admmap');  % below the img
set(gcf,'color','none'); set(gca,'color','none');
% save
if ~exist(args.base.save_dir, 'dir'), mkdir(args.base.save_dir); end
save_path = [args.base.save_dir,'/',args.base.fig_name];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
% close all;