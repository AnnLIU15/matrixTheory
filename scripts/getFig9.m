%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% args loading
args = ReadYaml('../configs/fig9.yaml');
rgb_fig = imread(args.base.img_path);
mask = imread(args.base.mask_path);
mask = (mask(:,:,1) ~= 0);
rgb_fig = double(rgb_fig);
[m,n,dim] = size(rgb_fig);
masked_fig = rgb_fig .* mask;



%%  ADMM
cur_time = datetime('now');
fprintf("Start TNNR-ADMM: %s ->",datestr(cur_time));chosen_algo = "admm";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
fprintf(" %.3fs\n",getMSecDiff(cur_time));
if size(X_l_recon_rank_list,1) == 1
    admm_recon = reshape(X_l_recon_rank_list,m,n,dim);
end
%% APGL
cur_time = datetime('now');
fprintf("Start TNNR-AGPL: %s ->",datestr(cur_time));chosen_algo = "apgl";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
fprintf(" %.3fs\n",getMSecDiff(cur_time));
if size(X_l_recon_rank_list,1) == 1
    apgl_recon = reshape(X_l_recon_rank_list,m,n,dim);
end
%% ADMMAP
cur_time = datetime('now');
fprintf("Start TNNR-ADMMAP: %s ->",datestr(cur_time));chosen_algo = "admmap";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
fprintf(" %.3fs\n",getMSecDiff(cur_time));
if size(X_l_recon_rank_list,1) == 1
    admmap_recon = reshape(X_l_recon_rank_list,m,n,dim);
end
%% baseline config
obs_p = sum(mask,'all')/m/n;
max_iter = args.base.max_iter; tol = args.base.tol;r0 = args.min_R;
%% SVT
tau = (m*n)^args.base.SVT_power;
step = 1.2 * obs_p;
cur_time = datetime('now');
fprintf("Start SVT: %s ->",datestr(cur_time));
SVT_recon_R = SVT(masked_fig(:,:,1), mask, tau, step, max_iter, tol);
SVT_recon_G = SVT(masked_fig(:,:,2), mask, tau, step, max_iter, tol);
SVT_recon_B = SVT(masked_fig(:,:,3), mask, tau, step, max_iter, tol);
SVT_recon = cat(3,SVT_recon_R,SVT_recon_G,SVT_recon_B);
fprintf(" %.3fs\n",getMSecDiff(cur_time));
%% SVP
cur_time = datetime('now');
fprintf("Start SVP: %s ->",datestr(cur_time));
step = 1/(1+args.base.delta_2k); %\delta_{2k}<1/3
SVP_recon_R = SVP(masked_fig(:,:,1),mask,step,r0,max_iter,tol);
SVP_recon_G = SVP(masked_fig(:,:,2),mask,step,r0,max_iter,tol);
SVP_recon_B = SVP(masked_fig(:,:,3),mask,step,r0,max_iter,tol);
SVP_recon = cat(3,SVP_recon_R,SVP_recon_G,SVP_recon_B);
fprintf(" %.3fs\n",getMSecDiff(cur_time));
%% OptSpace
tau = args.base.opt_tau;
cur_time = datetime('now');
fprintf("Start Optspace: %s ->",datestr(cur_time));
Opt_recon_R = optspace(masked_fig(:,:,1), mask, r0,tau, max_iter, tol);
Opt_recon_G = optspace(masked_fig(:,:,2), mask, r0,tau, max_iter, tol);
Opt_recon_B = optspace(masked_fig(:,:,3), mask, r0,tau, max_iter, tol);
Opt_recon = cat(3,Opt_recon_R,Opt_recon_G,Opt_recon_B);
fprintf(" %.3fs\n",getMSecDiff(cur_time));
%% plot
figure;
ha = tight_subplot(2,4,[.01 .01],[0.001 0.001],[.001 .001]) ;

axes(ha(1)); 
imshow(rgb_fig./255, [], 'border', 'tight')
title('(a) Original image','FontSize',6)
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(2)); 
imshow(masked_fig./255, [], 'border', 'tight')
title('(b) Masked image','FontSize',6)
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(3)); 
imshow(clip(round(SVT_recon),0, 255)./255,[], 'border', 'tight');    % show the recovered image
title('(c) SVT','FontSize',6);  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(4)); 
imshow(clip(round(SVP_recon),0, 255)./255,[], 'border', 'tight');    % show the recovered image
title('(d) SVP','FontSize',6);  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(5)); 
imshow(clip(round(Opt_recon),0, 255)./255,[], 'border', 'tight');    % show the recovered image
title('(e) Opt','FontSize',6);  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(6)); 
imshow(clip(round(admm_recon),0, 255)./255,[], 'border', 'tight');    % show the recovered image
title('(f) tnnr-admm','FontSize',6);  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(7)); 
imshow(clip(round(apgl_recon),0, 255)./255,[], 'border', 'tight');    % show the recovered image
title('(g) tnnr-apgl','FontSize',6);  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(8)); 
imshow(clip(round(admmap_recon),0, 255)./255,[], 'border', 'tight');    % show the recovered image
title('(h) tnnr-admmap','FontSize',6);  % below the img
set(gcf,'color','none'); set(gca,'color','none');
%% save
if ~exist(args.base.save_dir, 'dir'), mkdir(args.base.save_dir); end
save_path = [args.base.save_dir,'/',args.base.fig_name];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
close all;