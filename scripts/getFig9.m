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
datetime("now")
chosen_algo = "admm";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);

if size(X_l_recon_rank_list,1) == 1
    admm_recon = reshape(X_l_recon_rank_list,m,n,dim);
end
imshow(clip(admm_recon,0,255) ./ 255, 'border', 'tight');
datetime("now")
chosen_algo = "apgl";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
if size(X_l_recon_rank_list,1) == 1
    apgl_recon = reshape(X_l_recon_rank_list,m,n,dim);
end
chosen_algo = "admmap";
X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
if size(X_l_recon_rank_list,1) == 1
    admmap_recon = reshape(X_l_recon_rank_list,m,n,dim);
end
%% baseline
max_iter = args.base.max_iter; tol = args.base.tol;r0 = args.min_R;
tau = (m*n)^args.base.SVT_power;
step = 1.2 * args.base.obs_p;
SVT_recon_R = SVT(masked_fig(:,:,1), mask, tau, step, max_iter, tol);
SVT_recon_G = SVT(masked_fig(:,:,2), mask, tau, step, max_iter, tol);
SVT_recon_B = SVT(masked_fig(:,:,3), mask, tau, step, max_iter, tol);
SVT_recon = cat(3,SVT_recon_R,SVT_recon_G,SVT_recon_B);
step = 1/args.base.obs_p/sqrt(max_iter);
SVP_recon_R = SVP(masked_fig(:,:,1),mask,step,r0,max_iter,tol);
SVP_recon_G = SVP(masked_fig(:,:,2),mask,step,r0,max_iter,tol);
SVP_recon_B = SVP(masked_fig(:,:,3),mask,step,r0,max_iter,tol);
SVP_recon = cat(3,SVP_recon_R,SVP_recon_G,SVP_recon_B);
Opt_recon_R = OPTSPACE(masked_fig(:,:,1), r0);
Opt_recon_G = OPTSPACE(masked_fig(:,:,2), r0);
Opt_recon_B = OPTSPACE(masked_fig(:,:,3), r0);
Opt_recon = cat(3,Opt_recon_R,Opt_recon_G,Opt_recon_B);
%% plot
figure;
ha = tight_subplot(4,2,[.035 .3],[0.1 0.1],[.01 .01]) ;

axes(ha(1)); 
imshow(rgb_fig./255, [], 'border', 'tight')
title('(a) Original image')
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(2)); 
imshow(masked_fig./255, [], 'border', 'tight')
title('(b) Masked image')
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(3)); 
imshow(clip(SVT_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
title('(c) SVT');  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(4)); 
imshow(clip(SVP_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
title('(d) SVP');  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(5)); 
imshow(clip(Opt_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
title('(e) Opt');  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(6)); 
imshow(clip(admm_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
title('(f) tnnr-admm');  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(7)); 
imshow(clip(apgl_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
title('(g) tnnr-apgl');  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
axes(ha(8)); 
imshow(clip(admmap_recon,0,255) ./ 255, 'border', 'tight');    % show the recovered image
title('(h) tnnr-admmap');  % below the img
% set(gcf,'color','none'); set(gca,'color','none');
% save
if ~exist(args.base.save_dir, 'dir'), mkdir(args.base.save_dir); end
save_path = [args.base.save_dir,'/',args.base.fig_name];
print(gcf, '-dpdf',[save_path,'.pdf'])
print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
% close all;