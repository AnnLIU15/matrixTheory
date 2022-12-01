%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% config loading
config = ReadYaml('../configs/fig4.yaml');
if ~exist(config.base.save_dir, 'dir'), mkdir(config.base.save_dir); end
rgb_fig = imread(config.base.img_path);
gray_fig = rgb2gray(rgb_fig);
[m,n] = size(gray_fig);
mask = zeros(m,n);
mask(randperm(m*n, round(config.base.obs_p*m*n))) = 1;
gray_fig_d = double(gray_fig);
[result, X_rec] = admmap_pic(gray_fig_d,mask,config);
masked_fig = gray_fig_d .* mask;

X_best_rec = result.X_best_rec;

figure;
subplot(131)
imshow(gray_fig, [], 'border', 'tight')
title('(a) Original image')
set(gcf,'color','none'); set(gca,'color','none');
subplot(132)
imshow(masked_fig, [], 'border', 'tight')
title('(b) Masked image')
set(gcf,'color','none'); set(gca,'color','none');
subplot(133);
X_best_rec = clip(X_best_rec, 0, 255);
imshow(X_best_rec ./ 255, 'border', 'tight');    % show the recovered image
xlabel('(c) recovered image');
set(gcf,'color','none'); set(gca,'color','none');
% if ~exist(config.base.save_dir, 'dir'), mkdir(config.base.save_dir); end
%% initial the output mat


%% saving the results

% save_path = [config.base.save_dir,'/',config.base.fig_name];
% print(gcf, '-dpdf',[save_path,'.pdf'])
% print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
% close all;