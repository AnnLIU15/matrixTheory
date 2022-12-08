%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));
% run `addpath(genpath(cd));` firstly if nedded

%% Load origin figure
ori_img = imread('../assets/mountain.jpeg');

% setting save dir
ori_img_d = double(ori_img)/255;
fullpath = mfilename('fullpath'); 
[path,name]=fileparts(fullpath);
path = strcat(path, "/output/");
fprintf("The scripts will output the result to %s,\n we will del margin(inkscape) and mv it to assets",path);
dot_list = 1:100;
figure;
ha = tight_subplot(2,3,[.01 .01],[0.001 0.001],[.001 .001]) ;
axes(ha(1)); 
imshow(ori_img_d,[])
title('Origin','FontSize',6)
rank_list = [0, 1, 3, 5, 12, 15];
[m,n,dim_max] = size(ori_img_d);
tmp_image = zeros(m,n,dim_max);
for idx = 2:6
    %% get singular value of each dim
    axes(ha(idx));
    for dim = 1:dim_max
        [U_r, S_r, V_r] = svds(ori_img_d(:,:,dim),rank_list(idx));
        tmp_image(:,:,dim) = U_r * S_r * V_r';
    end
    imshow(clip(tmp_image,0,1),[])
    title(['rank=',num2str(rank_list(idx))],'FontSize',6)
end
set(gcf,'color','none'); % background color -> none
set(gca,'color','none'); % axis color -> none

%% saving the results
if ~exist(path, 'dir'), mkdir(path); end
name = strcat(path,name);
print(gcf, '-depsc',strcat(name,'.eps'))
print(gcf,'-dsvg',strcat(name,'.svg') )
close all;  % we will del margin and let it to assets

color_list = {'red','green','blue'};
figure;
ha = tight_subplot(2,4,[.1 .1],[0.1 0.1],[.1 .1]) ;
obs_p = 0.5;
[m,n,dim] = size(ori_img_d);
mask_noise = zeros(m, n);
mask_noise(randperm(m*n, round(obs_p*m*n))) = 1;

mask_str = imread("../assets/string.png");
mask_str = ~mask_str;
[row, col] = size(mask_str);
if row ~= m || col~=n
    tmp_mask = mask_str;
    mask_str = ones(m,n);
    begin = round(([m,n] - [row, col])/2);
    mask_str(begin(1)+1:begin(1)+row, begin(2)+1:begin(2)+col) = tmp_mask;
end
mask_block = imread("../assets/windows-mask.png");
mask_block = (mask_block(:,:,1) ~= 0);
[row, col] = size(mask_block);
if row ~= m || col~=n
    tmp_mask = mask_block;
    mask_block = ones(m,n);
    begin = round(([m,n] - [row, col])/2);
    mask_block(begin(1)+1:begin(1)+row, begin(2)+1:begin(2)+col) = tmp_mask;
end
ori_img_d_noise = ori_img_d .* mask_noise;
ori_img_d_str   = ori_img_d .* mask_str;
ori_img_d_block = ori_img_d .* mask_block;
axes(ha(1)); imshow(ori_img_d,[]); title('Origin','FontSize',6)
axes(ha(2)); imshow(ori_img_d_noise, []);title('random mask','FontSize',6)
axes(ha(3)); imshow(ori_img_d_str  , []);title('string mask','FontSize',6)
axes(ha(4)); imshow(ori_img_d_block, []);title('block mask','FontSize',6)
for idx = 1:3
    %% get singular value of each dim
    axes(ha(idx + 4))
    svd_sigma = svds(ori_img_d(:,:,idx),100);
    svd_sigma_noise = svds(ori_img_d_noise(:,:,idx),100);
    svd_sigma_str = svds(ori_img_d_str(:,:,idx),100);
    svd_sigma_block = svds(ori_img_d_block(:,:,idx),100);


    plot(dot_list, svd_sigma,       '-sc', ...
         dot_list, svd_sigma_noise, '-om', ...
         dot_list, svd_sigma_str,   '-dg', ...
         dot_list, svd_sigma_block, '-^k','MarkerSize',2);
    title(['(',char(96+idx),') ',color_list{idx},' channel']);
    legend('origin','noise mask','string mask','block mask');
    legend('boxoff');
    xlabel('the i-th singluar value'); ylabel('Magnitude');
    set(gcf,'color','none'); % background color -> none
    set(gca,'color','none','xtick',0:10:100); % axis color -> none
end

