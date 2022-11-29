%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% loading the configuration of experiment 
config = ReadYaml('configs/ADMM.yaml');
if ("ADMM" ~= config.algo{1}) || (size(config.algo,2) ~= 1) 
    error("not satified the requirement, algo name=%s, algo size=%d",...
        config.algo{1},size(config.algo,2));
else
    save_dir = ['./output/',config.algo{1},'/'];
    % setting save dir
    if ~exist(save_dir, 'dir'),mkdir(save_dir); end
    tol = config.vars.tol;                         % stopping criteria
    maxIter = config.vars.maxIter;                 % maximum allowable iterations
    mask_samp_rate = config.vars.mask_samp_rate;
end

%% Load origin figure
ori_img = imread('assets/RGB_figure.jpg');
ori_img_d = double(ori_img);

%% Normalize R, G, and B channels, respectively
img_R = ori_img_d(:,:,1); img_G = ori_img_d(:,:,2); img_B = ori_img_d(:,:,3);
min_R = min(img_R(:)); Io_R = img_R-min_R; img_R = Io_R/max(Io_R(:));
min_G = min(img_G(:)); Io_G = img_G-min_G; img_G = Io_G/max(Io_G(:)); 
min_B = min(img_B(:)); Io_B = img_B-min_B; img_B = Io_B/max(Io_B(:));

%% Random mask
[nx,ny,~] = size(ori_img_d);
mask = zeros(nx,ny);
chosen = randperm(nx*ny,round(mask_samp_rate*nx*ny));
mask(chosen) = 1 ;

%% Masked image
mask_R = img_R.*mask; mask_G = img_G.*mask; mask_B = img_B.*mask;
mask_image = cat(3,mask_R,mask_G,mask_B);

%% TNNR-ADMM
beta = 1; rank_r = 1;
TNNR_recon_R = TNNR_ADMM(mask_R,mask,beta,rank_r,maxIter,tol);
TNNR_recon_G = TNNR_ADMM(mask_G,mask,beta,rank_r,maxIter,tol);
TNNR_recon_B = TNNR_ADMM(mask_B,mask,beta,rank_r,maxIter,tol);
TNNR_image = cat(3,TNNR_recon_R,TNNR_recon_G,TNNR_recon_B);


%% Experimental results
figure; imshow(ori_img,[]); title('Original image','FontSize',15,'FontName','Times New Roman'); 
figure; imshow(mask_image,[]); title(['Masked image (sampling rate = ', num2str(mask_samp_rate),')'],'FontSize',15,'FontName','Times New Roman'); 
figure; imshow(TNNR_image, 'border', 'tight','initialmagnification','fit');

%% save the Experimental results
set(gca,'position',[0 0 1 1])
fig_pos = [config.algo{1},'_',num2str(mask_samp_rate),...
    '_',num2str(beta),'_',num2str(rank_r),'_'...
    ,num2str(maxIter),'_',num2str(tol)];
print(gcf, '-dpdf',[save_dir,fig_pos,'.pdf'])
print(gcf,[save_dir,fig_pos,'.jpeg'] ,'-djpeg','-r300')
close all;