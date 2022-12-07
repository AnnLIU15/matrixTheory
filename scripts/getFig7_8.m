%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));
%% args loading
yaml_str = '../configs/fig7_8.yaml';
args = ReadYaml(yaml_str);
datetime.setDefaultFormats('default',"yyyy/MM/dd HH:mm:ss:SSS")
for idx = 1:2
    %% fig loading
    fprintf("Start Fig. %d: %s\n",idx+6,datestr(datetime('now')));
    rgb_fig = imread(args.base.img_path{idx});
    mask = imread(args.base.mask_path);
    mask = ~mask;
    rgb_fig = double(rgb_fig);
    [m,n,dim] = size(rgb_fig);
    [row, col] = size(mask);
    if row ~= m || col~=n
        tmp_mask = mask;
        mask = ones(size(rgb_fig(:,:,1)));
        begin = round(([m,n] - [row, col])/2);
        mask(begin(1)+1:begin(1)+row, begin(2)+1:begin(2)+col) = tmp_mask;
    end
    missing = ~mask;
    masked_fig = rgb_fig .* mask;
    best_admm_psnr = 0; best_apgl_psnr = 0; best_admmap_psnr = 0;
    best_svt_psnr = 0; best_svp_psnr = 0; best_opt_psnr = 0;
    best_admm = zeros(m,n,dim); best_apgl = zeros(m,n,dim); best_admmap = zeros(m,n,dim); best_svt =  zeros(m,n,dim); best_svp =  zeros(m,n,dim); best_opt =  zeros(m,n,dim);
    best_r = zeros(6,1);
    %%  ADMM
    cur_time = datetime('now');
    fprintf("\tStart TNNR-ADMM: %s ->",datestr(cur_time));chosen_algo = "admm";
    X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
    fprintf(" %.3fs\n",getMSecDiff(cur_time));
    for rank_idx = 1:size(X_l_recon_rank_list,1)
        admm_recon = reshape(X_l_recon_rank_list(rank_idx,:,:,:),m,n,dim);
        [~, admm_psnr] = PSNR(rgb_fig, admm_recon, missing, 0);
        if best_admm_psnr < admm_psnr
            best_admm_psnr = admm_psnr;
            best_admm = admm_recon;
            best_r(4) = rank_idx;
        end
    end
    
    %% APGL
    cur_time = datetime('now');
    fprintf("\tStart TNNR-AGPL: %s ->",datestr(cur_time));chosen_algo = "apgl"; % best in 10
    X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
    fprintf(" %.3fs\n",getMSecDiff(cur_time));
    for rank_idx = 1:size(X_l_recon_rank_list,1)
        apgl_recon = reshape(X_l_recon_rank_list(rank_idx,:,:,:),m,n,dim);
        [~, apgl_psnr] = PSNR(rgb_fig, admm_recon, missing, 0);
        if best_apgl_psnr < apgl_psnr
            best_apgl_psnr = apgl_psnr;
            best_apgl = apgl_recon;
            best_r(5) = rank_idx;
        end
    end

    %% ADMMAP
    cur_time = datetime('now');
    fprintf("\tStart TNNR-ADMMAP: %s ->",datestr(cur_time));chosen_algo = "admmap";
    X_l_recon_rank_list = tnnr_recon(masked_fig, mask, chosen_algo, args);
    fprintf(" %.3fs\n",getMSecDiff(cur_time));

    for rank_idx = 1:size(X_l_recon_rank_list,1)
        admmap_recon = reshape(X_l_recon_rank_list(rank_idx,:,:,:),m,n,dim);
        [~, admmap_psnr] = PSNR(rgb_fig, admmap_recon, missing, 0);
        if best_admmap_psnr < admmap_psnr
            best_admmap_psnr = admmap_psnr;
            best_admmap = admmap_recon;
            best_r(6) = rank_idx;
        end
    end

    %% baseline config
    obs_p = sum(mask,'all')/m/n;
    max_iter = args.base.max_iter; tol = args.base.tol;
    %% SVT
    tau = (m*n)^args.base.SVT_power;
    step = 1.2 * obs_p;
    cur_time = datetime('now');
    fprintf("\tStart SVT: %s ->",datestr(cur_time));
    SVT_recon_R = SVT(masked_fig(:,:,1), mask, tau, step, max_iter, tol);
    SVT_recon_G = SVT(masked_fig(:,:,2), mask, tau, step, max_iter, tol);
    SVT_recon_B = SVT(masked_fig(:,:,3), mask, tau, step, max_iter, tol);
    best_svt = cat(3,SVT_recon_R,SVT_recon_G,SVT_recon_B);
    fprintf(" %.3fs\n",getMSecDiff(cur_time));
    [~, best_svt_psnr] = PSNR(rgb_fig, best_svt, missing, 0);

    %% SVP AND OptSpace
    cur_time = datetime('now');
    fprintf("\tStart SVP and OptSpace: %s ->",datestr(cur_time));% best rank_idx = 20;
    for rank_idx = args.min_R:args.max_R
        step = 1/(1+args.base.delta_2k); %\delta_{2k}<1/3
        SVP_recon_R = SVP(masked_fig(:,:,1),mask,step,rank_idx,max_iter,tol);
        SVP_recon_G = SVP(masked_fig(:,:,2),mask,step,rank_idx,max_iter,tol);
        SVP_recon_B = SVP(masked_fig(:,:,3),mask,step,rank_idx,max_iter,tol);
        SVP_recon = cat(3,SVP_recon_R,SVP_recon_G,SVP_recon_B);
        [~, SVP_psnr] = PSNR(rgb_fig, SVP_recon, missing, 0);
        if best_svp_psnr < SVP_psnr
            best_svp_psnr = SVP_psnr;
            best_svp = SVP_recon;
            best_r(2) = rank_idx;
        end
        Opt_recon_R = optspace(masked_fig(:,:,1), mask, rank_idx, args.base.opt_tau,  max_iter, tol);
        Opt_recon_G = optspace(masked_fig(:,:,2), mask, rank_idx, args.base.opt_tau,  max_iter, tol);
        Opt_recon_B = optspace(masked_fig(:,:,3), mask, rank_idx, args.base.opt_tau,  max_iter, tol);
        Opt_recon = cat(3,Opt_recon_R,Opt_recon_G,Opt_recon_B);
        [~, Opt_psnr] = PSNR(rgb_fig, Opt_recon, missing, 0);
        if best_opt_psnr < Opt_psnr
            best_opt_psnr = Opt_psnr;
            best_opt = Opt_recon;
            best_r(3) = rank_idx;
        end
    end
    fprintf(" %.3fs\n",getMSecDiff(cur_time));
    %% plot
    figure;
    ha = tight_subplot(2,4,[.01 .01],[0.001 0.001],[.001 .001]) ;
    disp(best_r')
    axes(ha(1)); 
    imshow(rgb_fig./ 255, [], 'border', 'tight')
    title('(a) Original image','FontSize',6)
    % set(gcf,'color','none'); set(gca,'color','none');
    axes(ha(2)); 
    imshow(masked_fig./ 255, [], 'border', 'tight')
    title('(b) Masked image','FontSize',6)
    % set(gcf,'color','none'); set(gca,'color','none');
    axes(ha(3)); 
    imshow(clip(best_svt,0, 255)./255,[], 'border', 'tight');    % show the recovered image
    title(['(c) SVT PSNR=',sprintf("%.2f",best_svt_psnr)],'FontSize',6);  % below the img
    % set(gcf,'color','none'); set(gca,'color','none');
    axes(ha(4)); 
    imshow(clip(best_svp,0, 255)./255,[], 'border', 'tight');    % show the recovered image
    title(['(d) SVP PSNR=',sprintf("%.2f",best_svt_psnr)],'FontSize',6);  % below the img
    % set(gcf,'color','none'); set(gca,'color','none');
    axes(ha(5)); 
    imshow(clip(best_opt,0, 255)./255,[], 'border', 'tight');    % show the recovered image
    title(['(e) OptSpace PSNR=',sprintf("%.2f",best_opt_psnr)],'FontSize',6);  % below the img
    % set(gcf,'color','none'); set(gca,'color','none');
    axes(ha(6)); 
    imshow(clip(best_admm,0, 255)./255,[], 'border', 'tight');    % show the recovered image
    title(['(f) tnnr-admm PSNR=',sprintf("%.2f",best_admm_psnr)],'FontSize',6);  % below the img
    % set(gcf,'color','none'); set(gca,'color','none');
    axes(ha(7)); 
    imshow(clip(best_apgl,0, 255)./255,[], 'border', 'tight');    % show the recovered image
    title(['(g) tnnr-apgl PSNR=',sprintf("%.2f",best_apgl_psnr)],'FontSize',6);  % below the img
    % set(gcf,'color','none'); set(gca,'color','none');
    axes(ha(8)); 
    imshow(clip(best_admmap,0, 255)./255,[], 'border', 'tight');    % show the recovered image
    title(['(h) tnnr-admmap PSNR=',sprintf("%.2f",best_admmap_psnr)],'FontSize',6);  % below the img
    set(gcf,'color','none'); set(gca,'color','none');
    if ~exist(args.base.save_dir{idx}, 'dir'), mkdir(args.base.save_dir{idx}); end
    save_path = [args.base.save_dir{idx},'/',args.base.fig_name{idx}];
    print(gcf, '-dpdf',[save_path,'.pdf'])
    print(gcf,[save_path,'.jpeg'] ,'-djpeg','-r300')
end
% close all;


