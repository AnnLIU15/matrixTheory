%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));
m = 100;
n = 200;
R = 10;
max_iter = 1000;
tol = 1e-3;

para.outer_iter = 10;     % maximum number of iteration
para.outer_tol = 1e-3;     % tolerance of iteration

para.admm_iter = 30;    % iteration of the ADMM optimization
para.admm_tol = 1e-4;    % epsilon of the ADMM optimization
para.admm_rho = 5e-2;    % rho of the ADMM optimization     beta

para.apgl_iter = 30;    % iteration of the APGL optimization
para.apgl_tol = 1e-4;    % epsilon of the APGL optimization
para.apgl_lambda = 5e-2; % lambda of the APGL optimization

para.admmap_iter = 30;  % iteration of the ADMMAP optimization
para.admmap_tol = 1e-4;  % epsilon of the ADMMAP optimization
para.admmap_rho = 2;     % rho0 of the ADMMAP optimization
para.admmap_kappa = 1e-3;% kappa of the ADMMAP optimization
para.admmap_beta = 1; % beta0 of the ADMMAP optimization
para.admmap_betamax = 1e10; % max beta of the ADMMAP optimization
%% different noise levels and different observed ratios
samp_rate = 0.9;
for ii = 1:10
    s = 0.1 * ii;
    % Generating matrix
    M = randn(m, R) * randn(R, n);
    Z = randn(m, n);
    % Generate synthetic data
    X = M + s * Z;
    mask = zeros(m, n);
    mask(randperm(m*n, round(samp_rate*m*n))) = 1;
    missing = ones(size(mask)) - mask; 
    X = X .* mask;
    % SVT
    tau = sqrt(m*n); 
    step = 1.2 * samp_rate; 
    SVT_recon = SVT(X, mask, tau, step, max_iter, tol);

    % SVP
    step = 1/samp_rate/sqrt(max_iter); 
    SVP_recon = SVP(X,mask,step,R,max_iter,tol);

    % Optspace
    Opt_recon = OPTSPACE(X, R);

    % TNNR-admm
    admm_recon = admm_mat(M, X, mask, R, para);
    
    % TNNR-apgl
    [apgl_recon, ~] = apgl_mat(M, X, mask, R, para);

    % TNNR-admmap
    [admmap_recon, ~] = admmap_mat(M, X, mask, R, para);
    
    SVT_erec(ii) =      norm(vec((M-SVT_recon).*missing), 'fro');
    SVP_erec(ii) =      norm(vec((M-SVP_recon).*missing), 'fro');
    Opt_erec(ii) =      norm(vec((M-Opt_recon).*missing), 'fro');
    admm_erec(ii) =     norm(vec((M-admm_recon).*missing), 'fro');
    apgl_erec(ii) =     norm(vec((M-apgl_recon).*missing), 'fro');
    admmap_erec(ii) =   norm(vec((M-admmap_recon).*missing), 'fro');
   
end

% plot and save results
s = 0.1:0.1:1;

plot(s, SVT_erec, '-sc', s, SVP_erec, '-om', s, Opt_erec, '-dg', ...
    s, admm_erec, '-^k', s, apgl_erec, '->b', s, admmap_erec, '-+r');
xlim([0.1 1]);
set(gcf,'color','none'); % background color -> none
set(gca,'color','none'); % axis color -> none
xlabel('Noise level');
ylabel('Total reconstruction error');
legend('SVT', 'SVP', 'Optspace', 'admm', 'apgl', 'admmap','Location','northwest');
legend('boxoff');

% saving the results
save_dir = '../output/fig2/';
% setting save dir
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
print(gcf, '-dpdf',[save_dir,'/fig2_',num2str(samp_rate),'.pdf'])

%% different observed ratios on matrices with different ranks
clear SVT_ecer SVP_erec Opt_erec admm_erec apgl_erec admmap_erec
s = 0.5;
samp_rate = 0.6;
for R = 1:15  % rank 
    % Generating matrix
    M = randn(m, R) * randn(R, n);
    Z = randn(m, n);
    % Generate synthetic data
    X = M + s * Z;
    mask = zeros(m, n);
    mask(randperm(m*n, round(samp_rate*m*n))) = 1;
    missing = ones(size(mask)) - mask; 
    X = X .*mask;
    % SVT
    tau = sqrt(m*n); 
    step = 1.2 * samp_rate; 
    SVT_recon = SVT(X, mask, tau, step, max_iter, tol);

    % SVP
    step = 1/samp_rate/sqrt(max_iter); 
    SVP_recon = SVP(X,mask,step,R,max_iter,tol);

    % Optspace
    Opt_recon = OPTSPACE(X, R);

    % TNNR-admm
    admm_recon = admm_mat(M, X, mask, R, para);
    
    % TNNR-apgl
    [apgl_recon, ~] = apgl_mat(M, X, mask, R, para);

    % TNNR-admmap
    [admmap_recon, ~] = admmap_mat(M, X, mask, R, para);
    
    SVT_erec(R) = norm(vec((M-SVT_recon).*missing), 'fro');
    SVP_erec(R) = norm(vec((M-SVP_recon).*missing), 'fro');
    Opt_erec(R) = norm(vec((M-Opt_recon).*missing), 'fro');
    admm_erec(R) = norm(vec((M-admm_recon).*missing), 'fro');
    apgl_erec(R) = norm(vec((M-apgl_recon).*missing), 'fro');
    admmap_erec(R) = norm(vec((M-admmap_recon).*missing), 'fro');
    
%     SVT_erec(R) = norm((M-SVT_recon).*missing, 'fro');
%     SVP_erec(R) = norm((M-SVP_recon).*missing, 'fro');
%     Opt_erec(R) = norm((M-Opt_recon).*missing, 'fro');
%     admm_erec(R) = norm((M-admm_recon).*missing, 'fro');
%     apgl_erec(R) = norm((M-apgl_recon).*missing, 'fro'); 
%     admmap_erec(R) = norm((M-admm_recon).*missing, 'fro');
end

%% plot and save result
R = 1:15;

plot(R, SVT_erec(1:15), '-sc', R, SVP_erec(1:15), '-om', R, Opt_erec(1:15), '-dg', R, admm_erec(1:15), '-^k', R, apgl_erec(1:15), '->b', R, admmap_erec(1:15), '-+r');
xlim([1 15]);
xlabel('Rank');
set(gcf,'color','none'); % background color -> none
set(gca,'color','none','xtick',1:15); % axis color -> none
ylabel('Total reconstruction error');
legend('SVT', 'SVP', 'Optspace', 'admm', 'apgl', 'admmap','Location','northwest');
legend('boxoff');

% saving the results
save_dir = '../output/fig3/';
% setting save dir
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
print(gcf, '-dpdf',[save_dir,'/fig3_',num2str(samp_rate),'.pdf'])
