function [ mse, psnr ] = PSNR( X_full, X_rec, missing )
%--------------------------------------------------------------------------
% Xue Shengke, Zhejiang University, April 2017.
% Contact information: see readme.txt.
%
% Hu et al. (2013) TNNR paper, IEEE Transactions on PAMI.
% First written by debingzhang, Zhejiang Universiy, November 2012.
%--------------------------------------------------------------------------
%     compute PSNR and reconstruction error for the recovered image and
%     original image
% 
%     Inputs:
%         X_full           --- original image
%         X_rec            --- recovered image
%         missing          --- index matrix of missing elements
% 
%     Outputs: 
%         erec             --- reconstruction error
%         psnr             --- PSNR (Peak Signal-to-Noise Ratio)
%--------------------------------------------------------------------------

X_rec = clip(X_rec, 0, 255);

[~,~, dim] = size(X_rec);
mse = norm(vec((X_full-X_rec).*missing))^2 / dim /nnz(missing);
% MSE = (total mean squard error)/3T
psnr = 10 * log10(255^2 / mse);
% 10 \times \log10(\frac{255^2}{MSE}) 
end