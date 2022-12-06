function [ mse, psnr ] = PSNR( X_full, X_rec, missing, type)
%--------------------------------------------------------------------------
%  10 \times \log10(\frac{255^2}{MSE}) normalization ->  10 \times \log10(\frac{1}{MSE})
%--------------------------------------------------------------------------
%     compute PSNR and reconstruction error for the recovered image and
%     original image
%
%     Inputs:
%         X_full           --- original image
%         X_rec            --- recovered image
%         missing          --- index matrix of missing elements
%         type             --- 0 -> 255(ori) 1 -> 1(normalization)
%
%     Outputs:
%         mse              --- reconstruction mse error
%         psnr             --- PSNR (Peak Signal-to-Noise Ratio)
%--------------------------------------------------------------------------
if type == 0
    max_val = 255;
elseif type == 1
    max_val = 1;
else
    error("PSNR -- unsupport type = %d",type)
end
X_rec = clip(X_rec, 0, max_val);

[~,~, dim] = size(X_rec);
mse = norm(vec((X_full-X_rec).*missing))^2 / dim /nnz(missing);

% MSE = (total mean squard error)/3T
psnr = 10 * log10(max_val^2 / mse);
% 10 \times \log10(\frac{255^2}{MSE}) normalization ->  10 \times \log10(\frac{1}{MSE})
end