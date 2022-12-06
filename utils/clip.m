function out_mat = clip(data, min_val, max_val)
%--------------------------------------------------------------------------
%  clip data to range (min_val, max_val)
%--------------------------------------------------------------------------
%     compute PSNR and reconstruction error for the recovered image and
%     original image
%
%     Inputs:
%         data           --- original data
%         min_val        --- min value
%         max_val        --- max value
%
%     Outputs:
%         out_mat        --- clip data
%--------------------------------------------------------------------------
    tmp_min = min_val;
    min_val = min(min_val, max_val);
    max_val = max(tmp_min, max_val);
    out_mat = max(min_val,min(data,max_val));
end