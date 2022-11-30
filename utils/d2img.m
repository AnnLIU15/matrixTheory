function dst_img = d2img(ori_img_d)
    dim_info = size(ori_img_d);
    mat_dimension = length(dim_info);
    dst_img = zeros(dim_info,'uint8');
    if mat_dimension == 2
        [max_val, min_val] = maxmin(ori_img_d);
        tmp_img = (ori_img_d - min_val)/(max_val - min_val)*255;
        dst_img = uint8(tmp_img);
    elseif mat_dimension == 3
        for idx = 1:dim_info(3)
            [max_val, min_val] = maxmin(ori_img_d(:,:,idx));
            tmp_img = (ori_img_d(:,:,idx) - min_val)/(max_val - min_val)*255;
        dst_img(:,:,idx) = uint8(tmp_img);
        end
    else
        error("current version doesn't support dimension = %d",mat_dimension);
    end
end