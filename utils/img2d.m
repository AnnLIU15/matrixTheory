function ori_img_d = img2d(ori_img)
    dim_info = size(ori_img);
    mat_dimension = length(dim_info);
    ori_img_d = zeros(dim_info);
    ori_img = double(ori_img);
    if mat_dimension == 2
        %% Normalize the one channel image(uint8)
        [max_val, min_val] = maxmin(ori_img);
        tmp_img = (ori_img - min_val)/(max_val - min_val+1e-10);
        ori_img_d = double(tmp_img);
    elseif mat_dimension == 3
        %% Normalize R, G, and B channels, respectively
        for idx = 1:dim_info(3)
            [max_val, min_val] = maxmin(ori_img(:,:,idx));
            tmp_img = (ori_img(:,:,idx) - min_val)/(max_val - min_val+1e-10);
            ori_img_d(:,:,idx) = double(tmp_img);
        end
    else
        error("current version doesn't support dimension = %d",mat_dimension);
    end
end