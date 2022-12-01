function out_mat = clip(data, min_val, max_val)
    tmp_min = min_val;
    min_val = min(min_val, max_val);
    max_val = max(tmp_min, max_val);
    out_mat = max(min_val,min(data,max_val));
end