function syn_mat = getSynthetic(config_base,sigma,begin)
    noise_mat = sigma * randn(config_base.m, config_base.n);
    r0,ratio_p
    mat_Ml = randn(config_base.m,config_base.r0);   % m x r0
    mat_Mr = randn(config_base.r0,config_base.n);   % r0 x n
    mat_M = mat_Ml * mat_Mr; % rank r0
    mat_Z = randn(config_base.m,config_base.n);
    mat_B_ori = mat_M + sigma * mat_Z;
    mask = zeros(config_base.m, config_base.n);
    chosen = randperm(config_base.m, config_base.n,...
        round(obs_ratio_p*config_base.m, config_base.n));
    mask(chosen) = 1 ;
end