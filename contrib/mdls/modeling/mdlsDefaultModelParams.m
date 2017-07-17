function params = mdlsDefaultModelParams(stats)
    if exist('stats','var')
        if iscell(stats)
            params = cell(size(stats));
            for c=1:length(stats)
                params{c} = mdlsDefaultModelParams(stats{c});
            end
            return
        else
            params = stats;
        end
    else
        params = struct();
    end
    params.reg_mode = 2; % minimization of the Lagrangian
    params.lambda = 0.1;
    params.theta = 30.0;
    params.kappa = 2.5;
    params.beta = 0.05;
    params.reg_type = 'l1';
    params.lambda_min = 10;
    params.lambda_max = 150;
    params.lla_iter = 5;
    params.L = 0; % maximum number of nonzero elements per reconstructed vector
    params.positive = false;
    params.project = false;
    params.l2err = 8^2*1.15^2*(1/255)^2; % for average 1 pixel level on
                                         % an 8x8 patch
end
