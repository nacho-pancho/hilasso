function prefix=mdlsGetModelPrefix(params)
    if nargin > 0
        if isfield(params,'D0') && ~isempty(params.D0)
            params.M = size(params.D0,1);
            params.K = size(params.D0,2);
        end
        if ~isfield(params,'mu0')
            params.mu0 = 0;
        end
        if ~isfield(params,'xmu0')
            params.xmu0 = 0;
        end
        if ~isfield(params,'base_name')
            params.base_name = 'unnamed';
        end
        prefix=sprintf('%s/%s-M%04d-k%04d-mu%g-xmu%g',...
                       params.output_dir,...
                       params.base_name,...
                       params.M, ...
                       params.K, ...
                       params.mu0, ...
                       params.xmu0);
    else
        params = struct();
        params.output_dir = '.';
        params.base_name  = 'noname';
        params.M          = 0;
        params.K          = 0;
        params.mu0        = 0;
        params.xmu0       = 0;
        params.lambda     = 0.1;
        prefix=params;
    end
end
