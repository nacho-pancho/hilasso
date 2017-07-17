function X=mdlsMassageData(X,params)    
    if nargin < 2
        params = struct();
        params.minv = 0;
        params.maxv = 1;
        params.normalize = 1;
        params.remove_dc = 1;
    end
    %
    % no arguments: return default parameters
    %
    if nargin < 1        
        X=params;
        return;
    end
    
    [M,N] = size(X);
    minX = min(X(:));
    maxX = max(X(:));
    rangeX = maxX - minX;
    if (params.minv ~= minX)
        X = X - (minX-params.minv);
    end
    if (params.maxv-params.minv) ~= rangeX
        X = X * ((params.maxv-params.minv)/rangeX);
    end
    if params.remove_dc
        X = X - repmat(mean(X),M,1);
    end
    if params.normalize
        aux = sqrt(sum(X.*X));
        aux2 = find(aux==0);
        if ~isempty(aux2)
            aux(aux2) = 1;
        end
        X = X.*(repmat(1./aux,M,1));
    end
end
