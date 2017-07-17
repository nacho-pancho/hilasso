%
% Learn a Gaussian Mixture Model. Uses netlab gmm* routines.
%
% INPUT
% X ............ (nxp) n d-dimensional training samples 
% nc ........... (1x1) how many clusters to learn
% 
% OUTPUT
% model ........ trained GMM. 
%
function model = learn_gmm(X, nc, covar_model)
    fprintf('Gaussian Mixture Model (diag covar)\n');
    if ~exist('covar_model','var')
        covar_model = 'diag';
    end
    dim = size(X,2);
    %
    % create structure
    %
    fprintf('GMM: creating structure...\n');
    model = gmm(dim, nc, covar_model);
    options = foptions;
    options(14) = 10;
    %
    % initialize GMM
    %
    fprintf('GMM: initializing...\n');
    model = gmminit(model,X,options);
    %
    % train GMM
    %
    fprintf('GMM: training...\n');
    [model,opt2,errlog] = gmmem(model,X,options);
    %
    % compute posterior probabilities
    %
end
