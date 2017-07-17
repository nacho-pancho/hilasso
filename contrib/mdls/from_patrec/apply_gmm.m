%
% Perform classification using a Gaussian Mixture Model. Uses netlab gmm* routines.
%
% INPUT
% X ............ (nxp) n d-dimensional training samples 
% model ........ (1x1) already computed GMM mode.
% t ............ (nx1) true labels.
%
% OUTPUT
% e ............ (1x1) classification error
% cm ...........       confusion matrix
%
function gt = apply_gmm(X,model,t)
    P = gmmpost(model,X)';
    [dum,gt] = max(P);
    cm = mdlsConfusionMatrix(gt,t);
    gt = mdlsFixUnsupervisedLabels(unique(t), cm, gt);
end

