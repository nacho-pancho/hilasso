% function [label,cm]=mdlsClassifySoft(D,X,mode,lambda [,max_nz,true_label] )
%
% Category: core
%
% Purpose: Perform soft classification using a sparseland model based on a set
%          of dictionaries (given in D). Given the representation energies
%          obtained with each dictionary for a given patch, {E_1,...,E_c}
%          a probability is assigned to each class i=1...c as p(i) = e^{-E_i}/sum_j(e^{-E_j})
%
% Description:
%
% Input:
% D  .......... cell of C (MxK)  dictionaries.  M is the dimension of the feature 
%               vectors, K is the number of atoms in each dictionary (same) and
%                C is the number of classes.
%
% X ........... (MxN) N M-dimensional feature vectors.
% mode ........ Each sample Xj is assigned to the class whose dictionary produces
%               the minimum cost function value. The cost function is selected by
%               this variable. The valid values are:
%               'fit' ..... Cost(c) = ||X-Dc*Ac||
%               'reg' ..... Cost(c) = ||Ac||_1 
%               'lasso' ... Cost(c) = ||X-Dc*Ac||+lambda*||Ac||_1
%
% lambda ...... Depending on mode, this can be the maximum allowed reconstruction
%               error (for mode 'reg'), the maximum allowed ||A||_1 
%               (for mode 'fit') or the lasso penalty.
%
% true_label .... If specified, show error rate and confusion matrix.
%
% Output:
% labels .......... (1xN) The soft labels assitned to each sample
% cm ............. Soft confusion matrix, if labels were given
%
function [est_label,cm]=mdlsSoftClassify(D, X, mode, lambda, true_label)
    if nargout > 1
        cm =[];
    end
    if ~iscell(D)
        error('D must be a cell of size 1xC.');
    end
    class_labels = sort(unique(true_label));
    C = length(class_labels);
    Nt = size(X,2);
    F = mdlsDiscriminant(D,X,mode,lambda);
    est_label = exp(-F);
    est_label = est_label./repmat(sum(est_label),C,1);
    if (nargout > 1) && exist('true_label','var') && ~isempty(true_label)
        cm = mdlsSoftConfusionMatrix(est_label,true_label);
    end
end