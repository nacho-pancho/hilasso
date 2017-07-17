% function [label,cm]=mdlsClassify(D,X,mode,lambda [,max_nz,true_label],ct )
%
% Category: core
%
% Purpose: Perform classification using a sparseland model based on a set
%          of dictionaries (given in D)
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
% coherence_thres .... If in (0,1), this will ignore, when computing the
%                      discriminant function F, the coefficients a_k such
%                      that the corresponding atom k in the dictionary is
%                      too coherent (angle greater than the threshold) to
%                      some atom from other class.
%
% Output:
% label .......... (1xN) The assigned label for each sample
% cm ............. Confusion matrix, if labels were given
%
function [est_label,cm]=mdlsClassify(D,X,mode,lambda, true_label,ct)
    if nargout > 1
        cm =[];
    end
    if ~iscell(D)
        error('D must be a cell of size 1xC.');
    end
    if ~exist('ct','var')
        ct=0;
    end
    class_labels = sort(unique(true_label));
    C = length(class_labels);
    Nt = size(X,2);
    F = mdlsDiscriminant(D,X,mode,lambda,ct);
    [fm,est_label] = min(F);
    for t=1:Nt
        est_label(t) = class_labels(est_label(t));
    end
    if size(true_label,1) > 1
        true_label = true_label';
    end
    if (nargout > 1) && exist('true_label','var') && ~isempty(true_label)
        cm = mdlsConfusionMatrix(est_label,true_label);
    end
end
