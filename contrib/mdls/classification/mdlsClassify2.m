% function [label,cm]=mdlsClassify(D,X,mode,lambda [,max_nz,true_label],ct )
%
% Category: core
%
% Purpose: Perform classification using a sparseland model based on a set
%          of dictionaries (given in D)
%
% Description:
%
%
% INPUT
%
% X ....... matrix of data samples ordered as coloumns.
% D ....... dictionary of atoms.
% params .. model parameters. This is a combination of the standard model
%           parameters (see help mdlsDefaultModelParams) with some
%           specific ones:
%              ct ........ coherence threshold: atoms whose inner
%                          product with an atom from some other
%                          dictionary is higher than this value 
%                          will not be considered in the
%                          reconstruction.
%              disc_mode ..discriminative function mode. Can be
%              'reg','fit' or 'sum'.
% true_label .... If specified, show error rate and confusion matrix.
% verbose . numeric verbosity level. 0 means nothing.
%
% OUTPUT
%
% label .......... (1xN) The assigned label for each sample
% cm ............. Confusion matrix, if labels were given
%
function [est_label,cm]=mdlsClassify2(X,D,true_label,model_params,ct,verbose)
    if ~exist('model_params','var')
        model_params = mdlsDefaultModelModel_Params(); % lambda,kappa,beta and all
                                           % that
        params.disc_mode = 'sum';
    end
    if ~exist('ct','var')
        ct = 0;
    end
    if ~exist('verbose','var')
        verbose = 0;
    end
    if nargout > 1
        cm =[];
    end
    if ~iscell(D)
        error('D must be a cell of size 1xC.');
    end
    class_labels = sort(unique(true_label));
    C = length(class_labels);
    Nt = size(X,2);
    F = mdlsDiscriminant2(X,D,model_params,verbose);
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
