% function F = mdlsDiscriminant2(D,X,params)
%
% Category: core
%
% Purpose: Perform classification using a sparseland model based on a set
%          of dictionaries (given in D)
%
% Description:
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
% verbose . numeric verbosity level. 0 means nothing.
%
% OUTPUT
%
% label .......... (1xN) The assigned label for each sample
%
function F = mdlsDiscriminant2(X,D,model_params,ct,verbose)

    if ~exist('model_params','var')
        model_params = mdlsDefaultModelParams();
        ct = 0;
    end
    if ~exist('verbose','var')
        verbose = 0;
    end

    if ~iscell(D)
        error('D must be a cell of size 1xC.');
    end
    if ~exist('ct','var')
        ct = 0;
    end
    %
    % ct >= 1 is trivial
    %
    if abs(ct) >= 1
        ct = 0;
    end

    if ct == 0
        atom_filter_mode = 0;
    elseif ct > 0
        atom_filter_mode = 1;
    else % ct < 0
        atom_filter_mode = 2;
        ct = -ct;
    end

    C = length(D);
    Nt = size(X,2);
    F = single(zeros(C,Nt));
    if ct > 0
        [blacklist,whitelist] = mdlsBlackList(D,ct);
    end
    %
    % main loop
    %
    for c=1:C
        Dc = D{c};        
        if iscell(model_params)
            model_paramsc = model_params{c};
        else
            model_paramsc = model_params;
        end
        if atom_filter_mode == 2
            Dc = Dc(:,whitelist{c});
        end
        max_nz = min(size(X,1)-1, size(Dc,2));
        Ac = mdlsSparseCoding(X, Dc, model_paramsc);
        if atom_filter_mode == 1
            Ac(blacklist{c},:) = 0;
        end
        F(c,:) = mdlsModelEnergy(X,Dc,Ac,model_paramsc);
    end
end