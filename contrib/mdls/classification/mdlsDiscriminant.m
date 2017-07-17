% function F = mdlsDiscriminant(D,X,mode,lambda,max_nz)
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
% max_nz ........ Maximum number of nonzero elements in reconstruction. 
%               By default sparsity is implicitly controlled by lambda,
%               but this can be fixed to a maximum for speed purposes.
% true_label .... If specified, show error rate and confusion matrix.
% 
% coherence_thres .... If in (0,1), this will ignore, when computing the
%                      discriminant function F, the coefficients a_k such
%                      that the corresponding atom k in the dictionary is
%                      too coherent (angle greater than the threshold) to
%                      some atom from other class.
%
% Output:
% label .......... (1xN) The assigned label for each sample
%
function F = mdlsDiscriminant(D,X,mode,lambda,ct)
    if ~exist('params','var')
        params = struct();
        params.mode = mode;
        params.lambda = lambda;
        params.ct = ct;
    end

    if ~iscell(D)
        error('D must be a cell of size 1xC.');
    end
    if ~exist('params.ct','var')
        params.ct = 0;
    end
    %
    % ct >= 1 is trivial
    %
    if abs(params.ct) >= 1
        params.ct = 0;
    end

    if params.ct == 0
        disc_mode = 0;
    elseif params.ct > 0
        disc_mode = 1;
    else % params.ct < 0
        disc_mode = 2;
        params.ct = -params.ct;
    end

    C = length(D);
    Nt = size(X,2);
    F = single(zeros(C,Nt));
    if params.ct > 0
        %
        % build concatenation of rival dictionaries
        %
        Dr = cell(1,C);
        for c=1:C
            Dr{c} = [];
            for c2=[(1:c-1) (c+1:C)]
                Dr{c} = [Dr{c} D{c2}];
            end
        end
        %
        % build black list of atoms that are too coherent 
        %
        blacklist = cell(1,C);
        whitelist = cell(1,C); 
       for c=1:C
            aux = Dr{c}'*D{c};
            % aux (1xKc) = max coherence of each atom in D{c} with
            % an atom from any other dictionary.
            aux = max(aux); 
            Kc = size(D{c},2);
            blacklist{c} = find(aux >= params.ct);
            whitelist{c} = setdiff(1:Kc,blacklist{c});
            fprintf('Class %d has %d blacklisted atoms\n',c,length(blacklist{c}));
        end        
    end
    %
    % main loop
    %
    for c=1:C
        Dc = D{c};        
        if disc_mode == 2
            Dc = Dc(:,whitelist{c});
        end
        max_nz = min(size(X,1)-1, size(Dc,2));
        if isequal('reg',params.mode)
            switch (disc_mode) 
              case {0,2}
                Ac = scLasso(X, Dc, max_nz, params.lambda, 1);
                F(c,:) = full( sum(abs( Ac )) );
              case 1
                Ac = scLasso(X, Dc, max_nz, params.lambda, 1);
                F(c,:) = full( sum(abs( Ac(whitelist{c},:) )) );
            end
            clear Ac;

        elseif isequal('fit',params.mode)
            Ac =[];
            switch (disc_mode)
              case 0
                Ac = scLasso(X, Dc, max_nz, params.lambda, 0);
              case 2
                Ac = scLasso(X, Dc(:,whitelist{c}), max_nz, params.lambda, 0);
            end
            Ec = single(X - Dc*Ac);
            F(c,:) = sum(Ec.^2); 
            clear Ec Ac;

        elseif isequal('lasso',params.mode)
            Ac = [];
            switch(disc_mode)
              case {0,2}
                Ac = scLasso(X, Dc, max_nz, params.lambda, 2);
                F(c,:) = full( sum(abs( Ac )) );
              case 1
                Ac = scLasso(X, Dc, max_nz, params.lambda, 2);
                F(c,:) = full( sum(abs( Ac(whitelist{c},:) )) );
            end
            Ec = single(X - double(Dc)*Ac); 
            F(c,:) = params.lambda*F(c,:) + 0.5*sum(Ec.^2); 
            clear Ac Ec;
        end
    end
end