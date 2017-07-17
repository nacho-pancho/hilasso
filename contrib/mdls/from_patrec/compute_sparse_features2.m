%
% D ............. (1xc cell) c dictionaries
% X ............. (mxn) data samples as columns
% mode .......... sparse coding mode: 'reg','fit','lagrangian'
% lambda ........ regularization parameter (role depends on mode)
% feature_mode .. this is a string which tells how to build the
%                 features out of the sparse coding result.
%                 It has the form xxx+xxx. The first xxx applies
%                 to how the reconstruction residual is added to the
%                 feature vector, and the second xxx to the
%                 rec. coefficients. Both can have one of three values:
%                 xxx = vec : include full vector as part of the feature
%                 xxx = sum : include a scalar which is the (abs) sum of the
%                             elements.  
%                 xxx = nul : exclude as a feature
%                       VEC : same as vec but preserving the signs.
%
% Example: for "vec+sum", the resulting (m+1)-dimensional  
%                         feature vector is (note the squared elements)
%                         (e_1^2 e_2^2 ... e_m^2 ||a||_1)
%              "sum+sum", the resulting 2-dimensional feature
%                         vector is (||e||^2_2 + ||a||_1)
%
% OUTPUT
% F ............ (pxn) n vectors with features of dimension p
% w0 ........... (cxp) linear weights corresponding to the d
%                       generative-based discirminant function:
%                       when each w is multiplied by F, we get the 
%                       reconstruction energy of each sample.
%
function  [F,w0] = compute_sparse_features2(D,X,mode,lambda,feature_mode, ...
                                      verbose)
    if ~exist('verbose','var')
        verbose = 0;
    end
    %
    %  dimension of feature vectors
    %
    p = 0;
    n = size(X,2);
    m = size(X,1);
    C = length(D);
    pc = zeros(1,C);
    apc = zeros(1,C);
    epc = zeros(1,C);

    
    vlambda = lambda;
    e_mode = feature_mode(1:3);
    a_mode = feature_mode(5:7);
    for c=1:C
        switch e_mode
          case {'vec','VEC'}
            ep = m;
          case 'sum'
            ep = 1;
          case 'tot'
            ep = 1;
          case 'nul'
            ep = 0;
        end
        epc(c) = ep;
        switch a_mode
          case {'vec','VEC'}
            ap = size(D{c},2);
          case 'sum'
            ap = 1;
          case 'nul'
            ap = 0;
          otherwise
            ap = 0;
        end
        apc(c) = ap;
        pc(c) = ep + ap;
        p = p + pc(c);
    end
    %
    % create w0
    %
    w0 = [];
    for jj=1:length(vlambda)        
        lambda = vlambda(jj);
        w0t = zeros(C,p);
        for c=1:C
            ci     = 1+ sum( pc(1:(c-1)) );
            cf     = ci + pc(c) - 1;
            w0t(c,ci:cf) = [ -0.5*ones(1,epc(c)) -lambda*ones(1,apc(c)) ];
        end
        w0 = [w0 w0t];
    end
    if verbose > 0
        fprintf(['Will generate %d feature vectors of dimension %d using ' ...
        '%d dictionaries\n'],n,size(w0,2),C);
    end

    max_nz = m-1;
    %
    % compute energies
    %
    pi = 1;
%    F  = zeros(p,n); % this one can be huge...
    for jj=1:length(vlambda)
        
        lambda = vlambda(jj);
        for c=1:C
            if verbose
                fprintf('Computing features for dictionary %d\n',c);
            end

            switch mode
                case 'reg',
                    A = scLasso(X,D{c},max_nz,lambda,1);
                case 'fit',
                    A = scLasso(X,D{c},max_nz,lambda,0);
                case 'lagrangian',
                    A = scLasso(X,D{c},max_nz,lambda,2);
            end
            E  = X  - D{c}*A;
            if ~isequal(e_mode,'VEC')
                E  = E.*E;
            end
            if ~isequal(a_mode,'VEC')
                A  = abs(A);
            end

            switch e_mode
              case {'vec','VEC'}
                pf = pi + m - 1;
                F(pi:pf,:) = E;
              case 'sum'
                pf = pi;
                F(pi:pf,:) = sum(E);
              case 'tot'
                pf = pi;
                F(pi:pf,:) = 0.5*sum(E) + lambda*sum(A);
              case 'nul'
                pf = pi - 1;
            end
            pi = pf + 1;

            switch a_mode
              case {'vec','VEC'}
                pf = pi + size(D{c},2) - 1;
                F(pi:pf,:) = A;
              case 'sum'
                pf = pi;
                F(pi:pf,:) = sum(A);
              case 'nul'
                pf = pi - 1;
              otherwise
                pf = pi - 1;
            end
            pi = pf + 1;
        end  % for each dictionary
    end
end