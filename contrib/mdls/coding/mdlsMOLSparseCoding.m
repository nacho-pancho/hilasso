%
% function A = mdlsMOLSparseCoding(X,D,A0,params)
%
% Name: mdlsMOLSparseCoding
%
% Category: core function
%
% Description: sparse coding based on MOL prior penalty
%
% Input:
% X ........ data (n x N)
% D ........ dictionary (n x K)
% A0 ....... initial coefficients (K x N). If a scalar of value a is passed
%            instead, A0 will be computed using the penalized least squares
%            solution with penalty a.
% params ... see PDSparseModeling for details
%
% Output:
% A ........ updated reconstruction coefficients
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
function A = mdlsMOLSparseCoding(X,D,params,A0)
    format short
    K=size(D,2);
    N=size(X,2);
    n = size(X,1);
    if  nargin < 4
         % if no initial A is given, start with  LS estimates
        A0 = mdlsLS(X,D,1e-3);
    end
    sigma2 = params.varError;
    kappa = params.kappaCoef;
    betta = params.betaCoef;
    aux=2.^(0:K-1);
    Dhat = D;
    A = A0;
    converged = 0;
    lambda = 2*sigma2*(kappa+1);
    fprintf('L=%d\tlambda=%8.6f\tbeta=%8.6f\titers=%d\ttol=%8.6f\n',...
            params.sparsity,lambda,betta,params.codingIters,params.codingTolerance);
    for j=1:N
        fprintf('j=%d... ',j);
        Aj = A(:,j);
        Xj = X(:,j);
        i = 1;
        change = realmax;
        while i <= params.codingIters
            Wj = (betta + abs(Aj));
            % reweighted Lasso
            Aprev = Aj;
            % Aj2 = scWeightedLasso(Xj,D,n,Wj,2);
            % using Lasso after change of variables
            % if this gives the same as Aj, then Aj is OK, but
            % otherwise what I've been doing is WRONG!
            % SHIT, it seems to be the case!!
            for k=1:K
                Dhat(:,k) = D(:,k)*Wj(k);
            end
            if params.debug
                %fprintf('iter=%d\nWj=',i-1);  Wj'
                %fprintf('DtD=\n');
                %Dhat'*Dhat
                %fprintf('Rdn=\n');
                %(Dhat'*Xj)'
            end            
            Aj = scLasso(Xj,Dhat,n,lambda,2).*Wj;
            %fprintf('Ahat=\n'); Aj
            % Aj = Aj.*Wj;
            %fprintf('A=\n'); Aj
            %fprintf('Aprev=\n'); Aprev
            %fprintf('Difference between Aj2 and Aj  %f\n',max(full(abs(Aj2 - Aj))));
            change = full(sum(abs(Aj-Aprev)));
            %fprintf('change=%f\n',change);
            if (change < params.codingTolerance)
                converged = converged + 1;
                break
            end
            i = i + 1;
        end
        A(:,j)=Aj;
    end
    fprintf('Converged in %d out of %d samples\n',converged,N);
end