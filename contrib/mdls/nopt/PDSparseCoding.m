%
% function A = PDSparseCoding(X,D,A0,params)
%
% Name: PDSparseCoding
%
% Category: core function
%
% Description: sparse coding for the PD problem.
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
function A = PDSparseCoding(X,D,A,params)
  %
  % the variance of the error is the parameter of the fitting term which
  % to simplify the problem gets factorized out, which means it multiplies
  % the value of Lambda
  %
  sigma2 = params.varError;
  %
  % this is the value of lambda equivalent to a first run of MOL
  %
  b = params.betaCoef;
  k = params.kappaCoef;
  if params.lambda == 0
    lambda = sigma2 * k / b;
  else
    lambda = params.lambda;
  end
  %fprintf('Equivalent L1 LAMBDA=%f\n',lambda);
  method = params.codingMethod;

  sparsity = params.codingSparsity;
  N = size(X,2);
  K = size(D,2);
  Ns = ceil(N/10);
  nd = 0;
  if isequal(method,'IS')
    %
    % -------- ITERATED SHRINKAGE
    %
    c = 0.1*norm(D'*D);
    %c = 1;
    %fprintf('IS constant c=%f\n',c);
    A = LeastSquaresSparseCoding(X,D,1e-1,2);
    A = IteratedShrinkage(X,D,A,c,lambda,params.codingIters,params.codingTolerance);
  elseif isequal(method,'Lasso')
    %
    % -------- LASSO
    %
    A = mexLasso(X,D, sparsity,lambda,2);
  else % end LASSO
    %
    % -------- MOL
    %
    if ~params.codingIncremental
      if isequal(params.codingSubMethod,'IS')
        c = 10*norm(D'*D);
        A = LeastSquaresSparseCoding(X,D,lambda,2);
        %A = IterativeShrinkage(X,D,A,c,lambda,params.codingIters,params.codingTolerance);
        A = IteratedShrinkage(X,D,A,c,lambda,params.codingIters,params.codingTolerance); % C implementation
      elseif isequal(params.codingSubMethod,'Lasso')
        A = mexLasso(X,D, sparsity,lambda,2);
      else
        error('Unknown coding submethod.');
      end
      i = 1;
      change = realmax;
      while (i <= params.codingIters) && (change > params.codingTolerance)
        A0 = A;
        W = sigma2*(k+1)./(b + abs(A));
        if isequal(params.codingSubMethod,'IS')
          for j = 1:N
            if (mod(j,Ns) == 0)
              nd = nd + 1;
              fprintf(1,'.');
            end
            A(:,j) = IterativeShrinkage(X(:,j),D,A(:,j),c,W(:,j),params.codingIters,params.codingTolerance);
          end
        elseif isequal(params.codingSubMethod,'Lasso')
          A = WeightedLasso(X,D,sparsity,W,2);
        else
          error('Unknown coding submethod.');
        end
        change = full(sum(sum(abs(A-A0))));
        if params.debug == 1
          mW = full(max(W(:)));
          avgAA = full(sum(sum(abs(A))))/N/K;
          fprintf(1,'||A||=%f max(W)=%f change: %f\n',avgAA,mW,change);
        end
        i = i + 1;
      end % coding iterations
    else % incremental mode
      if isequal(params.codingSubMethod,'IS')
        c = 1.1*norm(D'*D);
        A = LeastSquaresSparseCoding(X,D,lambda,2);
        A = IterativeShrinkage(X,D,A,c,lambda,params.codingIters,params.codingTolerance);
      elseif isequal(params.codingSubMethod,'Lasso')
        A = mexLasso(X,D, sparsity,lambda,2);
      else
        error('Unknown coding submethod.');
      end
      for j = 1:N % process each sample at a time
        Aj = A(:,j);
        Xj = X(:,j);
        i = 1;
        change = realmax;
        while (i <= params.codingIters) && (change > params.codingTolerance)
          A0 = Aj;
          Wj = 5*sigma2*(k+1)./(b + abs(Aj));
          Aj = WeightedLasso(Xj,D,sparsity,Wj,2);
          change = full(sum(sum(abs(Aj-A0))));
          i = i + 1;
        end % coding iterations
        A(:,j) = Aj;
      end  % for each sample
    end % end incremental mode
  end % end if MOL
  for i = 1:nd
    fprintf('\b');
  end
end % function
