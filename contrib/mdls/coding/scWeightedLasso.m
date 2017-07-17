% function A = scWeightedLasso(X,D,L,lambda,W,mode,iters,tol)
%
% Purpose: solve the weighted-Lasso sparse model
% ||X-DA|| + s^2*log phi(A) where phi(A) is a weighted Lasso with
% difeferent weights per col.
%
%
%
% Inputs:
%  X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%  D ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%  L ....... Sparsity level: indicates maximum number of nonzero entries in the representation of each column of X.
%            This means that each column of A must have no more than s nonzero entries.
%  lambda .. Penalty term, equivalent to sigma^2*(kappa+1)
%  W ....... Kx1 : weight corresponding to each row, usually the ML
%            estimate of the Laplacian parameter for each row. NOTE
%            DIFFERENCE WITH mexWeightedLasso, where W is for each ELEMENT.
%  mode .... Selects minimization problem to be solved:
%
%   0) min_A ||X-D*A||_2 s.t. phi(A_i) <= lambda
%   1) min_A phi(A) s.t. ||X_i - D*A_i|| <= lambda
%   2) min_A ||x-D*A|| + lambda*phi(A)
%
%  iters .... Maximum number of weighted LASSO iterations
%  tol ...... Minimum change in solution required in order to proceed with iterations.
%
% Outputs:
%  A ... Decomposition coefficients, argument of minimization.
%
function A = scWeightedLasso(X,D,L,lambda,W,mode)
  if nargin < 5
    error('Must specify at least X,D,L,lambda, beta');
  end
  params = struct();
  M = size(X, 1);
  if (L > 0) && (L <= M)
      params.L = L;
  else
      params.L = M;
  end
  if ~exist('mode','var')
      params.mode = 1;
  else
      params.mode = mode;
  end
  N = size(X, 2);
  K = size(D, 2);
  params.lambda = lambda;

  %
  % to solve ||x-D*a||_2^2 + ||w.*a||_1
  %
  % we do the change of variables
  %
  % a' = w.*a -> a = a'./w
  % 
  % and solve
  % 
  % ||x-D'a'||_2^2 + ||a'||_1
  %
  % with D' = D diag(1./w)
  % 
  % and then revert a = a'./w
  %
  W = 1 ./ (W+eps);

  if size(W,2) == 1
      WD = spdiags(W, 0, K, K);
  elseif size(W,2) ~= N
      error('W must be KxN or Kx1');
  end
  A = spalloc(K, N, params.L);
  matlabpool
  parfor (j=1:N)
      if size(W,2) == N
          WD = spdiags(W(:,j), 0, K, K);
      end
      tmpA = mexLasso(X(:,j), D*WD, params);
      A(:,j) = tmpA .* W(:,j);
  end
  matlabpool close
end
