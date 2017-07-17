%
% Function: scRowLasso
%
% Purpose: solve a weighted Lasso problem where each row of A is given a
% different lambda.
% 
% Usage: [A]=scRowLasso(X,D,L,T,W,mode)
%
% Inputs:
%  X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%  D ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%  L ....... Sparsity level: indicates maximum number of nonzero entries in the representation of each column of X.
%            This means that each column of A must have no more than s nonzero entries.
%  T ....... Penalty term. Role depends on the value of mode (See below).
%  W ....... Weights for each row
%  mode .... Selects minimization problem to be solved:
%
%   0) min_A ||X-D*A||_2^2 s.t. ||W.*A_i||_1 <= T
%   1) min_A ||W.*A||_1 s.t. ||X_i - D*A_i|| <= T
%   2) min_A ||x-D*A||_2^2 + T*||W.*A||_1
%
%
% Outputs:
%  A ... Decomposition coefficients, argument of minimization.
%
function A=scRowLasso(X,D,L,T,W,mode)
  if nargin < 5
    error('Must specify at least X,D,L,T,W');
  end
  if ~exist('mode','var')
      mode = 1;
  end
  params=struct();
  N = size(X,2);
  K = size(D,2);
  M = size(X,1);
  maxL = min(K,M);

  if (L > 0) && (L < maxL)
      L = L;
  else
      L = maxL;
  end
  lambda = T;
  WD = spdiags(1./W, 0, K, K);
  A = scLasso(X, D*WD, L, lambda, mode);
  A = WD*A;
end
