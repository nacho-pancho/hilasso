%
% Function: scMOL (using mexWeightedLasso)
%
% Purpose: solve the unconstrained 'Mixture Of Laplacians' sparse model ||X-DA|| + s^2*(kappa+1) sum(log(|A_i|+beta)
% 
% Usage: [A]=scMOL(X,D,s,lambda,beta,mode[,iters,tol])
%
% Inputs:
%  X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%  D ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%  L ....... Sparsity level: indicates maximum number of nonzero entries in the representation of each column of X.
%            This means that each column of A must have no more than s nonzero entries.
%  lambda .. Penalty term, equivalent to sigma^2*(kappa+1)
%  beta .... MOL scale parameter.
%  mode .... Selects minimization problem to be solved:
%   0) min_A ||X-D*A||_2 s.t. ||A_i||_1 <= lambda
%   1) min_A ||A||_1 s.t. ||X_i - D*A_i|| <= lambda
%   2) min_A ||x-D*A|| + lambda*||A||
%
%  iters .... Maximum number of weighted LASSO iterations
%  tol ...... Minimum change in solution required in order to proceed with iterations.
%
% Outputs:
%  A ... Decomposition coefficients, argument of minimization.
%
function A=scMOL(X,D,L,lambda,betta,mode,iters,tol)
  if nargin < 5
    error('Must specify at least X,D,L,lambda,beta');
  end
  params=struct();

  if (L > 0) && (L < size(X,1))
      params.L = L;
  end
  if ~exist('mode','var')
      params.mode = 1;
  else
      params.mode = mode;
  end
  if ~exist('iters','var')
      iters = 1; % one-step estimation a la Zou
  end
  if ~exist('tol','var')
      tol = size(X,1)*1e-4; 
  end

  params.lambda = lambda;
  
  N = size(X,2);
  K = size(D,2);
  if size(betta,1) == 1
      betta = repmat(betta,K,N);
  else
      betta = repmat(betta,1,N);
  end
  %
  % first iteration: weighted with lambda' = lambda/beta
  %
  W = 1./betta;
  A = mexLassoWeighted(X, D, W, params);  
  %
  % subsequent with lambda / (betta + |a|)
  %
  params.lambda = lambda;
  for i = 1:iters
      W = 1 ./ ( betta + abs(A) );
      A = mexLassoWeighted(X, D, W, params);
  end
end
