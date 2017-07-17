%
% Function: scMOL2 (using mexWeightedLasso)
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
function A=scMOL2(X,D,L,lambda,kappa,betta,mode,iters,tol)
  if nargin < 5
    error('Must specify at least X,D,L,lambda,beta');
  end
  params=struct();
  N = size(X,2);
  K = size(D,2);
  M = size(X,1);
  maxL = min(K,M);

  if (L > 0) && (L < maxL)
      params.L = L;
  else
      params.L = maxL;
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

  %
  % first iteration: weighted with lambda' = lambda/beta
  %
  if params.mode == 2
      % the value of lambda passed here should be 2*sigma^2*(kappa+1)
      % but on the first iteration we want to use a 2*sigma^*theta where
      % theta is the mean value of the Laplacian parameter under the MOL
      % distribution, which is kappa/beta, so we approximate it by
      % 2*sigma^2*(kappa+1)/lambda, which actually corresponds to theta
      % being the MODE of the Gamma distribution.
      %
      params.lambda = lambda*kappa/(kappa+1)/betta;
  else
      params.lambda = lambda;
  end
  A = mexLasso(X, D, params);  
  %
  % subsequent with lambda / (betta + |a|)
  %
  params.lambda = lambda;
  % process in chunks to avoid out of memory
  N2 = min(1000,N);
  for i = 1:iters
      for j = 1:N2:(N-N2+1)
          idx = j:(j+N2-1);
          W = betta ./ ( betta + abs(A(:,idx)) );
          A(:,idx) = mexLassoWeighted(X(:,idx), D, W, params);
      end
  end
end
