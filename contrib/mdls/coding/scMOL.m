%
% Function: scMOL (using mexLasso + change of variables)
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
function A=scMOL(X,D,L,lambda,kappa,betta,mode,iters,tol)
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
      params.lambda = lambda / betta;
  else
      params.lambda = lambda;
  end
  A = mexLasso(X, D, params);  

  params.lambda = lambda;
  %
  % subsequent with w = lambda / (betta + |a|)
  %
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
  % since w = 1/(|a|+beta), we have 1/w = (|a|+beta)
  %
  matlabpool
  parfor (j=1:N)
      for i=1:iters
          Wj = abs(A(:,j)) + betta;
          WD = spdiags(Wj, 0, K, K);
          tmpA = mexLasso(X(:,j), D*WD, params);
          A(:,j) = tmpA .* Wj;
      end % iters
  end % samles
  matlabpool close
end
