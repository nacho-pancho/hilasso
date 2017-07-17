%
% Function: scJOE (using mexWeightedLasso)
%
% Purpose: solve the unconstrained 'Mixture Of Laplacians' sparse model
% ||X-DA|| + s^2*log phi(A) where phi(A) is the JOE prior.
%
% The technique used is the LLA approximation, which consists of solving
% a sequence of convex linear approximations to the problem, each of which boils
% down to a weighted lasso where the weights are determined by phi'(A) evaluated
% at the current iterate A_0.
% 
% For the IID Joe model we have phi'(x) = (-log p(x))' = 
%       [-log 1/x - log ( ae^{-bx} - be^{-bx} )/( e^{-ax} - e^{-bx} )]' =
%             1/x + ( ae^{-bx} - be^{-bx} )/( e^{-ax} - e^{-bx} ) 
%
% whose limit when x -> 0+ is (b^2-a^2)/(b-a). Thus we can define at in
% x=0 safely using this value.
%
% the first iterate is computed using the unweighted Lasso and a penalty
% that corresponds to the average value of the Laplacian parameter under
% the constrained Jeffreys prior (which is 1/ln(b/a)*1/x in the interval
% [a,b]), that is, (b-a)/ln(b/a).
%
% in the MAP model, this corresponds to the initial lambda being 2*sigma^2*(b-a)/ln(b/a)
%
% Usage: [A]=scJOE(X,D,s,lambda,beta,mode[,iters,tol])
%
%
% Inputs:
%  X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%  D ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%  L ....... Sparsity level: indicates maximum number of nonzero entries in the representation of each column of X.
%            This means that each column of A must have no more than s nonzero entries.
%  lambda .. Penalty term, equivalent to sigma^2*(kappa+1)
%  a ....... lower end of interval for JOE prior.
%  b ....... upper end of interval for JOE prior.
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
function A = scJOE(X,D,L,lambda,a,b,mode,iters,tol)
  if nargin < 5
    error('Must specify at least X, D, L, lambda, a, b');
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

  N = size(X,2);
  K = size(D,2);
  M = size(X,1);

  params.lambda = lambda;
  %
  % first iteration: unweighted 
  %
  if a > b
      error('a must be less than b');
  elseif a==b
      warning('a==b means this is ordinary L1 regularization (using lambda=lambda*a)');
      params.lambda = lambda * a;
      A  = mexLasso(X, D, params);
      return;
  else
      C0 = (b^2-a^2) / (b-a);
  end
  params.lambda = lambda * C0;
  A  = mexLasso(X, D, params);
  matlabpool
  parfor (j=1:N)
      for i = 1:iters
          w = C0*ones(K,1);
          nzi = find(A(:,j));
          Anz = A(:,j);
          eaAj      = exp(-a*Anz);
          ebAj      = exp(-b*Anz);
          tmp       = abs(eaAj - ebAj);
          w(nzi)  = 1./(1./Anz + (a*eaAj - b*ebAj) ./ (abs(eaAj - ebAj) + ...
                                                       eps));
          % change of variables and LASSO
          WD = spdiags(w, 0, K, K);
          tmpA = mexLasso(X(:,j), D*WD, params);
          A(:,j) = tmpA .* w;          
      end % iters
  end % samples
  matlabpool close
end
