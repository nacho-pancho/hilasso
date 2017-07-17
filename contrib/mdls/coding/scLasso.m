%
% Wrapper for using Julien's mexLasso as if we were calling the old scLasso
%
% function [A]=scLasso(X,D,s,lambda,mode,ols,pos)
% 
%   Purpose: 
%     Solves various L1-regularized problems (see 'mode') using LARS.
% 
%     See B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, "Least angle regression," Annals of Statistics, 32(2):407--499, 2004.
%  
%   Inputs:
%    X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%    D ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%    s ....... Sparsity level: indicates maximum number of nonzero entries in the representation of each column of X.
%              This means that each column of A must have no more than s nonzero entries.
%    lambda .. Penalty term
%    mode .... Selects minimization problem to be solved:
%     0) min_A ||X-D*A||_2 s.t. ||A_i||_1 <= lambda
%     1) min_A ||A||_1 s.t. ||X_i - D*A_i|| <= lambda
%     2) min_A 0.5*||x-D*A||_2^2 + lambda*||A||
%    pos ...... Look for a positive solution (only positive coefficients)
%
function A=scLasso(X,D,L,lambda,mode,pos)
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

if exist('mode','var')
  params.mode = mode;
end
if exist('pos','var')
  params.pos = pos;
end
params.lambda = lambda;
A=mexLasso(X,D,params);
end
