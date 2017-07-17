%
% Wrapper for using Julien's mexLasso as if we were calling the old scLasso
%
% function [A]=scOMP(X,D,L,l2err)
% 
%   Purpose: 
%     Solves various L1-regularized problems (see 'mode') using LARS.
% 
%     See B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, "Least angle regression," Annals of Statistics, 32(2):407--499, 2004.
%  
%   Inputs:
%    X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%    D ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%    L ....... Sparsity level: indicates maximum number of nonzero entries in the representation of each column of X.
%              This means that each column of A must have no more than s nonzero entries.
%    eps ..... Maximum L2 error in the reconstruction of each vector.
%
function A=scOMP(X,D,L,l2err)
params=struct();
if (L > 0) && (L < size(X,1))
   params.L = L;
else
   params.L = size(X,1);
end
if exist('l2err','var')
  params.eps = l2err;
else
  params.eps = eps;
end
A=mexOMP(X,D,params);
end
