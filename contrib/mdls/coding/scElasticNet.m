%
% Wrapper for using Julien's mexLasso as if we were calling the old scLasso
%
% function [A]=scElasticNet(X,D,s,lambda1,lambda2)
% 
%   Purpose: 
%     Solves the elastic net problem
%
%     min_A 0.5*||x-D*A|| + lambda1*||A||_1 + lambda2*||A||_2
% 
%   Inputs:
%    X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%    D ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%    lambda1 .. L1 penalty term
%    lambda2 .. L2 penalty term
%
function A=scElasticNet(X,D,lambda1,lambda2)
  params=struct();
  params.lambda1 = lambda1;
  params.lambda2 = lambda2;
  params.lambda3 = 0;
  params.accelerated = false;
  A=mexLasso(X,D,params);
end
