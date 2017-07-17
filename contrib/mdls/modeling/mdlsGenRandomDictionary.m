function D = mdlsGenRandomDictionary(n,K,seed)
% 
% D = GenerateRandomDictionary(n,K,varargin)
%
% Name: GenerateRandomDictionary
%
% Category:
%
% Description: Create an initial dictionary with random gaussian entries
%
% Input:
% n ........ dimension of the atoms in the dictionary.
% K ........ number of atoms in the dictionary.
%
% Output:
% D ........ the dictionary as a (n x K) matrix.
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
  if ~exist('seed','var')
      seed = 0;
  end
  if seed == 0
      seed = n*100+K;
  end
  randn('state',seed);
  D=randn(n,K);
  for l=1:K
    D(:,l)=D(:,l)-mean(D(:,l)); % zero mean
    D(:,l)=D(:,l)/norm(D(:,l)); % unit norm
  end
end