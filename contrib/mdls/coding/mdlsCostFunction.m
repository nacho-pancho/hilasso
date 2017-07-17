%
% function [F,R,C,T]=mdlsCostFunction(X,D,A,params)
%
% Name: CostFunction
%
% Category: core function
%
% Description: computes the cost function of the modeling problem. 
%
% Input:
% X ........ data to be fitted
% D ........ dictionary (n x K)
% A ........ coefficients (K x N)
% params ... model parameters. Struct with at least these members:
%    sigma ..... noise std. dev. assuming Gaussian corrupted X
%    lambda .... L1 regularization parameter
%    mu ........ self incoherence term
%    eta ....... normalization term
%    xmu ....... if multiple dictionaries are learned, cross-coherence
%                term (not implemented yet)
%
% Output:
% cost ..... struct containing the following terms:
%   T ........ Total cost (only mandatory one)
%   F ........ fitting term
%   R ........ coefficient regularization term
%   C ........ dictionary coherence penalty term
%   N ........ atom norm term
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Idn$
%
function cost=mdlsCostFunction(X, D, A, params)
  cost = struct();
  K = size(D,2);
  %
  % fitting term: ||X - DA||_F^2 
  %
  t = X - D*A ; 
  cost.F =  0.5 * sum(sum(t.*t));
  %
  % dictionary regularization term: ||D'D - I||_F^2
  %
  if params.mu > 0
      t = D'*D - eye(K); 
      cost.C = params.mu * sum(sum(t.*t));  
  else
      cost.C = 0;
  end
  %
  % additional term for alternate cost function: C_n \sum_{k}{ (||D_k||^2 - 1)^2 }
  % where Cn is a model parameter which should be big.
  %
  if params.eta > 0
      t = sum(D.*D) - 1; 
      N = params.eta * sum(t.*t);
  else
      N = 0;
  end
  %
  % log penalty
  %
  b = params.betaCoef;
  k = params.kappaCoef;
  if isequal(params.codingMethod,'Lasso')
      R = params.thetaCoef*sum( abs(A(:)) );
  else
      R = (k+1)*sum( log(abs(A(:))+b) );
  end
  if issparse(A)
    R = full(R);    
  end
  cost.T = cost.F + cost.R +  cost.C + cost.N;
end
