%
% function [F,R,C,T]=CostFunction(X,D,A,params)
%
% Name: CostFunction
%
% Category: core function
%
% Description: computes the cost function of the modeling problem. 
%              not necessarily efficient for evaluation, but
%
% Input:
% X ........ data to be fitted
% D ........ dictionary (n x K)
% A ........ coefficients (K x N)
% params ... model parameters
%
% Output:
% T ........ Total cost (only mandatory one)
% F ........ fitting term
% R ........ coefficient regularization term
% C ........ dictionary coherence penalty term
% N ........ atom norm term
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Idn$
%
function [T,F,R,C,N]=CostFunction(X,D,A,params)
  K = size(D,2);
  %
  % fitting term: ||X - DA||_F^2 
  %
  t = X - D*A; F = sum(sum(t.*t));
  %
  % dictionary regularization term: ||D'D - I||_F^2
  %
  t = D'*D - eye(K); C = sum(sum(t.*t));  
  %
  % additional term for alternate cost function: C_n \sum_{k}{ (||D_k||^2 - 1)^2 }
  % where Cn is a model parameter which should be big.
  %
  t = sum(D.*D) - 1; N = sum(t.*t);
  %
  % log penalty
  %
  b = params.betaCoef;
  k = params.kappaCoef;
  %R = sum(sum( (k+1)*log(abs(A)+b) ));
  R = sum(sum(abs(A)));
  if issparse(A)
    R = full(R);    
  end
  T = params.varError * F + R + params.mu * C + params.eta * N;
end
