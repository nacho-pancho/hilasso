%
% Estimate parameters of a Laplacian  distribution.
%
% A Laplace distribution of parameters mu,lambda has a continuous density function 
% f(x)=lambda/2*e^(-lambda*|x|).
%
% This function estimates the parameters lambda and mu given a
% vector X. If X is a matriz, it gives a set of estimations for each row in
% X.
%
% inputs: 
% X ........ date samples
%
% outputs:
% lambda ... scale parameter of the Laplacian
% mu ....... centroid of Laplacian
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%

function [lambda, mu]=mdlsLaplacianFit(X)
[m,n]=size(X);
if (m > 1) && (n == 1)
    X=X';
    [m,n]=size(X);
end
if nargout < 2
  mu=zeros(m,1);
else
  if m == 1
      mu=median(X);
  else
      mu=median(X,2);
  end
end
S=sum(abs(X-mu*ones(1,n)),2);
lambda=n./S;


  
