%
% Estimate parameters of a Bernoulli-Gamma-Laplace (BGL) distribution.
%
% A random variable X distributed according to the BGL distribution
%
% This function estimates the parameters theta, k and b given a
% vector X. If X is a matriz, it gives a set of estimations for each row in
% X.
%
% inputs: 
%
% outputs:
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [theta, k, b]=mdlsBerMOLFit(X,k0,method)
if ~exist('k','var')
    k0=-1;
end
if ~exist('method','var')
    method='moments';
end

[m n]=size(X);
if (m>1) && (n == 1)
    X=X';
    [m,n]=size(X);
end

k=zeros(m,1); b=zeros(m,1); theta=zeros(m,1);
for i=1:m
  nzi=find(X(i,:) ~= 0);
  n1=length(nzi);
  NZX=full(X(i,nzi));
  clear nzi;
  theta(i)=n1/n;
  if n1 == 0
      warning(sprintf('%d:no nonzero elements!\n',i));
      k(i)=1; b(i)=realmax;
  else
      [k(i) b(i)]=mdlsMOLFit(NZX,k0,method);
  end
end
    
  
