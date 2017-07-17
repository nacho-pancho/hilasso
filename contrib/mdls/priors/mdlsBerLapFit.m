%
% Estimate parameters of a Bernoulli-Laplace (BL) distribution.
%
% A random variable X distributed according to the BL distribution
% is modeled as X=Z*Y, where Z ~Bernoulli(theta) and Y ~ Laplace(mu,lambda)
% this is a continuos R.V. with a density function which is discontinuous
% at zero, having a Kronecker delta of magnitude 1-theta. also f(X|Z=1) 
% is a continuos density distribution of the form lambda/2*e^(-lambda*|X|).
%
% This function estimates the parameters theta, lambda and mu given a
% vector X. If X is a matriz, it gives a set of estimations for each row in
% X.
%
% inputs: 
% X ......... a d x n matrix representing n d-dimensional vectors
%
% outputs:
%
% theta ..... Bernoulli parameter
% lambda .... Laplacian scale
% mu ........ Laplacian mean
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [theta, lambda, mu]=mdlsBerLapFit(X)
[m,n]=size(X);
if (m > 1) && (n == 1)
    X=X';
    [m,n]=size(X)
end
lambda=zeros(m,1); mu=zeros(m,1);
for i=1:m
    nzi=find(X(i,:) ~= 0);
    n1=length(nzi);
    NZX=full(X(i,nzi));
    clear nzi;
    theta(i)=n1/n;
    [lambda(i,1), mu(i,1)]=mdlsLaplacianFit(NZX);
end


    
  
