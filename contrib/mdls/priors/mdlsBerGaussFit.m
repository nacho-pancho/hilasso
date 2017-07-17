%
% Estimate parameters of a Bernoulli-Gaussian(BG) distribution.
%
% A random variable X distributed according to the BG distribution
% is modeled as X=Z*Y, where Z ~Bernoulli(theta) and Y ~
% Gaussian(mu,lambda).
% This is a continuos R.V. with a density function which is discontinuous
% at zero, having a Kronecker delta of magnitude 1-theta. also f(X|Z=1) 
% is a continuos density distribution of the form 1/sqrt(2*pi*sigma^2)*e^(-(X/sigma)^2).
%
% This function estimates the parameters theta, sigma and mu given a
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
function [theta, sigma, mu]=mdlsBerGaussFit(X)
S=sum(X.*X,2);
n=size(X,2);
n0=sum(X == 0,2);
n1=n-n0;
%
% very reliable because the EXACT value 0 has probability theta due only to
% the Bernoulli. The continuos distribution has probability 0 for
% any measure 0 point in R.
%
theta=n1./n;
m=size(X,1);
  %
  % as the probability of x given that Z=1 is a Gaussia of parameters mu,sigma, 
  % their estimators are  just the usual ML estimators f
if nargout < 3
  % assume mean in Laplacian is 0
  sigma=sqrt(S./(n1-1));
  return;
else
  mu=zeros(m,1);
  sigma=zeros(m,1);
  for i=1:m
    apos=X(i,find(X(i,:)));
    mu(i)=mean(apos);
    sigma(i)=sqrt(var(apos));
  end
end


    
  
