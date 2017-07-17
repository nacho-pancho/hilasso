% 
% draw samples from a BG(theta,sigma,mu) distribution
% returns a sample from X=Z*S*L where
% Z ~ Ber(theta)
% S ~ Ber(1/2) is the sign of x
% L ~ Laplace(lambda) is distributed as a one-sided laplacian of parameter
% lambda
%
% inputs: 
% N ........ number of samples to draw
% theta .... Bernoulli parameter, i.e., p(Z=1)
% sigma ... shape parameter of the Gaussian (std. dev.)
% mu ....... centroid of Gaussian
%
% outputs:
% x ........ a vector of N samples from the Laplacian(mu,lambda)
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%function x=mdlsBerGaussSample(N,theta,sigma,mu)
%
if nargin < 4
    mu=0;
end
x=sigma*randn(1,N);
x(rand(1,N) >= theta)=0;

