% 
% draw samples from a BL(lambda,p) distribution.
% returns a sample from X=Z*L where
% Z ~ Ber(theta)
% L ~ Laplace(lambda) is distributed as a Laplacian(mu,lambda), i.e.,
% f(L)=lambda/2*e^(-lambda*|X-mu|)
%
% 
% inputs: 
% N ........ number of samples to draw
% theta .... Bernoulli parameter, i.e., p(Z=1)
% lambda ... shape parameter of the Laplacian
% mu ....... centroid of Laplacian
%
% outputs:
% x ........ a vector of N samples from the Laplacian(mu,lambda)
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function x=mdlsBerLapSample(N,theta,lambda,mu)
if nargin < 4
    mu=0;
end
x=zeros(1,N);
for i=1:N
    % Bernoulli part
    Z=rand();
    if Z >= theta
        continue;
    end
    S=rand();
    %
    % sign
    %
    if S<0.5
        S=-1;
    else
        S=1;
    end
    %
    % magnitude
    % for a Laplacian
    %
    % F(l)=1-e^-(lambda*l)
    %
    % choose L so that L=F^(-1)(a), where a ~ Uni[0,1]
    A=rand();
    x(i)=S*(-log(1-A)/lambda)+mu;
end
