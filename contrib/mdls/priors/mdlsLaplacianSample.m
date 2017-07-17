% 
% Draw N samples from a Laplacian(mu,lambda) distribution.
% returns a sample from X=Z*L where L is distributed as a Laplacian(mu,lambda), i.e.,
% f(L)=lambda/2*e^(-lambda*|L-mu|)
% 
% inputs: 
% N ........ number of samples to draw
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
function X=mdlsLaplacianSample(N,lambda,mu)
if nargin < 4
    mu=0;
end
X=2*(rand(1,N)>0.5)-1;
A=rand(1,N);
X=X.*(-log(1-A)/lambda) + mu;
