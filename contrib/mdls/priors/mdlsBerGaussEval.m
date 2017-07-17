%
% purpose:
% evaluate the BG distribution at the given points.
%
% inputs:
% theta .... Bernoulli parameter, i.e., p(Z=1)
% sigma ... shape parameter of the Gaussian (i.e., the std. dev.)
% mu ....... centroid of Laplacian
%
% outputs:
% fx ....... the distribution evaluated at x
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function fx=mdlsBerGaussEval(x,theta,sigma,mu)
if nargin < 3
    mu=0;
end
 fx=theta*pdf('norm',x,mu,sigma);
 i0=find(x==0);
 if ~isempty(i0)
     fx(i0)=fx(i0)+(1-theta);
 end


