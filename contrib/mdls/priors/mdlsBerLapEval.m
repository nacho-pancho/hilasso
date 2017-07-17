%
% purpose:
% evaluate the BL distribution at the given points.
%
% inputs:
% x ........ evaluation point(s)
% theta .... Bernoulli parameter, i.e., p(Z=1)
% lambda ... shape parameter of the Laplacian
% mu ....... centroid of Laplacian
%
% outputs:
% fx ....... the distribution evaluated at x
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function fx=mdlsBerLapEval(x,theta,lambda,mu)
if nargin < 4
    mu=0;
end
 fx=theta*lambda/2*exp(-lambda*abs(x-mu));
 i0=find(x==0);
 if ~isempty(i0)
     fx(i0)=fx(i0)+(1-theta);
 end


