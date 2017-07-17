%
% purpose:
% evaluate the Gaussian distribution at the given points.
%
% inputs:
% x ........ evaluation point(s)
% sigma2 ... variance
% mu ....... mean
%
% outputs:
% fx ....... the distribution evaluated at x
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function fx=mdlsGaussEval(x,sigma2,mu)
if nargin < 3
    mu=0;
end
 fx=1/sqrt(2*pi*sigma2)*exp((-0.5/sigma2)*(x-mu).^2);
end

