%
% purpose:
% evaluate the Laplacian distribution at the given points.
%
% inputs:
% x ........ evaluation point(s)
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
function fx=mdlsLaplacianEval(x,lambda,mu)
if nargin < 3
    mu=0;
end
 fx=lambda/2*exp(-lambda*abs(x-mu));
end

