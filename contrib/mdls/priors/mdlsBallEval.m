%
% Evaluate the "ball" continuous probability density function.
%
% This models the marginal distribution of a any single coordinate of a vector
% when the vectors are uniformly distributed onto the surface of a
% d-dimensional unit sphere centered at 0.
% Note: This turns out to be a nice prior for the marginal distribution of dictionary coefficients.
%
% inputs: 
% a ........ argument to the density function, can be scalar, vector or
% matrix. The returned value p will have the same size.
% d ........ dimension of the Sphere
%
% outputs:
% p ........ the PDF evaluated at the points specified in a
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function p=mdlsBallEval(a,d)
  p=1/sqrt(pi)*gamma(d/2)/gamma((d-1)/2)*(1-a.^2).^(d/2-1.5);
  p=p.*(a>-1).*(a<1);
  
