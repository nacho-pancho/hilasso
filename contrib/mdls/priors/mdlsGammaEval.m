%
% Evaluate the Gamma probability density function,
%
%                  x^(kappa-1)*beta^kappa*e^(-beta*x)
% Gamma(x;k,b)= - ---------------
%                 Gamma(kappa)
%   
% inputs: 
% a ........ argument to the density function, can be scalar, vector or
% matrix. The returned value p will have the same size.
% k ........ k, the shape parameter of the Gamma prior (k>=1)
% b ........ beta, the scale parameter of the Gamma prior.
%
% outputs:
% f ........ the PDF evaluated at the points specified in a
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function f=mdlsGammaEval(x,k,b)
    f= (b^k/gamma(k)) * x.^(k-1) .* exp(-b*x);
 end
 
