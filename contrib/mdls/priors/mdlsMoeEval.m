%
% Evaluate the Mixture of Exponentials (MOE) probability density function, which is
% Gamma(k,b) prior.
%
%              k        1
% f(x;k,b)= - ------------------
%              b (x/b + 1)^(k+1)
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
function f=mdlsMoeEval(x,k,b,sym)
    if ~exist('sym','var')
        sym = false;
    end
    if ~sym
        f = (b^k * k) ./ (x + b).^(k+1);
    else
        f = (0.5 * b^k * k) ./ (abs(x) + b).^(k+1);
    end
end
 
