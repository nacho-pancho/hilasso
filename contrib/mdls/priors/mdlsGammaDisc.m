%
% Evaluate the discretized Gamma probability mass function
% corresponding to the continuous Gamma density function:
%
%                  x^(kappa-1)*beta^kappa*e^(-beta*x)
% Gamma(x;k,b)= - ---------------
%                 Gamma(kappa)
%   
% a set of *equally spaced* bin centers x is given, and the
% returned values correspond, for each x, to the integral of
% the Gamma density function between x-dx/2,x+dx/2 (when x is not an extreme point) or
% the corresponding infinity otherwise. This would be the expected
% histogram value for histograms centered at x.
%
% inputs: 
% x ........ argument to the density function, can be scalar, vector or
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
function H=mdlsGammaDisc(x,k,b)
    F0=0;
    dx=(x(2)-x(1))/2;
    H=zeros(1,length(x));
    for i=1:(length(x)-1)
        xi=x(i);
        F1 = gammainc(b*(xi+dx),k);
        H(i) = F1 - F0;
        F0 = F1;
    end
    H(end)=1-F0;
end

 
