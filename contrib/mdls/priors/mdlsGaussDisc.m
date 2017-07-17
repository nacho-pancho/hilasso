%
% Evaluate the discretized Gaussian probability mass function
% corresponding to the continuous Gaussian density function with
% variance sigma2
%   
% a set of *equally spaced* bin centers x is given, and the
% returned values correspond, for each x, to the integral of
% the density function between x-dx/2,x+dx/2 (when x is not an extreme point) or
% the corresponding infinity otherwise. This would be the expected
% histogram value for histograms centered at x.
%
% inputs: 
% x ........ argument to the density function, can be scalar, vector or
%            matrix. The returned value p will have the same size.
% sigma2 ... variance of the Gaussian
%
% outputs:
% f ........ the PDF evaluated at the points specified in a
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function H=mdlsGaussDisc(x,sigma2,mu)
    if ~exist('mu','var')
        mu = 0;
    end
    F0=0;
    dx=(x(2)-x(1))/2;
    H=zeros(1,length(x));
    scale=1/sqrt(2*sigma2);
    for i=1:(length(x)-1)
        xi=x(i);
        xe = xi+dx-mu;
        if xe < 0
            F1 = 0.5*(1-erf(-xe*scale));
        else
            F1 = 0.5*(1+erf(xe*scale));
        end
        H(i) = F1 - F0;
        F0 = F1;
    end
    H(end)=1-F0;
end

 
