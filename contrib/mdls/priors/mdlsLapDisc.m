%
% Evaluate the discretized Laplacian probability mass function
% corresponding to the continuous Laplacian density function:
%
% Laplacian(x;lambda)=(lambda/2)*exp(-lambda*|x|)
%   
% a set of *equally spaced* bin centers x is given, and the
% returned values correspond, for each x, to the integral of
% the  density function between x-dx/2,x+dx/2 (when x is not an extreme point) or
% the corresponding infinity otherwise. This would be the expected
% histogram value for histograms centered at x.
%
% inputs: 
% x ........ argument to the density function, can be scalar, vector or
%            matrix. The returned value p will have the same size.
% lambda ... parameter of the laplacian
%
% outputs:
% f ........ the PDF evaluated at the points specified in a
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function H=mdlsLapDisc(x,lambda,mu)
    if ~exist('mu','var')
        mu = 0;
    end
    if (lambda==Inf)
        H=double(x==0); % delta at 0
        return;
    end
    F0=0;
    dx=(x(2)-x(1))/2;
    H=zeros(1,length(x));
    for i = 1:length(x)
        xi = x(i);
        xr = xi + dx;
        xl = xi - dx;
        if xr < mu
            if i == 1
                H(i) = 0.5*exp(-lambda*(mu-xr));
            else
                H(i) = 0.5*(exp(-lambda*(mu-xr)) - exp(-lambda*(mu-xl)));
            end
        else
            if i == length(x)
                H(i) = 0.5*exp(-lambda*(xl-mu));
            else 
                if xl > 0
                    H(i) = 0.5*(exp(-lambda*(xl-mu)) - exp(-lambda*(xr-mu)));
                else
                    H(i) = 1 - 0.5*exp(-lambda*(xr-mu)) - 0.5*exp(-lambda*(mu-xl));
                end
            end
        end
    end
end

 
