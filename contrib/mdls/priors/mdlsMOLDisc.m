%
% Evaluate the discretized Mixture of Laplacians (MOL) probability mass function
% corresponding to the continuous MOL density function
%
% MOL(x;k,b)=0.5*(k/b)*(|x|/b+1)^(-k-1)
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
% kappa .... shape parameter of MOL
% betta .... scale parameter of MOL
%
% outputs:
% f ........ the PDF evaluated at the points specified in a
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function H=mdlsMOLDisc(x,kappa,betta,mu)
    if ~exist('mu','var')
        mu = 0;
    end
    if (betta==Inf)
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
                H(i) = 0.5*((mu-xr)/betta + 1)^(-kappa);
            else
                H(i) = 0.5*(((mu-xr)/betta + 1)^(-kappa) - ((mu-xl)/betta + 1)^(-kappa));
            end
        else
            if i == length(x)
                H(i) = 0.5*((xr-mu)/betta + 1)^(-kappa);
            else
                if xl > 0
                    H(i) = 0.5*( ((xl-mu)/betta + 1)^(-kappa) - ((xr-mu)/betta + 1)^(-kappa));
                else % xl is negative, xr is positive
                    H(i) = 1 - 0.5*((mu-xl)/betta + 1)^(-kappa) - 0.5*((xr-mu)/betta + 1)^(-kappa);
                end
            end
        end
    end
end

 
