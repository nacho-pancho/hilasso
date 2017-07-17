%
% Evaluate the Jeffreys of Exponential (Joe) prior
%
%                1     1
% Joe(x;k,b)= ------- --- [ e^(-ax) - e^(-bx) ]
%             ln(b/a)  x 
% 
% inputs: 
% x ........ argument to the density function, can be scalar, vector or
%            matrix. The returned value p will have the same size.
% a ........ beginning of support interval for x 
% b ........ end of support interval for x
%
% outputs:
% f ........ the PDF evaluated at the points specified in a
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function f=mdlsJoeEval(x,a,b,sym)
    C = 1/log(b/a);
    if exist('sym','var') 
        if sym
            x = abs(x);
            C = C/2;
        end
    else
        sym = ~isempty(x<0);
        C = C/2;
    end
    f = C*(exp(-a*x)-exp(-b*x))./x;
    f(x == 0) = C*(b-a);
end
 
