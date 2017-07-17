%
% Evaluate the Jeffreys of Laplacian (JOL) prior
%
%                1     1
% JOL(x;k,b)= ------- --- [ e^(-a|x|) - e^(-b|x|) ]
%             2ln(b/a)|x| 
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
function f=mdlsJolEval(x,a,b)
    x = abs(x);
    f = (.5/log(b/a))*(exp(-a*x)-exp(-b*x))./x;
    f(x == 0) = 0.5*(b-a)/log(b/a);
end
 
