%
% Evaluate the Discretized Jeffreys of Exponential (Joe) prior
%
% Given a set of quantization points, the function returns the 
% probability mass function associated to a random variable which is
% a discretized version of a continuous-valued r.v. having a JOE density:
%
%                1     1
% Joe(x;k,b)= ------- --- [ e^(-ax) - e^(-bx) ]
%             ln(b/a)  x 
% 
% The integral involves the use of the Exponential Integral (Matlab expint()) function.
% This is explained next:
%
% To evaluate \int_{x0}^{x1}{x^{-1}e^{-ax}} we do it as
% F_a(x0,x1) = \int_{x0}^{inf}{x^{-1}e^{-ax}} - \int_{x1}^{inf}{x^{-1}e^{-ax}} 
%          = expint(a*x0)
% 
% Now, the integral of the PDF in the interval [x_0,x_1] is
% (ln(b/a))^(-1) * ( F_a(x_0,x_1) - F_b(x_0,x_1) )
%
% inputs: 
% x ........ points where the discrete distribution is to be evaluated.
%            the function assumes uniform quantization, that is, that the
%            points in x are equispaced and increasing.
% a ........ beginning of support interval for x 
% b ........ end of support interval for x
% sym ...... use symmetric distribution (mixture of laplacians)
% outputs:
% f ........ the PDF evaluated at the points specified in a
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function f=mdlsJoeDisc(x,a,b,sym)
    if ~exist('sym','var')
        sym = false;
    end
    N = 1/log(b/a);
    dx = (x(2)-x(1))/2;
    x1 = x - dx;
    x2 = x + dx;
    x2(end) = realmax/(10*b);    
    if ~sym
        x1(1) = 0;        
        f = expint(a*x(2:end)) - expint(a*x(1:end-1)) + ... 
            expint(b*x(2:end)) - expint(b*x(1:end-1));
        f = f*N;
    else
        N = N/2;
        x1(1) = -realmax/(10*b);
        ineg = max(find(x2 < 0));
        ipos = (ineg+2):length(x);
        izero = ineg+1; % special treating
        ineg = 1:ineg;
        %
        % positive values
        %
        f(ipos) = N*(expint(b*x2(ipos)) - expint(b*x1(ipos)) + ... 
                     expint(a*x1(ipos)) - expint(a*x2(ipos)));
        %
        % negative values
        %
        f(ineg) = N*(expint(-b*x1(ineg)) - expint(-b*x2(ineg)) + ... 
                     expint(-a*x2(ineg)) - expint(-a*x1(ineg)));
        %
        % at zero
        %
        f(izero) = 1 - N*(-expint(b*x2(izero))-expint(-b*x1(izero))+ ...
                           expint(a*x2(izero))+expint(-a*x1(izero)));
    end
end
 
