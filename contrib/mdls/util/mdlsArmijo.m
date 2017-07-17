%
% function [x,s] = mdlsArmijo(f, x0, g0, d0, a, b, s0, m)
%
% Updates the solution to a minimization problem of the form x = arg min f(X)
% using the Armijo rule. This is an approximate line search algorithm
% which, given a descent direction d0, a function f and its gradient g0, obtains 
% a step size so that the decrease in the value of the target function is
% 'good enough', meaning it is a above a fixed percentage of the decrease
% predicted by a linear approximation of the function around the current
% iterate, that is
%
% s : ( f(x0) - f(x0 + s0*b^m * d0) ) >= -a * s0 * b^m * <g0,d0>
%
% (note that the second term is positive if d0 is a descent direction,
% that is, <g0,d0> < 0).
%
% input:
%
% f ........ (function handle) for a function of the form f(x), x is m x 1
% x0 ....... (m x 1) current iterate
% g0 ....... (m x 1) gradient at current iterate
% d0 ....... (m x 1) descent direction at current iterate
%     a ........ (0,1] 'satisfaction' factor. The greater required, the 'tougher'
%                the algorithm will be with the final step. Default: 0.1
%     b ........ (0,1) geometric decrease of s at each Armijo
%                iteration. Default: 0.5
%     s0 ....... initial step. Default 1.0
%     m ........ maximum number of Armijo iterations. If this number is
%                reached, s = 0 is returned indicating convergence was
%                achieved. Default: m = 20.
%
% output:
%
% x ......... updated iterate
% s ......... final stepsize
%
function [x,s] = mdlsArmijo(f, x0, g0, d0, a, b, m, s0)
    % 
    % parse input
    %
    if nargin < 5
        a = 0.1;
    end
    if nargin < 6
        b = 0.5;
    end
    if nargin < 7
        m = 20;
    end
    if nargin < 8
        s0 = 1;
    end
    %    
    % algorithm
    %
    f0 = f(x0); 
    s = s0;
    acceptable_decrease = -a*g'*d;
    for j=1:m
        x = x0 + s*d0;
        if  f0 - f(x) >= acceptable_decrease
            return;
        end
        % no luck: reduce step and acceptable decrease
        s = s * b;
        acceptable_decrease = acceptable_decrease * b;
    end
    %
    % failed miserably
    %
    s = 0;
end