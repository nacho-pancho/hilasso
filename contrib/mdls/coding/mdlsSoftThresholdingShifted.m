function out = SoftThresholdingShifted(in,mu,c)
%
% out = SoftThresholding(in,mu)
%
% Name: SoftThresholding
%
% Category: core function
%
% Description: Variant of the SoftThresholding algorithm where
% the target function to minimize is 0.5(x-x0)+\lambda|x-c|
% where c is a shift from 0.
% The solution generalizes the original soft thresholding operator 
% in the following way:
%
% \begin{equation}
%   S_\mu^c (x)= \left\{
%   \begin{array}{cl}
%     x + \frac{\mu}{2} & \quad\mbox{if } x\leq c-\frac{\mu}{2}\\
%     c & \quad\mbox{if } |x-c|<\frac{\mu}{2} \\
%     x - \frac{\mu}{2} & \quad\mbox{if } x\geq c+\frac{\mu}{2} 
%   \end{array}
%   \right.
% \end{equation}
%
% Input: 
% in ....... input signal
% mu ....... \mu parameter in the filter. May be scalar or the same size as in
% c ........ shift (optional, 0 if not specified)
%
% Output:
% out ...... filtered signal
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
if nargin < 3
 c = 0
end

mu2=mu/2; clear mu;
out = in;
% 0 & \quad\mbox{if } |x|<\frac{\mu}{2} \\
out(abs(out-c) < mu2) = c;
% x - \frac{\mu}{2} & \quad\mbox{if } x\geq \frac{\mu}{2} 
ind = find (out >= c+mu2);
out(ind) = out(ind) - mu2;
% x + \frac{\mu}{2} & \quad\mbox{if } x\leq -\frac{\mu}{2}\\
ind = find (out <= c-mu2);
out(ind) = out(ind) + mu2;

end

