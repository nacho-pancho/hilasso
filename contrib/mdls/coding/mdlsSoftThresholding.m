function out = SoftThresholding(in,mu)
%
% out = SoftThresholding(in,mu)
%
% Name: SoftThresholding
%
% Category: core function
%
% Description: Filter the signal with the nonlinear (soft) thresholding
% \begin{equation}
%   S_\mu (x)= \left\{
%   \begin{array}{cl}
%     x + \mu & \quad\mbox{if } x\leq - \mu\\
%     0 & \quad\mbox{if } |x|< \mu \\
%     x - \mu & \quad\mbox{if } x\geq \mu
%   \end{array}
%   \right.
% \end{equation}
%
% Input: 
% in ....... input signal
% mu ....... \mu parameter in the filter
%
% Output:
% out ...... filtered signal
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
  out = in;
  % 0 & \quad\mbox{if } |x|<\frac{\mu}{2} \\
  out(abs(out) < mu) = 0;
  % x - \frac{\mu}{2} & \quad\mbox{if } x\geq \frac{\mu}{2}
  ind = find (out >= mu);
  if length(mu) > 1
    out(ind) = out(ind) - mu(ind);
  else
    out(ind) = out(ind) - mu;
  end
  % x + \frac{\mu}{2} & \quad\mbox{if } x\leq -\frac{\mu}{2}\\
  ind = find (out <= -mu);
  if length(mu) > 1
    out(ind) = out(ind) + mu(ind);
  else
    out(ind) = out(ind) + mu;
  end

end

