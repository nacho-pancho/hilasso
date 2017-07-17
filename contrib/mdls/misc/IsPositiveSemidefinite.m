% 
% [isPD,iota] = IsPositiveSemidefinite(M)
%
% Name: IsPositiveSemidefinite
%
% Category: auxiliar function
%
% Description: Check if the matrix M is positive semidefinite
%
% Input:
% M ........ Input matrix
%
% Output:
% isPD ..... boolean, true if M is positive semidefinite, false if not.
% iota ..... constant to add for in case is not positive semidefinite, M +
%            iota*eye(size(M,1)
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
function [isPD,iota] = IsPositiveSemidefinite(M)
% Using chol is preferable to using eig for determining positive
% definiteness. [see doc chol].
[R,p] = chol(M);
if isequal(p,0)
  isPD = true;
  iota = 0;
else
  isPD = false;
  eigenvalues = eig(M);
  iota = -min(eigenvalues) + 1e-1;
end
end % function
