function spar = mdlsSparseness(A,thres)
if ~exist('thres','var')
  thres = 0;
end
% 
% sparseness = Sparseness(A)
%
% Name: Sparseness
%
% Category: 
%
% Description: Measure the sparseness of the coefficient's matrix.
%
% Input: 
% A ........ coefficient's matrix
%
% Output:
% spar ..... sparseness level
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%

[K,N] = size(A);
spar = sum(sum(abs(A)>thres))/K/N*100;

end % function
