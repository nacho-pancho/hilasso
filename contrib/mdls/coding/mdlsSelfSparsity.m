%
% function A=selfSparsity(D,metric,sparsity,lambda)
%
% purpose:
% measure some structural properties of a dictionary of atoms:
% 1) self-sparsity: how sparsely can each atom be represented in
% terms of the other atoms in the dictionary?  This is given as a
% k x k matrix similar to alpha but where the diagonal elements
% are 0 (inhibit self-contribution,  otherwise the trivial
% solution (eye) would be the optimal reconstruction).
%
% input:
% D .......... Dictionary to be studied
% metric ..... Metric to use in penalty for decomposition: 0 (OMP), 1 (mine) or
%              2 (LASSO). 
% sparsity ... Desired sparsity level (forced).
% lambda ..... Penalty term.
%
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function A=mdlsSelfSparsity(D,metric,sparsity,lambda)

if nargin ~= 4
  error('Must specify all parameters\n');
end
[d,k]=size(D);

if sparsity <= 0
  sparsity = k;
end

A=zeros(k,k);

for kk=1:k
  Y=D(:,kk);
  D2=D(:, [1:(kk-1) (kk+1):k] );
  if metric == 0
    a=mexOMP(Y,D2,sparsity,0); % 0 reconstruction error
  elseif metric == 1
    a=sparseCodingL1(Y,D2,sparsity,lambda);
  elseif metric == 2
    a=mexLasso(Y,D2,sparsity,lambda,2); % lamb 
  else
    error('Invalid metric specified.');
  end
  %
  % put coefficients back in A
  % the place of the 'self' coefficient (the diagonal) is left
  % with a 0.
  %
  %a=full(a);
  if kk > 1
    A(1:(kk-1),kk)=a(1:(kk-1));
  end
  if kk < k
    A((kk+1):k,kk)=a(kk:end);
  end
end
