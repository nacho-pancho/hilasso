%
% purpose: compute the Kullback-Leibler divergence (also known as Divergence or Relative
% Entropy or Kullback-Leibler distance) between two probability
% distributions.
% If p and q are matrices, then a column vector is returned for the KL
% distance of each corresponding row in p and q.
%
% inputs:
% p ..... 'reference' distribution. Note that the KL distance is NOT! symmetric!
% q ..... second distribution. 
% base .. base of logarithm. Defaults to 2.
%
% outputs:
% 
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function D=KLdiv(p,q,base)
if nargin < 3
  base=2;
end
main_dim=2;
if size(p,1) > size(p,2)
  main_dim=1;
end

  % define 0 log 0 = 0 and so skip zero entries in p:
  posidx=p~=0;
  p=p(posidx);
  q=q(posidx);
  D=sum(p.*log(p./q),main_dim);
  D=D/log(base);
