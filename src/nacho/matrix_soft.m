%
% Vector thresholding function.
% y = \arg \min_z (1\2) ||z-x||_2^2 + \tau ||z||_2
%
function Y = matrix_soft(X,tau)
  thenorm = norm(X,'fro');
  Y = max([thenorm-tau,0])/(thenorm+eps)*X;
end
