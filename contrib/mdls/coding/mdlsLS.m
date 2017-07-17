% function A = mdlsLS(X,D,lambda,mode)
%
% Category: core function
%
% Description: compute the regularized Least Squares estimation of reconstruction
%              coefficients A such that X=D*A.
%
% Input:
% X ............... (nxN) N patches of dimension n
% D ............... (nxK) Dictionary of K codewords
% lambda .......... (1x1) Regularization parameter (penalty)
% mode (opt) ...... 1 = Ridge (squared penalty), 2 = pseudoinverse
%
% Output:
% A ............... (KxN) Reconstruction coefficients.
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
function A = mdlsLS(X,D,lambda,mode)
  if ~exist('mode','var')
    mode=1;
  end
  K = size(D,2);
  N = size(X,2);
  A = zeros(K,N);
  if mode == 1
    % forced regularization
    DD = inv(D'*D + lambda*eye(K))*D';
    A=DD*X;
  else
    % 
    [u,s,v] = svd(D'*D);
    d = diag(s); 
    d2 = zeros(1,K);    
    for i = 1:K      
      if (abs(d(i)) > lambda)
        d2(i) = 1/d(i);
      else
        break;
      end
    end
    DD = v*diag(d2)*u';
    A = DD*D'*X;
  end
end
