%
% function alpha=sparseCodingL1(X,D,sparsity,lambda)
% 
% purpose:
% 
% given data X and dictionary D, obtain the optimal solution to
%
% ||X-D*alpha||_1 + lambda*||alpha||_1
%
% inputs:
% X ......... d x n matrix representing n d-dimensional vectors
% D ......... d x k matrix representing a dictionary of k d-dim. atoms.
% sparsity .. maximum number of nonzero coefficients in
%             representation. 0 means no sparsity constrains.
% lambda .... penalty term for L1 norm of alpha
% thres ..... hard thresholding threshold of coefficients to weed
%             out small values. 0 means no thresholding.
% eta ....... quasi-norm auxiliary term. Default is 1e-2
% method .... optimization method. Defined values are:
%             0) eta is the value fixed by the user.
%             1) eta is adapted at each iteration so that
%                the pseudo-norm approaches L1. (default).
%                The initial value is the one given.
%             3) same as 1, and lambda is also recomputed as the
%             ML estimator of alpha.
%
% implemented).
% et ........ Energy threshold. If change between iterations
% falls below this value for three iterations, stop.
% maxiter ..... max number of approximation steps. By default 20. 
%
% outputs:
% alpha ..... d x n matrix with the estimated reconstruction
% coefficients for X using dictionary D.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function alpha=mdlsSparseCodingL1(X,D,sparsity,lambda,thres, eta, method, et, steps)
[m,k]=size(D);
[m2,n]=size(X);
alpha=zeros(k,n);
if m2 ~= m
  error('Dictionary and data are of incompatible dimensions\n');
end
if ~exist('sparsity','var')if steps == 0
  steps=30;
end

  sparsity=0; 
end
if ~exist('lambda','var')
  lambda=1; % reasonable??
end
if ~exist('et','var')
  et=0.1;
end
if ~exist('thres','var')
  thres=0; 
end
if ~exist('method','var')
  method=1;
end
eta_e=1e-2;
eta_a=1e-2;
if exist('eta','var')
  eta_e=eta;
  eta_a=eta;
end
if ~exist('steps','var')
  steps=20;
end

if sparsity <= 0
  sparsity = k;
end
E0=Inf;
J = 1;
stalled=0;
while 1 
  [alpha,e]=mdlsSparseCodingL1Step(X,D,alpha, sparsity,lambda,thres, eta_e,eta_a);
  Ea=sum(sum(abs(alpha)))/(k*n);
  Ee=e/(m*n);
  if bitand(method,1)
    eta_a = Ea/100;
    eta_e = Ee/100;
  end
  if bitand(method,2) 
    lambda = 1/Ea;
  end
  E=Ee+lambda*Ea;
  dE=E0-E; E0=E;
  fprintf('iter=%d E(error)=%f E(alpha)=%f E(tot)=%f dE=%f lambda=%f eta_e=%f eta_a=%f\n',J,Ee,Ea,E,dE,lambda,eta_e,eta_a);
  J=J+1;
  if J > steps
    fprintf('Maximum iterations reached.\n');
    return;
  end
  if dE < 0
    fprintf('Minimum reached.\n');
    return;
  end
  if dE < et
    stalled=stalled+1;
    if stalled >= 3
      fprintf('No further improvement.\n');
      return;
    end
  else
    stalled=0;
  end      
end


