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
% steps ..... desired number of approximation steps. 0 (default)
% for automatic.
% eta ....... quasi-norm auxiliary term. Default is 1e-5
% method .... optimization method, not used yet (only fixed point
% implemented).
%
% outputs:
% alpha ..... d x n matrix with the estimated reconstruction
% coefficients for X using dictionary D.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [alpha,e]=mdlsSparseCodingL1Step(X,D,alpha, sparsity,lambda,thres, eta_e,eta_a, method)

if ~exist('lambda','var')
  lambda=1; % reasonable??
end
if ~exist('thres','var')
  thres=1e-10; % reasonable??
end
if ~exist('method','var')
  method=0;
end
if ~exist('eta_e','var')
  eta_e=1e-5;
end
if ~exist('eta_a','var')
  eta_a=eta_e;
end

if sparsity <= 0
  sparsity = k;
end

[m,k]=size(D);
[m2,n]=size(X);
if m2 ~= m
  error('Dictionary and data are of incompatible dimensions\n');
end

%
% pseudo-norm term: should be as small as possible 
% without causing the condition of A to grow too much
%
e=0;
for j=1:n
  alpha_j=alpha(:,j);
  E=X(:,j)-D*alpha_j;
  e=e+sum(abs(E));
  E=1./sqrt(E.*E+eta_e);
  D2=D.*(E*ones(1,k));
  A=D2'*D; % singular since D is d x k, d < k,  and D2'*D is k x k
  A=A+lambda*diag(1./sqrt(alpha_j.*alpha_j+eta_a));
  b=D2'*X(:,j);
  %fprintf('j=%05d: cond(A)=%f max(A)=%f min(A)=%f max(b)=%f min(b)=%f\n',...
  %   j,cond(A),max(max(A)),min(min(A)),max(b),min(b));
  alpha_j=A\b;
  % optional thresholding
  % sparsity
  if sparsity < k 
    [dummy,sortidx]=sort(abs(alpha_j),'descend');
    alpha_j(sortidx(sparsity+1:end))=0; % keep only /sparsity/ nonzero values.
  end
  alpha(:,j)=alpha_j;
end
%
% do the thresholding only at the end
%
if thres > 0
  to_discard = abs(alpha) < thres;  % logical index
  alpha(to_discard)=0;
end
