%
% function [Xr,A,v,Xo,Ao,vo] = HiLassoColMethod2(Y,D,A0,lambda1,lambda2,lambda,...
%                                               tol,H,max_iter,c)
%
% INPUT
% Y ............ (MxN) mixed signal to analyze
% D ............ (1xNC cell) each element is an MxKi dictionary. KT := \sum_i{Ki}
% A0 ........... (KTxN) initial coefficients. If empty, start with 0.
% lambda1 ...... (1x1)  L1 penalty
% lambda2 ...... (1x1)  L2 penalty
% lambda ....... (1x1)  L1 penalty for second stage lasso. 0 means no
%                 second stage (0.1)
% tol .......... (1x1) tolerance (1e-3)
% H ............ (MxN) missing data mask. 1 means missing sample
% max_iter ..... (1x1) maximum iterations (100)
% c ............ (1x1) augmented lagrangian parameter (1)
%
% OUTPUT
% Xr
% A
% v
% Xo
% Ao
% vo
%
function [Xr,A,v,Xo,Ao,vo] = HiLassoColMethod2(Y,D,A0,lambda1,lambda2,lambda,...
                                               tol,H,max_iter,c)

if ~exist('lambda','var')
    lambda = 0.1;
end

if ~exist('tol','var')
    tol = 0.001;
end

if ~exist('max_iter','var')
    max_iter = 100;
end

if ~exist('H','var')
    H = [];
end

if ~exist('c','var')
    c = 1;
end

if isempty(H)
    H = ones(size(Y));
end

stop_criterion = 5;

% Parameters
n = size(Y,1);
N = size(Y,2);
numSig = length(D);
sizeD = size(D{1},2);

% Construct combined dictionary
Do = [];
groups = [];

for i=1:numSig
    ng(i) = size(D{i},2);
    Do = [Do D{i}];
    groups = [groups i*ones(1,size(D{i},2))];
end

% Solve collaborative Hilasso
if isempty(A0)
    A0 = zeros(size(Do,2),N);
end
[A,obj,times,seq] = HiLassoCollaborative(Y,Do,A0,groups,lambda1,lambda2, ...
                                          tol,max_iter, stop_criterion ,H, ...
                                          c);
Ao = A;
if lambda > 0
  % compute OLS with the detected dictionaries
  recompute = 1;
  [Xr,A,v] = refineSolution(Y,A,D,Do,groups,lambda,H,recompute);

  % comute OLS with the detected active set
  recompute = 0;
  [Xo,Ao,vo] = refineSolution(Y,Ao,D,Do,groups,lambda,H,recompute);
  if lambda == 0
      Xr = Xo;
      v  = vo;
  end
else
  Xr = [];
  v  = [];
  vo = [];
  Xo = [];
end
end
