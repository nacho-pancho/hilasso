function [Xr,A,v,Xo,Ao,vo] = HiLassoColMethodW(Y,D,A0,lambda1,lambda2,lambda,...
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
%
% Solve collaborative Hilasso
%
%
% find weights for L2 term
%
Aols = Do'*inv(Do*Do')*Y;
ge   = group_energy(Aols,groups);
ge   = sqrt(mean(ge'));
mge  = sqrt(mean(ge));
lambda2 = lambda2*(mge./(ge+eps))
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
end

% comute OLS with the detected active set
recompute = 0;
[Xo,Ao,vo] = refineSolution(Y,Ao,D,Do,groups,lambda,H,recompute);
if lambda == 0
    Xr = Xo;
    v  = vo;
end
end
