function [Xr,A,v] = HiLassoMethod(Y,D,lambda1,lambda2,tol)

% Parameters
n = size(Y,1);
N = size(Y,2);
numSig = length(D);
sizeD = size(D{1},2);

% Construct combined dictionary
Do = [];
groups = [];
for i=1:numSig
    Do = [Do D{i}];
    groups = [groups i*ones(1,sizeD)];
end

if ~exist('tol','var')
    tol = 1e-6;
end
%keyboard
% Solve constrained L1
A = zeros(size(Do,2),N);
for i=1:N
    %A0 = inv(Do'*Do+0.1*eye(size(Do,2)))*Do'*Y(:,i);
    A0 = A(:,i);
    A(:,i) = HiLasso(Y(:,i),Do,A(:,i),groups,lambda1,lambda2,tol,300, 5);
    fprintf('.')
end
fprintf('\n')
%as = group_act_set(A);
if mean(sum(A~=0))<n
    [Xd,A] = compute_ols(Y,Do,A);
end

% Recovered signals
for i=1:numSig
    idx = find(groups == i);
     Xr{i} = D{i}*A(idx,:);
%    Xr{i} = compute_ols(Z,Do,A)
     v(i) = sum(sum((A(idx,:).^2)));
end
