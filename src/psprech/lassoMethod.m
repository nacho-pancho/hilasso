function [Xr,A,v] = lassoMethod(Y,D,lambda,H)

if ~exist('H','var')
    missing = 0;
else
    missing = 1;
end


% Parameters
n = size(Y,1);
N = size(Y,2);
numSig = length(D);

% Construct combined dictionary
Do = [];
groups = [];
for i=1:numSig
    Do = [Do D{i}];
    groups = [groups i*ones(1,size(D{i},2))];
end


% Solve constrained L1
% A=scLasso(Y,Do,length(groups),lambda,1);
if ~missing
    A=scLasso(Y,Do,length(groups),lambda,2);
    if mean(sum(A>0.001))<n
        [Xd,A] = compute_ols(Y,Do,A);
    end
else
    A = zeros(size(Do,2),N);
    for i=1:N
        L = min([ length(groups) size(Y,1)-1 ]);
        A(:,i)=scLasso(Y(H(:,i)==1,i),Do(H(:,i)==1,:),L,lambda,2);
        [Xd,A(:,i)] = compute_ols(Y(H(:,i)==1,i),Do(H(:,i)==1,:),A(:,i));
    end
end

    
% Recovered signals
for i=1:numSig
    idx = find(groups == i);
    Xr{i} = D{i}*A(idx,:);
    v(i) = full(sum(sum((A(idx,:).^2))));
end

