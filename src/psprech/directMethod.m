function [Xr,A] = directMethod(Y,D,sigma)


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


% Solve OMP
lambda = (1+sigma/6)*sigma^2;
 A=scOMP(Y,Do,size(Do,2),lambda);

 
% Recovered signals
for i=1:numSig
    idx = find(groups == i);
    Xr{i} = D{i}*A(idx,:);
end


