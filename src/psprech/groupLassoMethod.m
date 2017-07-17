function [Xr,A] = groupLassoMethod(Y,D,A0,lambda)


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

% Define SpaRSA functions & parameters
psi = @(x,tau) group_vector_soft(x,tau,groups);
phi = @(x) group_l2norm(x,groups);
stop_crit = 5;
tolerance = 0.001;
max_iter = 300;


% Solve constrained L1
if isempty(A0)
   A = zeros(size(Do,2),N);
else
   A = A0;
end
for i=1:N

    A(:,i)= SpaRSA(Y(:,i),Do,lambda,...
        'Psi',psi,...
        'Phi',phi,...
        'Monotone',1,...
        'Debias',0,...
        'StopCriterion',stop_crit,...
        'ToleranceA',tolerance, ...
        'MaxiterA',max_iter,...
        'Verbose',false,...
        'Initialization',A(:,i));
    fprintf('.')
end
fprintf('\n')


% Recovered signals
for i=1:numSig
    idx = find(groups == i);
     Xr{i} = D{i}*A(idx,:);
%    Xr{i} = compute_ols(Z,Do,A)
end  


    


