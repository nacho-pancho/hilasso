function [Xr,A,v,Xo,Ao,vo] = HiLassoColMethod(Y,D,lambda1,lambda2,lambda,tol,H,max_iter,c)

if ~exist('lambda','var')
    lambda = 0.1;
end

if ~exist('tol','var')
    tol = 0.0001;
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
A = zeros(size(Do,2),N);
[A,obj,times,seq] = HiLassoCollaborative(Y,Do,A,groups,lambda1,lambda2, ...
                                          tol,max_iter, stop_criterion ,H, ...
                                          c);
Ao = A;

% compute OLS with the detected dictionaries
recompute = 1;
[Xr,A,v] = refineSolution(Y,A,D,Do,groups,lambda,H,recompute);

% comute OLS with the detected active set
recompute = 0;
[Xo,Ao,vo] = refineSolution(Y,Ao,D,Do,groups,lambda,H,recompute);


end


function [Xr,A,v] = refineSolution(Y,A,D,Do,groups,lambda,H,recompute)


% Parameters
n = size(Y,1);
N = size(Y,2);
numSig = length(D);
sizeD = size(D{1},2);


for i=1:numSig
    nAtoms(i) = size(D{i},2);
end
nAtomsCum = [0 cumsum(nAtoms)];

% Agregado para hacer un lasso convencional con los grupos detectados como
% presentes.
Ao = A;
if recompute

    Dr = [];
    for i=1:numSig
        idx = find(groups == i);
        v(i) = sum(sum((A(idx,:).^2)));
        if v(i)> mean(sum(abs(A)))*0.005;
            Dr = [Dr D{i}];
        else
            v(i) = 0;
            A(nAtomsCum(i)+1:nAtomsCum(i+1),:) = 0;
        end
    end

    param.lambda = lambda*sqrt(size(Dr,2))/sqrt(mean(nAtoms));
    param.L = size(Y,1)-1;
    try
        a = mexLasso(Y,Dr,param);
    catch
        % para p
        a = pLasso(Y,Dr,param.lambda);
    end
    com = 0;
    for i=1:numSig
        if v(i)>0
            A(find(groups == i),:) = full(a((com+1):(com+nAtoms(i)),:));
            com = com+nAtoms(i);
        end
    end
    
    for i=1:N
    [Xd,A(:,i)] = compute_ols(Y(H(:,i)==1,i),Do(H(:,i)==1,:),A(:,i));
    end
    
else
    e = mean(abs(A(:)))*0.005;
    for i=1:N

        Dr = Do(:,A(:,i)>e);
        param.lambda = 0.15*sqrt(size(Dr,2))/sqrt(mean(nAtoms));
        try
            a = mexLasso(Y(:,i),Dr,param);
        catch
            % para p
            a = pLasso(Y(:,i),Dr,param.lambda);
        end

        A(A(:,i)>e,i) = a;
        [Xd,A(:,i)] = compute_ols(Y(H(:,i)==1,i),Do(H(:,i)==1,:),A(:,i));
    end

end

% Recovered signals
for i=1:numSig
    idx = find(groups == i);
     Xr{i} = D{i}*A(idx,:);
     v(i) = sum(sum((A(idx,:).^2)));
%    Xr{i} = compute_ols(Z,Do,A)
end

end

function a = pLasso(Y,D,lambda)


  for i=1:size(Y,2)

    a(:,i) = SpaRSA(Y(:,i),D,lambda,...
                'Monotone',1,...
                'Debias',0,...
                'StopCriterion',1,...
                'ToleranceA',0.00000001, ...
                'ToleranceD',0.000001, ...
                'MaxiterA',1000,...
                'Verbose',0);

    end


end

