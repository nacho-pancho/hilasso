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
    param.L = min([size(Y,1)-1,size(Dr,2)]);
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
    e = 0.001;
    for i=1:N

        Dr = Do(:,A(:,i)>e);
        param.lambda = 0.15*sqrt(size(Dr,2))/sqrt(mean(nAtoms));
        param.L = min([size(Y,1)-1,size(Dr,2)]);
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
