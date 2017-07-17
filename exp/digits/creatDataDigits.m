
function [Y,Do] = createDataDigits(D,N,ind,sigma,active)


rand('twister',23908123);
randn('state',12345678);


% ind = 1:10;
Do = [];
Dgroups = [];
groups = [];
X = [];

for i=1:length(ind)
    
    Do = [Do D{ind(i)}];
    
    idx = find(data.label == ind(i) -1);
    Dgroups = [Dgroups i*ones(1,size(D{ind(i)},2))];
    X = [X data.X(:,idx)];
    groups = [groups i*ones(1,length(idx))];
    
end


Y = zeros(size(X,1),N);

for i=1:N
    
    ii = randperm(length(active));
    
    gt = sort(active(ii(1:2)));
    
    g1 = find(groups == gt(1));
    g2 = find(groups == gt(2));

    
    ig1 = randperm(length(g1));
    ig2 = randperm(length(g2));
    
    It = X(:,g1(ig1(1))) + X(:,g2(ig2(1)));
    It = It/norm(It);
    Y(:,i) = It;

end
