
function [Y,X] = createDataDigits(data,N,active,k,sigma)

if ~exist('sigma','var')
    sigma = 0;
end

if ~exist('k','var')
    k = length(active);
end

g = length(unique(data.label));
n = size(data.X,1);

Y = zeros(n,N);
for i=1:g
    X{i} = zeros(n,N);
end

for i=1:N
    
    ii = randperm(length(active));
    
    h=1;
    It = zeros(n,1);
    for j=1:k
        g1 = find(data.label == active(ii(h)) );
        ig = randperm(length(g1));
        It = It + data.X(:,g1(ig(1)));
        X{active(ii(h))+1}(:,i) = data.X(:,g1(ig(1)));
        h=h+1;
    end
    
%     It = It/norm(It);
    Y(:,i) = It;

end

Y = Y + sigma*randn(n,N);