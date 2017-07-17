
function [Y,D,X,R,Ao] = createData(n,N,k,numD,sizeD,sigma,Nactive,ort)


if ~exist('Nactive','var')
  Nactive = numD;
end

if ~exist('ort','var')
  ort = 0;
end


% Select active groups
idxg = randperm(numD);
active = zeros(1,numD);
active(idxg(1:Nactive)) = 1;


% Initialize
for i=1:numD
    X{i} = zeros(n,N);
    D{i} = randn(n,sizeD); 
    A{i} = zeros(sizeD,N); 
    
    if ort
        D{i} = orth(D{i});
    else
        for j=1:sizeD
            D{i}(:,j)=D{i}(:,j)/norm(D{i}(:,j));
        end
    end
end

R = randn(n,N);
for j=1:N
    
    % signal 1
    for i=1:numD
        
        if active(i)
            idx = randperm(sizeD);
            d = zeros(1,sizeD);
            d(idx(1:k)) = rand(1,k)-0.5;
            X{i}(:,j) = D{i}*d';
            A{i}(:,j) = d'/norm(X{i}(:,j));
            X{i}(:,j) = X{i}(:,j)/norm(X{i}(:,j));
        end
        
    end
    
    R(:,j) = R(:,j)/norm(R(:,j));
   
end

Y = sigma*R;
for i=1:numD
    Y = Y+X{i};
end

 Ao = [];
for i=1:numD
    A{i} = sparse(A{i});
    Ao = [Ao;A{i}];
end


