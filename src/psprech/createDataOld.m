
rand('twister',23908123);
randn('state',56411578);

% sigma = 0.02;

n = 30;
sizeD = 60;
N = 500;
k=3;

D1 = randn(n,sizeD);
D2 = randn(n,sizeD);

for i=1:sizeD
    D1(:,i)=D1(:,i)/norm(D1(:,i));
end

for i=1:sizeD
    D2(:,i)=D2(:,i)/norm(D2(:,i));
end

D{1}=D1;
D{2}=D2;

X1 = zeros(n,N);
X2 = zeros(n,N);
R = randn(n,N);

for j=1:N
    
    % signal 1
    idx = randperm(sizeD);
    d = zeros(1,sizeD);
    d(idx(1:k)) = randn(1,k);             
    X1(:,j) = D1*d';
    X1(:,j) = X1(:,j)/norm(X1(:,j));
    
    % signal 2
    idx = randperm(sizeD);
    d = zeros(1,sizeD);
    d(idx(1:k)) = randn(1,k);
    X2(:,j) = D2*d';
    X2(:,j) = X2(:,j)/norm(X2(:,j));
    
    % noise
     R(:,j) = R(:,j)/norm(R(:,j));
    
end

Y = X1 + X2 + sigma*R;

X{1}=X1;
X{2}=X2;


