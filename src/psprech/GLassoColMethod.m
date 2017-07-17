
        
        
function [Xr,A,v,Xo,Ao,vo] = GLassoColMethod(Y,D,A0,lambda2,tol,max_iter)

if ~exist('tol','var')
    tol = 1e-5;
end

if ~exist('max_iter','var')
    max_iter = 100;
end

%if isempty(H)
    H = ones(size(Y));
    %end

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
    Do = [Do D{i}];
    groups = [groups i*ones(1,sizeD)];
end

% Solve collaborative Hilasso
if isempty(A0)
  A = zeros(size(Do,2),N);
else
  A = A0;
end  
[A,obj,times] = GLassoCollaborative(Y,Do,A,groups,lambda2, ...
                                    tol,max_iter, stop_criterion);

Ao = A;
reg = 1e-2;
if lambda2~=0
    fprintf('Back-projecting solution using OLS\n');
    for i=1:N        
        [Xd,A(:,i)] = compute_ols(Y(H(:,i)==1,i),Do(H(:,i)==1,:),A(:,i),reg);
    end
end

% Recovered signals
for i=1:numSig
    idx = find(groups == i);
     Xr{i} = D{i}*A(idx,:);
     v(i) = sum(sum((A(idx,:).^2)));
end

% Recovered signals
for i=1:numSig
    idx = find(groups == i);
     Xo{i} = D{i}*Ao(idx,:);
     vo(i) = sum(sum((Ao(idx,:).^2)));
end
