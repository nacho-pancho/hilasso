%
% x ........ data
% tau ...... L1 penalty
% fac12 .... ratio L2 penalty / L1 penalty
% groups ... group partition
%
function y = group_lam1l2lam2l1(x,tau,fac12,groups)
num_groups = max(groups);
if min(groups)~=1
   error(['Wrong group structure vector'])
end
y = zeros(size(x));
for i=1:num_groups
    
    thisgroup = find(groups == i);
    f = x(thisgroup);
    n = length(thisgroup);
    
    % first check if 0 is a solution of L2-regularized and L1-regularized problems
    %
    if  (norm(vector_soft(f,tau*fac12))>0) && (norm(soft(f,tau)) > 0)
        w = minAdmom(f,tau,fac12,1000,1e-4);
    else
        w = zeros(size(f));
    end
    
    y(thisgroup) = w;
end

function w = minAdmom(x,tau,fac12,max_iter,tol,c)

% Set optimization parameters
if ~exist('c','var')
    c = 1; % works good here
end

if ~exist('max_iter','var')
    max_iter = 500;
end
if ~exist('tol','var')
    tol = 1e-8;
end

%st intial points
p    = ones(size(x));
b    = rand(size(x));
beta = rand(size(x));

for iter=1:max_iter
    b0 = b;
    beta0 = beta;
    % update b
    % call the soft thresholding changing variable
    b = soft(x+c*beta-p,tau);
    b = b*(1/(c+1));
    
    % update beta
    % call the vector soft thresholding changing variable
    beta = vector_soft(p+c*b,tau*fac12);
    beta = beta*(1/c);
    
    % update lagrange multipliers
    %p = p + c*(b - beta);
    p = p + c*(b - beta);
    
    %dif = norm(b-beta)/norm(b+eps);
    dif = norm(beta-beta0)/(norm(beta)+eps);
    if dif < tol
        break;
    end
end
w = b;

