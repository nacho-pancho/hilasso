%
% solve the Hierarchical Lasso problem using SpaRSA + ADMOM
%
% y ........ m x n data vector
% D ........ m x p dictionary
% x0
% 
% ........ partition of the index set {1,2,...,p}
%
function [x_hilasso,obj,times,seq] = HiLasso(y,D,x0,groups,lambda1,lambda2,tol,max_iter, stop_crit)

%
% in this script, lambda2 is L1 penalty and lambda1 is L2 penalty
% but on calling, it is reversed, so we swap them.
%
aux = lambda2;
lambda2 = lambda1;
lambda1 = aux;

if ~exist('max_iter','var')
    max_iter = 200;
end
if ~exist('tol','var')
    tol = 1e-6;
end


% compatible with SpaRSA criterion id's
STOP_BY_COST = 1;
STOP_BY_ARG  = 5;

hR = @(x) D*x;
hRt = @(x) D'*x;


if ~exist('stop_crit','var')
    stop_crit = STOP_BY_COST;
end

if lambda1==0
  lambda1 = eps;
end
%keyboard
%psi = @(x,tau) group_lam1l2lam2l1(x,tau,lambda1/lambda2,groups);
fac12 = lambda1/lambda2; % L2/L1 penalty
psi = @(x,tau) hilasso_subproblem(x,groups,tau,fac12,500,1e-6,10);
phi = @(x) group_lam1l2lam2l1norm(x,fac12,groups);

[x_hilasso,x_debias_SpaRSA,obj,...
    times,debias,...
    mse,taus,seq]= ...
    SpaRSA(y,hR,lambda2,... % lambda2 is L1 penalty
    'Psi',psi,...
    'Phi',phi,...
    'Monotone',1,...
    'Debias',0,...
    'AT',hRt,...
    'StopCriterion',stop_crit,...
    'ToleranceA',tol, ...
    'MaxiterA',max_iter,...
    'Verbose',0,...
    'Initialization',x0);


% x_lasso = SpaRSA(y,hR,lambda2,...
%     'Monotone',1,...
%     'Debias',0,...
%     'AT',hRt,...
%     'StopCriterion',stop_crit,...
%     'ToleranceA',tol, ...
%     'MaxiterA',max_iter,...
%     'Verbose',0,...
%     'Initialization',x0);

