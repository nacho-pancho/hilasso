%
% Solve the collaborative Group Lasso problem using BCD
%
% X .......... m x n data vector
% D .......... m x p dictionary
% A0 ......... initial solution
% G .......... partition of the index set {1,2,...,p}
% lambda2 .... Group lasso penalty
% tol ........ tolereance
% max_iter ... maximum iterations (defaults to 500)
% stop_crit .. stop criterion 1:cost fun, 5:argument
% H .......... ????
%
function [A,obj,times] = ...
    GLassoCollaborative(X,D,A0,G,lambda2,tol, ...
    max_iter,stop_crit)


if ~exist('max_iter','var')
    max_iter = 500;
end
if ~exist('tol','var')
    tol = 1e-5;
end
if ~exist('stop_crit','var')
    stop_crit = 1;
end

[M,N] = size(X);
debias = 0;

hD = @(z) D*z;
hDt = @(z) D'*z;

psi = @(z,tau) col_group_vector_soft(z,tau,G);
phi = @(z) col_group_l2norm(z,G);

[A, Adeb, obj,...
      times, debias_start_SpaRSA, mse] = ...
    SpaRSA(...
        X,...
        hD,...
        lambda2,...
        'Psi',psi,...
        'Phi',phi,...
        'Monotone',1,...
        'Debias',debias,...
        'AT',hDt,... 
        'Initialization',A0,...
        'StopCriterion',stop_crit,...
        'ToleranceA',tol, ...
        'ToleranceD',tol, ...
        'MaxiterA',max_iter);

times = times(end);
