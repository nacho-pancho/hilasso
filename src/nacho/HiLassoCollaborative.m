%
% solve the Hierarchical Lasso problem using BCD
%
% X .......... m x n data vector
% D .......... m x p dictionary
% G .......... partition of the index set {1,2,...,p}
% lambda1 .... L1 penalty
% Flambda2 .... Base group penalty. It will be multiplied by sqrt(num
%               samples)*sqrt(dic size).
% tol ........ tolreance
% max_iter ... maximum iterations (defaults to 500)
% stop_crit .. stop criterion 1:cost fun, 5:argument
% H .......... for missing data, indicates missing samples with an 1
% c .......... ADMOM constant (advanced usage)
%
function [A,obj,times,seq] = HiLassoCollaborative(X,D,A0,groups,lambda1,Flambda2,tol, ...
    max_iter,stop_crit,H,c)

if ~exist('H','var')
    H = ones(size(X));
end

if ~exist('c','var')
    c = 0.2;
end

if ~exist('max_iter','var')
    max_iter = 500;
end
if ~exist('tol','var')
    tol = 1e-6;
end
save_seq = (nargout > 3);
%
% compatible with SpaRSA criterion id's
%
STOP_BY_COST = 1;
STOP_BY_ARG  = 5;
STOP_BY_ACT_SET = 6;

if ~exist('stop_crit','var')
    stop_crit = STOP_BY_COST;
end

times = zeros(1,max_iter);
obj   = zeros(1,max_iter);
seq   = zeros(length(A0),max_iter);

[n,m] =  size(X);
k = size(D,2);
G = max(groups);

P    = zeros(m,k);
if isempty(A0)
    B    = zeros(k,m);
    Beta = zeros(k,m);
else
    B = A0;
    Beta = A0;
end

ng = zeros(1,G);

for i=1:G;
    g{i} = (groups == i);
    ng(i) = sum(g{i});
    Bi = B(g{i},:);
    baux(i) = norm(Bi,'fro');
end
lambda2 = (Flambda2*sqrt(m)).*sqrt(ng);
for t=1:max_iter
    times(t) = cputime();
    r = X-D*B;
    r = r(:);
    obj(t)   = 0.5*r'*r + baux*lambda2' + lambda1*sum(abs(B(:)));
    lag(t)   = obj(t) +  0.5*c*sum(sum((B-Beta).^2));
    %fprintf('t=%3d : Lc=%f\t',t,obj(t));
    fprintf('t=%3d: f=%6f  Lc=%6f  ',t,obj(t),lag(t));
    
    B0 = B;
    Beta0 = Beta;
    %
    % *** update b
    %
    % need to solve a modified LASSO-like problem for each colum of X
    %
    % p must be a column vector in lasso_alt
    for i=1:m
        B(:,i) = lasso_alt(X( H(:,i)==1 ,i),D(H(:,i)==1,:),lambda1,P(i,:)',c,Beta(:,i),B0(:,i));
    end
    %
    % *** update beta
    %
    % call the vector soft thresholding for each group
    
    for i=1:G;
        
        
        Bi = B(g{i},:);
        bi = Bi(:);
        baux(i) = norm(bi);
        PiT = P(:,g{i})';
        pi = PiT(:);
        
        betai = vector_soft(pi+c*bi,lambda2(i));
        betai = betai*(1/c);
        
        Beta(g{i},:) = reshape(betai,ng(i),m);
        
    end
    %
    % update lagrange multipliers
    %
    P0 = P;
    P = P + c*B' - c*Beta';


    % Display values
%     difb = norm(Beta - B,'fro')/(norm(B,'fro')+eps);
    h = mdlsHammingDistance(B,B0);
    fprintf('HAM = %7.2f\t',h);
    gg = activeGroups(Beta,g);
    fprintf('ActGroups = %s\t',gg);
    difb = norm(P - P0,'fro')/(norm(P0,'fro')+eps);
    fprintf('||da||/||a|| = %f\n',difb);
    if t > 1
        switch (stop_crit)
            case STOP_BY_ARG,
%                 if abs(0.5*c*sum(sum((B-Beta).^2))) < tol
%                     fprintf('Converged (in value).\n');
%                     break;
%                 end
                if difb < tol
                    fprintf('Converged (in arg).\n');
                    break;
                end
            case STOP_BY_COST,
                if abs(obj(t)-obj(t-1))/abs(obj(t)) < tol
                    fprintf('Converged (in value).\n');
                    break;
                end
            case STOP_BY_ACT_SET,
              if h  < tol
                  fprintf('Converged (in act set).\n');
                  break;
              end
        end
    end
end
if t >= max_iter
    fprintf('Reached maximum iterations!!.\n')
end
%
%
%
if t < max_iter
    times = times(1:t);
    obj   = obj(1:t);
    seq   = seq(:,1:t);
end

%baux
A = abs(B).*sign(Beta);

end




function b = lasso_alt(x,D,lambda,p,c,beta,b0)

b = SpaRSAalt(x,D,lambda,p,c,beta,...
    'Monotone',1,...
    'Debias',0,...
    'StopCriterion',1,...
    'ToleranceA',0.0001, ...
    'ToleranceD',0.000001, ...
    'MaxiterA',200,...
    'Verbose',0, ...
    'Initialization',b0);
end


function gg = activeGroups(B,g)

gg = zeros(1,length(g));
for i=1:length(g)
    Bi = B(g{i},:);
    bi = Bi(:);
    gg(i) = norm(bi);
end
gg = mdlsShowGroupActivity(gg);
end



