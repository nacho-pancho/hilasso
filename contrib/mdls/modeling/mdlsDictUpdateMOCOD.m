%
% function D = MOCOD(data,D,A,params)
%
% Name: MOCOD
%
% Category: core function
%
% Description: update a dictionary using the MOCOD method
%
% Input:
% data ..... data structure
% D ........ dictionary (n x K)
% A ........ coefficients (K x N)
% params ... parameters
%
% Output:
% D ........ updated dictionary
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Idn$
%
function D = mdlsDictUpdateMOCOD(D,XAt,AAt,params)

%
% Use the sparse version of the matrix of coefficients
%
    if nargin == 0
        params = struct();
        params.mu = 0;
        params.eta = 0;
        params.max_iter = 10;
        params.min_change = 1e-3;
        params.debug = false;
        D = params;
        return;
    end
    mu = params.mu;
    eta = params.eta;
    [n,K] = size(D);
    Idn = eye(n); % nxn identity matrix
    IdK = eye(K); % KxK identity matrix
    
    D0 = D;
    iter = 1;
    dif = realmax;
    H = AAt;
    G = XAt;
    if (mu > 0)  || (eta > 0) % do fixed point iteration
        while iter <= params.max_iter && (dif > params.min_change)
            if params.eta > 0
                H = AAt +2*mu*(D0'*D0) + 2*eta*diag(sum(D0.*D0));
                G = (XAt+2*(mu+eta)*C*D0);
            else
                H = inv(AAt + 2*mu*(D0'*D0));
                G = XAt+2*mu*C*D0;
            end
            L = modchol(H);
            Y = L \ G';
            D = mdlsDictNormalize((L' \ Y)');
            dif = norm(D-D0,'fro')/norm(D0,'fro');
            if params.debug
                fprintf('||D-D0||=%f\n',dif);
            end
            D0=D;
            iter = iter + 1;
        end
    else
        % D = XAt*inv(AAt); using modified Cholesky factorization
        L = modchol(H);
        Y = L \ G';
        D = mdlsDictNormalize((L' \ Y)');
    end
end % function
