% function A=mdlsStagewiseCSC(X,D,qstep)
%
% This function implements a non-greedy approximate minimum-codelength 
% sparse coding of signals.
% Signals X are decomposed as D*A + R, where A in turn is represented as
% the product of three independent random variables: Z, a binary vector
% representing the support of A, S, the sign of A, and V, the absolute
% values.
% Each of the four components (R,Z,S,V) are encoded using universal codes
% that can be selected by the user
% (with the exception of S, which is always assumed Bernoulli(1/2)).
%
% INPUT
%
% Required:
%
% X ............... Data sample
% D ............... Dictionary
% qR .............. Quantization step for residual
% qV .............. Quantization step for coefficient magnitude
%
% Optional: (set any to [],  or simply do not specify, to use default value)
%
% (note: the logfun functions are of the generic form:  P =
%   Px_logfun(x (column vector), params (array), q (scalar) )
%
% PR_logfun ........ Function (name or handle) for -log P(R). 
%                    Defaults to universal Gaussian 
%                    (see gMDL model in Hansen and Yu 2001)
%
% PR_params ........ Array of parameters for P(R)
%                    Defaults to empty array.
%
% PZ_logfun ........ Function (name or handle) for -log P(support)
%                    Defaults to universal model for Bernoulli
%                    (see some paper)
%
% PZ_params ........ Array of parameters for -log P(support). 
%                    Defaults to empty array.
%
% PV_logfun ........ Function (name or handle) for -log P(coefficients)
%                    Defaults to universal Mixture of Exponentials  (MOE)
%                    (see Ramirez and Sapiro 2010)
%
% PV_params ........ Array of parameters for -log P(A)
%                    Defaults to empty array.
%
% OUTPUT
%
% A ........... Coefficients
% 
function [A,RPATH]=mdlsStagewiseCSC(X,D,qstep)
%
%-----------------------------------------------------
% Parse input
%-----------------------------------------------------
%
if ~exist('PR_logfun','var')
    PR_logfun = [];
end
if isempty(PR_logfun)
    PR_logfun = @quantized_iid_gaussian_mixture;
end
if ~exist('PR_params','var')
    % for universal, this is computed by ML and encoded at precision sqrt(M) 
    PR_params = []; 
end

if ~exist('PZ_logfun','var')
    PZ_logfun = [];
end
if isempty(PZ_logfun)
    PZ_logfun = @iid_bernoulli_mixture;
end
if ~exist('PZ_params','var')
    % for universal, no need to specify since we are using parameter-free
    % Jeffreys prior
    PZ_params = [];
end

if ~exist('PV_logfun','var')
    PV_logfun = [];
end
if isempty(PV_logfun) 
   PV_logfun = @quantized_iid_laplacian_mixture;
end
if ~exist('PV_params','var')
    % for universal, this is computed by ML and encoded at precision sqrt(M) 
    PV_params = []; 
end
%
%-----------------------------------------------------
% Initialization
%-----------------------------------------------------
%
[M,K] = size(D);
   N  = size(X,2);
DTD = D'*D;
DTX = D'*X;
DTR = DTX;
sDTR = sign(DTR);

A = zeros(K,1);
R = X; % initial residual
LR = PR_logfun(R,PR_params,qR);
LZ = PZ_logfun(R,PZ_params);
LV = PV_logfun(A,PR_params,qV); 
LS  = 0; % number of bit signs = size of active set.
L = LR + LZ + LV + LS;
RPATH=zeros(K,500);
%
%-----------------------------------------------------
% Model selection loop
%-----------------------------------------------------
%
% P(R) is monotonically decreasing in R so that
% at the beginning the first step goes in the direction
% of the maximally-correlated atom.
%
t = 0;
while L0 > L
    %
    % extremely inefficient initial implementation
    %
    A0 = A;
    R0 = R;
    sel_k = 1;
    L = L0;
    for k=1:K
        %
        % step i k-th atom's direction
        %
        stepk = qV*sDTR(k);
        R = R0 - stepk*D(:,k); 
        A(k) = A0(k) + stepk;
        %
        % re-evaluate codelength
        %
        LR = PR_logfun(R,PR_params,qR);
        Z  = double(A~=0);
        LZ = PZ_logfun(Z,PZ_params);
        LS = nnz(Z);             % number of sign  bits = size of active set.
        V  = abs(A);
        LV = PV_logfun(V(Z),PV_params,qV); % does not work if V is not IID!!
        Lk = LR + LZ + LV + LS;
        if Lk < L
            sel_k = k;
            L = Lk;
        end
    end
    stepk = qV*sDTR(sel_k);
    R = R0 -    stepk*D(:,sel_k); 
    A = A +     stepk*sDTR(k);
    DTR = DTR - stepk*sDTR(k)*DTD(:,sel_k);
    t = t + 1;
    RPATH(:,t) = A;
end
