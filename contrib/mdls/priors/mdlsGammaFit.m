%
% Estimate parameters of a Gamma distribution.
%
%
%
%
% inputs: 
% x ........ a vector of samples. Must be all non-negative!!
% kappa ... optionally fix kappa and solve only for beta, which
% gives a trivial estimate.
%
%
% kappa ...... the parameter kappa
% beta ...... the parameter beta
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [kappa,betta]=mdlsGammaFit(X,k0)
[m,n]=size(X);
if (m > 1) && (n == 1)
    X=X';
    [m,n]=size(X);
end

if nargin >= 2
    kappa=k0*ones(m,1);
    betta=(n*kappa)./sum(X,2);
else
    kappa=zeros(m,1);
    for i=1:m
        a=log(1/n*sum(X(i,:)));
        b=1/n*sum(log(X(i,:)));
        s=a-b;
        %
        % N-R to obtain Kappa. Taken from Wikipedia (Gamma Function)
        % initial guess
        %
        k=(3-s+sqrt((s-3)^2-24*s))/(12*s);
        %
        % explicit Newton step
        % psi(k) = digamma(k), psi(1,k) = trigamma(k)
        j=1;
        while j < 10
            k0=k;
            k=k-(log(k)-psi(k)-s)/(1/k-psi(1,k));
            dif=abs(k-k0);
            if dif < 1e-10
                break;
            end
            j = j + 1;
        end
        kappa(i)=k;
        betta(i)=((kappa(i))*n)/sum(X(i,:));
    end
end

