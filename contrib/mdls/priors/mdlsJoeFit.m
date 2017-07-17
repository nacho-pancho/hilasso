%
% Fit parameters for the Jeffreys of Exponential (Joe) prior
%
%                1     1
% Joe(x;k,b)= ------- --- [ e^(-ax) - e^(-bx) ]
%             ln(b/a)  x 
% 
% inputs: 
% x ......... training samples
% a0 ........ initial value for a
% b0 ........ initial value for b
%
% outputs:
% a,b ........ fitted distribution parameters

% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [a,b]=mdlsJoeFit(x,a0,b0,sym)
    if size(x,2) == 1 
        x = x';
    end
    if ~exist('do_newton', 'var')
        do_newton = false;
    end
    if ~exist('b0','var')         
        b0 = -1;
    end
    if ~exist('a0','var') 
        a0 = -1;
    end
    if ~exist('sym','var')
        sym = ~isempty(find(x<0));
    end
    if b0 <= 0
        b0 = 100;
    end
    if a0 <= 0
        a0 = .1;
    end
    if sym
        x = abs(x);
    end
    % we sort x for better numerical implementation
    x=fliplr(sort(x,2));
    a = a0;
    b = b0;
    av = zeros(1,20);
    bv = zeros(1,20);
    [m,n] = size(x);
    finished = false;
    for i=1:20
        fprintf('.');
        av(i) = a;
        bv(i) = b;
        
        iSabd = 1./(exp(-a*x) - exp(-b*x)+eps); % zeros may occur
        Sa = exp(-a*x).*x;
        Sb = exp(-b*x).*x;
        Sabn = exp(-(a+b)*x).*(x.^2);
        g = [  -n/log(b/a)/a + sum(Sa.*iSabd); ...
                n/log(b/a)/b - sum(Sb.*iSabd)];
        Sabab = sum(Sabn.*iSabd.^2);
        Hab = n/(a*b*log(b/a).^2) + Sabab;
        Hba = Hab;
        Haa = n*(log(b/a)-1)/(log(b/a)*a)^2 +  Sabab;
        Hbb = -n*(1+log(b/a))/(log(b/a)*b)^2 - Sabab;
        H = [ Haa , Hab ; Hba Hbb ];
        L = modchol(H);
        y = L\g;
        d = -L'\y;
        if find(isnan(d))
            error('Numeric error');
        end
        g'*d;
        ab = [a;b];
        a = a+d(1);
        if a <= 0
            fprintf('reprojecting a=%g\n',a);
            a = eps;
        end
        b = b+d(2);
        if b<a
            fprintf('reprojecting b=%g\n',b);
            b = 1.1*a; 
        end
        ab2=[a; b];
        %mdlsArmijo(@joe_loglik, ab, g, d);
        if norm(ab2-ab)/norm(ab) < 1e-3
            break;
        end
    end
    fprintf('\n');
    av(i) = a;
    bv(i) = b;
    av=av(1:i); bv=bv(1:i);
    if 0
        avi= min(av):(max(av)-min(av))/10:max(av);
        bvi = min(bv):(max(bv)-min(bv))/10:max(bv);
        me = zeros(length(avi),length(bvi));
        for i=1:length(avi)
            for j=1:length(bvi)
                me(i,j) = joe_loglik(x,avi(i),bvi(j));
            end
        end
        mdlsFigure('LogLik'); mesh(bvi',avi',me); hold on;
        for i=1:length(av)
            ll(i) = joe_loglik(x,av(i),bv(i));
        end
        plot3(bv,av,ll,'o-','LineWidth',2,'Color',[1 0.5 0]); hold off;
    end
end
 
function f=joe_loglik(x,a,b)
    %x = sort(x,'descend');
    n =  size(x,2);
    f = n*log(log(b/a))- sum( log((exp(-a*x)-exp(-b*x))./x) );
end