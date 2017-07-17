%
% Estimate parameters of a Gamma-Laplace (GAMLAP) distribution.
%
% This is the result of mixing a Laplacian distribution using a
% Gamma(lambda;k,b) prior for the Laplacian parameter lambda.
% This results in a distribution with the two "hyper-parameters"
% k,b, which correspond to those of the mixing distribution
% Gamma.
%
% This function estimates the parameters b and k  using the
% method of moments applied to the one-sided distribution of the absolute
% values of x. An ML estimation of these does not have a
% closed form (and I don't figure out yet how to do it
% numerically).
%
% Optionally, the value of k can be given in advance. This value controls
% the *shape* of the mixing Gamma distribution and one may fix it
% to define a desired shape of the mixing Gamma. Also, the
% estimation of both parameters is difficult and unreliable due to
% an almost flat direction in the surface of the likelihood function.
%  
% A value of k=3 to 5 is recommended for modeling responses to image patches
% for various reasons: first, it is the most common one. Second,
% the shape of the mixing Gamma resembles closely the empirical
% distribution  of the observed fitted Laplacian parameters.
%
% Finally, the mean of the distribution will also be estimated if
% a third output parameter is present.
%
% inputs: 
% x ........ a vector of samples.
% k ........ optionally fix k and solve only for b. If this value
%            is less than or equal of 0 (invalid), it is ignored
%            and computed.
% method ... method for obtaining the parameter estimates. This is
%            really only needed if k is not specified and the estimated k
%            is below 2.
%            Can be:
%            'moments' ... use method of moments. Only gives valid
%                          results for k >= 3. 
%            'mle' ....... maximum likelihood. This is more expensive.
%           
% outputs:
%
% k ...... the parameter k
% b ...... the parameter b
% u ...... the mean of X. If this parameter is requested, it is
%          computed, otherwise it is assumed 0.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [k,b,u]=mdlsMoeFit(X,k0,method,sym)
%
% estimation of mean is requested.
%
if nargout == 3
    u = mean(X);
    X = X - u;
end
if ~exist('sym','var')
    sym = ~isempty(find(X<0));
end
if ~exist('method','var')
    method = 'moments';
end
if ~exist('k0','var')
    k0=-1;
end
%
% symmetric evaluation
%
if sym
    X = abs(X);
end

if size(X,2)==1
    X=X';
end
[m,n]=size(X);
%
% initial guess
% 
mu1=sum(X,2)/n;
%
% if kappa is given, betta is very simple and robust
%
only_fit_betta = (k0 > 0);
if only_fit_betta
    b0=(k0-1)*mu1;
    k0=k0*ones(1,m);
else
    mu2=sum(X.*X,2)/n;
    k0=2*(mu2-mu1.^2)./(mu2-2*mu1.^2);
    b0=(k0-1).*mu1;
end
if isequal(method,'moments')
    k=k0; b=b0; return;
end

converged = 0;
opt_error = '';
armijo_stats = zeros(m,1);
if ~only_fit_betta
    %
    % fit both betta and kappa
    %
    for i=1:m
        ki=k0(i);
        bi=b0(i);
        tol = 1e-8;
        dif = realmax;
        max_iters = 20;
        j = 1;
        x = X(i,:);
        while j < max_iters
            %
            % since we want to minimize -f
            % we descend in the direction of -g
            % therefore the minus sign in the definition of g and H
            % gradient
            %
            fprintf('k=%f b=%f\n',ki,bi);
            xb = x+bi;
            S1 = sum(log(xb));
            S2 = sum(1./xb);        
            gk = n/ki + n*log(bi) - S1;
            gb = n*ki/bi -(ki+1)*S2;
            g  = -[gk; gb];
            %
            % Hessian
            %
            Hkk = -n/ki^2;
            Hbb = -n*ki/bi^2+(ki+1)*sum(xb.^(-2));
            Hkb =  n/bi - S2;
            Hbk = Hkb;
            H = -[Hkk, Hkb; Hbk, Hbb]; 
            eig(H)
            L = modchol(H);
            y = L\g;
            d = -(L'\y);
            ang = dot(d,g);
            if ang >= 0
                opt_error = 'Not a descent direction';
                break;
            end
            dif = norm(g);
            if dif < tol
                converged = converged+1;
                break;
            end
            % Step selected using Armijo rule
            armijo_sigma = 0.1; 
            armijo_beta = 0.25; 
            armijo_ok = 0;
            s = 1;
            k0 = ki; b0 = bi;
            f0 = -(n*log(k0)+n*k0*log(b0)-(k0+1)*sum(log(x+b0)));
            armijo_delta = -s*armijo_sigma*ang
            armijo_m = 1;
            while ~armijo_ok && (armijo_m < 20) 
                ki = k0 + s*d(1);
                bi = b0 + s*d(2);
                %
                % project onto feasible set: b>0, k>=2
                %
                bi = max(bi,1e-5);
                ki = max(ki,2);
                f = -(n*log(ki)+n*ki*log(bi)-(ki+1)*sum(log(x+bi)));
                df = f0 - f
                if df >= armijo_delta
                    armijo_ok = 1;
                else
                    armijo_delta = armijo_delta*armijo_beta;
                    s = s*armijo_beta;
                    armijo_m = armijo_m + 1;
                end
            end
            armijo_stats(i) = armijo_m;
            if ki < 2
                opt_error = 'Invalid solution (kappa<2).';
                break;
            end
            if  (f0 - f) < 1e-5
                converged = converged + 1;
                break;
            end
            fprintf('%d: f=%f\tk=%f\tb=%f\tdif=%f\ts=%f\n',j,f,ki,bi, ...
                    dif,s);
            j = j + 1;
        end
        if length(opt_error) > 1
            warning(sprintf(['Optimization error: %s. Failing back ' ...
                             'to moment based estimations.'],opt_error));
        else
            k(i)=ki;
            b(i)=bi;
        end
    end  % for each row of X
else
    %
    % only fit betta
    %
    fprintf('only fit betta\n');
    k=k0;
    for i=1:m
        bi=b0(i);
        ki=k0(i);
        tol = 1e-8;
        dif = realmax;
        max_iters = 10;
        j = 1;
        x = X(i,:);
        while j < max_iters
            %
            % since we want to minimize -f
            % we descend in the direction of -g
            % therefore the minus sign in the definition of g and H
            % gradient
            %
            fprintf('k=%f b=%f\n',ki,bi);
            xb = x+bi;
            S1 = sum(log(xb));
            S2 = sum(1./xb);        
            gb = n*ki/bi -(ki+1)*S2;
            g  = -[gb];
            %
            % Hessian
            %
            Hbb = -n*ki/bi^2+(ki+1)*sum(xb.^(-2));
            H = -[Hbb]; 
            d = -g/H;
            ang = d*g;
            if ang >= 0
                opt_error = 'Not a descent direction';
                break;
            end
            dif = norm(g);
            if dif < tol
                converged = converged+1;
                break;
            end
            % Step selected using Armijo rule
            armijo_sigma = 0.1; 
            armijo_beta = 0.25; 
            armijo_ok = 0;
            s = 1;
            k0 = ki; b0 = bi;
            f0 = -(n*log(k0)+n*k0*log(b0)-(k0+1)*sum(log(x+b0)));
            armijo_delta = -s*armijo_sigma*ang
            armijo_m = 1;
            while ~armijo_ok && (armijo_m < 20) 
                bi = b0 + s*d;
                %
                % project onto feasible set: b>0, k>=2
                %
                bi = max(bi,1e-5);
                f = -(n*log(ki)+n*ki*log(bi)-(ki+1)*sum(log(x+bi)));
                df = f0 - f
                if df >= armijo_delta
                    armijo_ok = 1;
                else
                    armijo_delta = armijo_delta*armijo_beta;
                    s = s*armijo_beta;
                    armijo_m = armijo_m + 1;
                end
            end
            armijo_stats(i) = armijo_m;
            if  (f0 - f) < 1e-5
                converged = converged + 1;
                break;
            end
            fprintf('%d: f=%f\tk=%f\tb=%f\tdif=%f\ts=%f\n',j,f,ki,bi, ...
                    dif,s);
            j = j + 1;
        end
        b(i) = bi;
    end  % for each row of X
end

fprintf('Average Armijo iterations: %5.2f\n',mean(armijo_stats));
