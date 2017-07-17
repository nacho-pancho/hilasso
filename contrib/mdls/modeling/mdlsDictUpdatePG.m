%
% function D = mdlsDictUpdatePG(D0,AAt,XAt,D2,params)
%
% Category: modeling
%
% Description: Update a dictionary using 
% approach, optionally scaling by a diagonalized Hessian or using a more
% expensive Two-Metrics update (full Hessian).
%
%  The cost function is defined as:
% 
% min_D f(D) = ||X-DA||_F^2 + \lambda||A||_1 + \mu||D'D-I||_F^2 +
% \xmu||D_2'D|| 
% 
% s.t. ||D_k|| = 1 k = 1 , ... , K
%
% This cost function accomodates the possibility of optimizing
% for internal incoherence in the dictionary D, or cross-incoherence
% with another dictionary (or set of dictionaries).
%
% The dictionary is updated on a per-column (atom) basis. The gradient
% w.r.t. a given column D_k is
%
% g_k(D) = 2*(D*A-X)*(A^k)' + 4*mu1*(D*D'-I)*D_k + 2*mu2*(D_2*D_2')*D_k
% 
% where A^k is the k-th row of A. The Hessian is
%
% H_k(D) = 2*A^k(A^k)'*I_m + 4*mu*[ D_k'*D_k*I_m + D_k*D_k' ] + 
%        2*xmu*D_2*D_2'
%
% Both the Gradient and Hessian are determined by D0, A*A' and X*A'
% and since the last two are generally much smaller than A and X alone,
% we take these as arguments.
%
% Input:
% D0 ....... current dictionary
% AAt ...... A*A' empirical correlation matrix for the representation
%            coefficients
% XAt ...... X*A' empirical cross-correlation between input and
%            representation
% D_2 ...... If D_2 = [] or mu2 = 0, the 'cross-incoherence term' is not
%            taken into account.
%
% params ... update parameters
%     mu ...... self-coherence penalty
%     xmu ...... cross-coherence penalty
%     amuco ..... average mutual coherence. This is the off-diagonal
%                 value to which we compare against in f. Defaults to 0.
%     max_iter    number of iterations of the algorithm
%     reg_delta   regularization delta value for inverting non-definite
%     positive Hessians.
%     armijo_a .. see mdlsArmijo help
%     armijo_b .. see mdlsArmijo help
%     armijo_s .. see mdlsArmijo help
%     armijo_m .. see mdlsArmijo help
%     scaling ... 
%        'none' ....... no scaling
%        'diag' ....... diagonalized Hessian
%        'full'   ....... two-metrics
%
% with no input, the argument returns the default parameters.
%
% Output:
% D ........ updated dictionary
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Idn$
%
function [D,varargout] = mdlsDictUpdatePG(D0,AAt,XAt,D2,params)
    if nargin==0
        params = struct();
        params.mu = 0;
        params.xmu = 0;
        params.max_iter = 1;
        params.amuco = 0;
        params.reg_delta = 1e-20;
        params.armijo_a = 0.1;
        params.armijo_b = 0.5;
        params.armijo_m = 10;
        params.armijo_s = 1;
        params.scaling = 'diag';
        params.debug = 0;
        params.debug_armijo = 0;
        params.positive = false;
        D = params;
        return;
    end
    D = D0;
    [M,K]=size(D);
    I_M = speye(M);
    I_K = speye(K);
    if params.debug
        fprintf('eig(AAt): min,max=(%g,%g)\n',min(diag(AAt)),max(diag(AAt)));
    end
    armijo_m_stats = zeros(K,1);
    stuck_stats = zeros(K,1);
    dCost = 0;
    % special treatment if no incoherence is required
    %
    % The algorithm is radically different if self-incoherence is sought
    % in the target dictionary, since the problem becomes highly coupled
    % and an update of each atom separately is not possible.
    %
    %
    % case 1 : no incoherence sought: this is plain L2-L1 minimization
    %
    if (params.mu == 0) && (params.xmu == 0)
%        fprintf('No coherence\n');
        for J=1:params.max_iter
            for k=1:K
                g_k = 2*(D*AAt(:,k) - XAt(:,k));
                D_k = D(:,k) - (0.5/max([AAt(k,k) params.reg_delta]))* ...
                      g_k;
                if params.positive
                    D_k = max(D_k,0);
                end
                if norm(D_k) == 0
                    fprintf('Ooops!');
                    D_k = randn(M,1);
                end
                D(:,k) = (1/norm(D_k))*D_k;
            end
            %
            % evaluate cost function
            %
            
        end
    elseif params.mu == 0 && ~isempty(D2)
%        fprintf('Only cross coherence, xmu=%g\n',params.xmu);
        %
        % crossed incoherence with another dictionary D2
        %
        DDt2 = D2*D2';
        for J=1:params.max_iter
            if params.debug
                fprintf('%4d/%4d',J,params.max_iter);
            end
            D0 = D;
            for k=1:K
                if params.debug
                    fprintf('k=%4d/%4d:',k,K);
                end
                D_k = D0(:,k);

                %
                % decesnt direction
                %
                % gradient:
                g_k = 2*(D0*AAt(:,k) - XAt(:,k)) + (2*params.xmu)*(DDt2)* ...
                      D_k;
                %
                % scaling:
                %
                if isequal(params.scaling,'full')
                    a = 2*AAt(k,k);
                    H_k = a*I_M + (2*params.xmu)*DDt2;
                    L_k = modchol(H_k);
                    y = L_k \ -g_k;
                    d_k = L_k' \ y;           
                elseif isequal(params.scaling,'diag')
                    a = max(AAt(k,k),params.reg_delta);
                    H_k = 2*a + (2*params.xmu)*diag(DDt2);
                    d_k = -g_k ./ H_k;
                else
                    d_k = -g_k;
                end
                %
                % step selection: Armijo
                %
                s = params.armijo_s;
                acceptable_decrease = -params.armijo_a*g_k'*d_k;
                if acceptable_decrease < 0
                    error('Not a descent direction');
                end
                dfs = zeros(1,params.armijo_m);
                ads = zeros(1,params.armijo_m);
                D0 = D;
                D_k0 = D0(:,k);
                % df is not computed explicitly but a fast computation
                % is done when only one atom changes.
                %
                %f0 = trace(-2*XAt*D0'+D0*AAt*D0') + ...
                %     params.xmu*trace(D0'*DDt2*D0)
                Ck = -XAt(:,k) + D0*AAt(:,k);
                Bk = D_k0' * DDt2;
                for j=1:params.armijo_m
                    D_k = D_k0 + s*d_k;    % move along d_k
                    if params.positive
                        D_k = max(D_k,0);
                    end
                    D_k = D_k* (1/norm(D_k)); % projection onto unit sphere
                    D(:,k) = D_k;
                    dD_k = (D_k-D_k0); 
                    %
                    % fast computation of df
                    %
                    df = -(2*Ck'*dD_k + AAt(k,k)*sum(dD_k.^2) + ...
                         params.xmu*((2*Bk + dD_k'*DDt2)*dD_k));
                    if params.debug_armijo 
                        fprintf('df=%g\ttarget=%g\n',df,acceptable_decrease);
                    end
                    if  df >= acceptable_decrease
                        break;
                    end
                    % no luck: reduce step and acceptable decrease
                    s = s * params.armijo_b;
                    acceptable_decrease = acceptable_decrease * params.armijo_b;
                end
                if df < 0
                    D(:,k) = D_k0;
                    df = 0;
                end
                dCost = dCost + df;
                if j >= params.armijo_m % no further significant improvement
                    break;
                    stuck_stats(j) = 1;
                end 
                armijo_m_stats(k) = j;
                if params.debug
                    fprintf('\b\b\b\b\b\b\b\b\b');
                end
            end % update of each atom
            if params.debug
                fprintf('\r');
            end
        end % descent iterations

        if params.debug
            fprintf('\n');
            fprintf('Average armijo iterations:%g\n',mean(armijo_m_stats));
            fprintf('Stuck :%g%%\n',100*mean(stuck_stats));
            fprintf('dCost=%g\n',dCost);
        end
    else
%        fprintf('Mutual incoherence, mu=%g xmu=%g\n',params.mu,params.xmu);
        % when mu > 0, 
        % because of the coupling 
        % we need to do the descent simultaneously on all D
        %
        armijo_m_stats=zeros(1,params.max_iter);
        for J=1:params.max_iter
            if params.debug
                fprintf('%4d/%4d',J,params.max_iter);
            end
            DtD = D'* D;
            DDt = D * D';
            %
            % descent direction
            %
            % gradient
            g_k = 2*(D*AAt - XAt);
            g_k = g_k + DDt*((4*params.mu)*D) - (4*params.mu)*D;
            if params.xmu > 0 && ~isempty(D2)
                DDt2 = D2*D2';
                g_k = g_k + (2*params.xmu)*(DDt2)*D;
            end
            %
            % Hessian (diagonalized)!
            % actually, I store it as a matrix of the size of D, and then do d_k = -G_k./H_k
            if isequal(params.scaling,'diag')
                tmp = 2*diag(AAt) + 4*params.mu*diag(DtD);
                tmp2 = (4*params.mu)*D.^2;
                H_k = repmat(tmp',M,1) + tmp2;
                if params.xmu > 0 && ~isempty(D2)
                    H_k = H_k + repmat( (2*params.xmu) * diag(DDt2), 1, K);
                end
                H_k = max( H_k , params.reg_delta );            
                d_k = -g_k ./ H_k;
            else % only other possible choice is no scaling
                d_k = -g_k;
            end
            %
            % step selection: Armijo
            %
            D0 = D;
            s = params.armijo_s;
            acceptable_decrease = -params.armijo_a*g_k'*d_k;
            if acceptable_decrease < 0
                error('Not a descent direction');
            end
            %E(D) = Tr[(X-DA)(X-DA)t] = Tr[XXt -XAtDt -DAXt + DAAtDt] = C + Tr[DA*(DA-X)t]
            %E(D+dD) = C + Tr[(D+dD)A((D+dD)A-X)t] 
            %E(D+dD) - E(D) = Tr(-2*XAt*dD'+2*D*AAt*dD'+dD*AAt*dD')
            f0 = trace(-2*XAt*D0'+D0*AAt*D0') + ...
                 params.mu*sum(sum( (D0'*D0-I_K).^2 ));
            if ~isempty(D2)
                f0 = f0 + params.xmu*sum(sum((D2'*D0).^2));
            end
            for j=1:params.armijo_m
                D = D0 + s*d_k;    % move along d_k
                ND = repmat(sqrt(sum(D.^2)),M,1);
                if params.positive
                    D = max(D,0);
                end
                D = D./ND;
                f = trace(-2*XAt*D'+D*AAt*D') + ...
                    params.mu*sum(sum( (D'*D-I_K).^2 ));
                if ~isempty(D2)
                    f = f + params.xmu*sum(sum((D2'*D).^2));                    
                end
                df = f0 - f;
                if params.debug_armijo 
                    fprintf('df=%g\ttarget=%g\n',df,acceptable_decrease);
                end
                if  df >= acceptable_decrease
                    break;
                end
                % no luck: reduce step and acceptable decrease
                s = s * params.armijo_b;
                acceptable_decrease = acceptable_decrease * params.armijo_b;
            end
            if params.debug_armijo
                mdlsFigure('Armijo');
                plot(1:params.armijo_m,ads,'.-',1:params.armijo_m,dfs,'.-');
                pause(1);
            end
            if df < 0
                D = D0;
                df = 0;
            end
            dCost = dCost + df;
            if j >= params.armijo_m % no further significant improvement
                break;
                stuck_stats(j) = 1;
            end 
            armijo_m_stats(J) = j;
            if params.debug
                fprintf('\b\b\b\b\b\b\b\b\b');
            end
        end % descent iterations (J)
        if params.debug
            fprintf('\n');
            fprintf('Average armijo iterations:%g\n',mean(armijo_m_stats));
            fprintf('Stuck :%g%%\n',100*mean(stuck_stats));
            fprintf('dCost=%g\n',dCost);
        end
    end % PG update for the case mu1 ~= 0
    %
    % if requested, return 'stuck' flag
    %
    if nargout > 1
        varargout = cell(1,nargout-1);
        varargout{1} = sum(stuck_stats);
    end
end % function 
