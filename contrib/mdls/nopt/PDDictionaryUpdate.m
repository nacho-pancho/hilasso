%
% function D = PDDictionaryUpdate(data,D,A,params)
%
% Name: PDDictionaryUpdate
%
% Category: core function
%
% Description: dictionary update for the PD problem.
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
function D = PDDictionaryUpdate(data,D,A,params)

  method = params.updateMethod;
  stepsizeRule = params.updateStepsizeRule;
  incremental = params.updateIncremental;

  %
  % Use the sparse version of the matrix of coefficients
  %
  if ~issparse(A)
    A = sparse(A);
  end
  X = data.X;
  mu = params.mu;
  [n,K] = size(D);
  Idn = eye(n); % nxn identity matrix
  IdK = eye(K); % KxK identity matrix
  %
  % Armijo's parameters
  %
  beta = 0.25;
  sigma = 1e-1;
  s = 1;

  C = 1;
  eta = params.eta;

  %
  % MOCOD-like method
  %
  if isequal(params.updateMethod,'MOCOD')
    % precompute some constants
    XAt = X*A';
    AAt = A*A';
    D0 = D;
    iter = 1;
    dif = realmax;
    rand('twister',101); % used by condest
    c=condest(AAt);
    if c > 1e3
      [isPD,iota]=IsPositiveSemidefinite(AAt);
      if ~isPD
        fprintf('warning: AAt badly conditioned (%f). Iota = %f\n',c,iota);
        AAt = AAt + iota*eye(size(AAt));
      end
    end
    if (mu > 0)  || (eta > 0) % do fixed point iteration
      while iter <= params.updateIters && (dif > params.updateTolerance)
        if params.eta > 0
          D = (XAt+2*(mu+eta)*C*D0)*inv(AAt + 2*mu*(D0'*D0) + 2*eta*diag(sum(D0.*D0)));
          % other variant
          % D = (XAt)*inv(AAt + 2*mu*(D0'*D0) + 2*eta*diag(sum(D0.*D0)) - 2*(mu+eta)*IdK );
        else
          D = (XAt+2*mu*C*D0)*inv(AAt + 2*mu*(D0'*D0));
        end
        dif = norm(D-D0,'fro')/(K*K);
        if params.debug
          fprintf('||D-D0||=%f\n',dif);
        end
        D0=D;
        iter = iter + 1;
      end
    else
      D = XAt*inv(AAt); % first update : MOCOD
    end
    % end MOCOD
  else if isequal(params.updateMethod,'MOCOD2')
      % precompute some constants
      num = X*A';
      den = A*A';
      if (mu > 0)
        den = den +  2*mu*(D'*D -IdK);
      end
      if (eta > 0)
        den = den + 2*eta*(diag(sum(D.*D)) - IdK);
      end
      rand('twister',101); % used by condest
      c=condest(den);
      if c > 1e3
        fprintf('warning: den badly conditioned (%f). Using pseudoinverse.\n',c);
        invden = pinv(den);
      else
        invden = inv(den);
      end
      D = num*invden;
    else
      %
      % Descent methods
      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Do the loop for each atom: Compute the new dictionary
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      tic;
      dotCount = 0;
      Dnew = D;
      armijoIters=zeros(1,K);
      for p = 1:K % randperm(K)
        if ~mod(p,25)
          fprintf('.');
          dotCount = dotCount + 1;
        end
        %
        % Select one subset of the data for speed up
        %
        Ap_ = A(p,:); % p-th row of the coefficients
        D_p = D(:,p); % p-th column of the dictionary
        %
        % hessian is computed only once per atom
        %
        if isequal(method,'ND')
          Hk = 2*(Ap_*Ap_')*Idn + 4*(D*D' + D_p*D_p' + (D_p'*D_p - C)*Idn);
          if params.eta > 0
            Hk = Hk + params.eta*(4*(2*D_p*D_p' + (D_p'*D_p-1)*Idn)); % additional term
          end
          [isPSD,iota] = IsPositiveSemidefinite(Hk);
          if ~isPSD
            LogFile(params,'write',sprintf('The hessian matrix is not Positive Semidefinite (atom:%d)\n',p));
            Sk = inv(Hk + iota*Idn);
          else
            Sk = inv(Hk);
          end
        end
        AAp = sparse(A*Ap_');
        XAp = sparse(X*Ap_');
        %
        % Prepare the first iteration and run ....
        %
        gk = -2*XAp + D*AAp + 4*mu*(D*D' - Idn)*D_p;
        if params.eta > 0
          gk = gk + params.eta*(4*(D_p'*D_p-1)*D_p);
        end
        normG = norm(gk);
        iterk = 1;
        Rp = X - (D*A - D_p*Ap_); % residual reconstruction without Dp
        while (iterk <= params.updateIters) && (normG > params.updateTolerance)
          %
          % Compute the descent direction, dk.
          %
          if isequal(method,'SD')
            %
            % Steepest Descent
            %
            dk = -gk;
          elseif isequal(method,'ND')
            %
            % Newton Direction
            %
            dk = -Sk*gk;
          elseif isequal(method,'CG')
            %
            % Conjugate Gradient
            %
            % Approximate the quadratic matrix Q by the hessian
            if isequal(iterk,1)
              dk = -gk;
            else
              % betak = (gk'*gk)/(gkp'*gkp);
              gkpn = gkp'*gkp;
              if gkpn < 1e-5
                fprintf(2,'CG: gkp is null. skipping\n');
                dk = -gk; % reset directions
              else
                betak = gk'*(gk - gkp)/(gkp'*gkp); % See Bertseka, eq. (1.172)
                dk = -gk + betak*dkp;
              end
            end
          end
          %
          % Update the atom
          %
          if isequal(stepsizeRule,'FIX')
            %
            % Fix stepsize
            %
            Dp1 = D_p + params.updateStep * dk;
          elseif isequal(stepsizeRule,'ARO')
            %
            % Armijo rule
            %
            fit = Rp - D_p*Ap_;
            %reg = D'*D_p - IdK(:,p);
            reg1 = D'*D - IdK;
            reg2 = sum(D.*D) - 1;
            iniCost = sum(sum(fit.*fit)) + mu*sum(sum(reg1.*reg1)) + eta*reg2*reg2';
            m = 0; goahead = true;
            tD = D; % tA = A;
            predDecrease = sigma*s*gk'*dk;
            %fprintf('Armijo: pred: %f\n',predDecrease);
            while goahead && (m < 10)
              Dp1 = D_p + beta^m*s*dk;
              tD(:,p) = Dp1;
              fit = Rp - Dp1*Ap_ ;
              reg1 = tD'*tD - IdK;
              reg2 = sum(tD.*tD) - 1;
              newCost = sum(sum(fit.*fit)) + mu*sum(sum(reg1.*reg1)) + eta*reg2*reg2';
              %fprintf('Armijo: iter=%2d pred=%f cost: old=%f new=%f\n',m,predDecrease,iniCost,newCost);
              deltaCost = newCost - iniCost;
              if (deltaCost <= predDecrease)
                goahead = false;
              else
                m = m + 1;
              end
              predDecrease = predDecrease*beta;
            end % while (Armijo)
            armijoIters(1,p)=m;
            %fprintf(1,'Total iterations: %d\n',m)
          elseif isequal(stepsizeRule,'EXT')
            %
            % Exact line search, only with CG
            %
            step = -(dk'*gk)/(dk'*Hk*dk);
            Dp1 = D_p + step * dk;
          end
          %
          % Prepare the next iteration
          %
          D_p = Dp1;
          % Insert the updated atom in the dictionary
          if incremental
            D(:,p) = D_p;
          else
            Dnew(:,p) = D_p;
          end
          iterk = iterk + 1;
          if isequal(method,'CG')
            dkp = dk; gkp = gk;
          end
          gk = -2*XAp + D*AAp + 4*mu*(D*D' - C*Idn)*D_p;
          if params.eta > 0
            gk = gk + params.eta*(4*(D_p'*D_p-1)*D_p);
          end
          normG = norm(gk);
        end % iterations for this atom
        drawnow
      end % for each atom
      if ~incremental
        D=Dnew;
      end
      for i=1:dotCount
        fprintf('\b');
      end
      fprintf('Median Armijo iterations: %d\n',median(armijoIters));
    end % end not MOCOD

  end % function
