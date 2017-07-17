% function [D,state] = mdlsLearnModel(data,params[,state])
% function def_params = mdlsLearnModel()
%
% This is a multi-purpose dictionary learning function capable of
% learning single and multi-class dictionaries while imposing
% coherence or cross-coherence on them.
% Without parameters, it returns a default set of parameters.
% With data and params, it starts anew from the provided learning data.
% If the state struct is provided, the learning state of the algorithm
% is initualized with it, so that a previous run can be continued.
%
% INPUT
% data   .......... Struct with learning data. For NC classes, X and Xt
%                   will be cells of size (1xNC) where each element is
%                   of the following format. If X and Xt  are matrices, 
%                   then a single dict is learned.
%   X ............. (MxN) training samples
%   Xt ............ (MxN') testing samples
%
% params ................ Struct with algorithm parameters:
%   params.D0 ............ initial dictionary(ies) (empty)
%                         If matrix, it is the same dict for all classes 
%                         If scalar, indicates number of atoms K in dict
%                         If cell, each element is the D0 for each class, 
%                         same as before, being either a matrix or a
%                         scalar.
%   params.max_iter ...... max number of iterations (100)
%   params.min_change .... minimum change in cost function 
%                         that justifies continuation of algorithm (1e-4)
%   params.batch_size .... If > 0, each iteration draws this number
%                         of samples from X for encoding (0)
%   params.test_size ..... If > 0, use a random subsample of this
%                         size for testing (0)
%   params.remember_factor For minibatch, how much weight is given to
%                         statistics learned from previous batches (0.8)
%   params.discard_unused_atoms if the frequency of usage of an atom is below
%                         this threshold (rel to 1/K), re-initialize it
%                         (1e-3)
%   params.discard_constant_patches  
%                         filter out patches whose variance is below
%                         this threshold (1e-3)
%   params.xval_step ..... Every how many interations x-validation is done
%   params.mu0 ....... ... coherence penalty (0):  
%                         mu*||DtD-I||_2^2, where mu is 
%                         mu0 * N / K^2
%   params.xmu0 .......... base value for the cross-coherence penalty:   
%                         xmu*||DtD||_2^2, where mu is computed as
%                         xmu = xmu0 * N / C / K^2, N = # of
%                         samples, K = # of atoms, C = # of classes.
%   params.mu_mode ....... How will we start imposing
%                         incoherence ([-1 50])
%                         If this is an integer rm,
%                         mu=mu0*heavyside(rm). This also affects xmu
%                         If this is two integers, it will produce a
%                         linear ramp from 0 to mu0 starting at the
%                         first number and ending at the second
%                         If it is two numbers and the first number
%                         is negative, it means that the starting time
%                         is automatically selected when the
%                         dictionaries are 'stuck'. Then the second
%                         number specifies the length of the ramp.
%
%   params.cache_dir .... Cache directory (cache/dict/)
%
%   params.model .............. Struct  specifying the regression model used:
%     params.model.reg_mode ... See scLasso: params.mode
%     params.model.reg_type ... 'l0','l1','moe','joe','l2','rwl2'. 'l1'
%                              is the default.
%     params.model.kappa ...... MOE shape parameter 
%     params.model.beta ....... MOE scale parameter
%
% state ..................... Learning state of previous run. Useful 
%                             for continuing. It also contains relevant
%                             useful statistics obtained through learning
%   state.AAt ............... empirical covariance matrix of coefficients
%   state.XAt ............... empirical cross-correlation of data and 
%                             coeffs.
%   state.rho ............... frequency of usage of each atom.
%   state.xrho .............. co-occurence matrix of atoms, that is, 
%                             abs(sgn(A))*abs(sgn(A))^T
%   state.N ................. effective number of samples processed.
%
% OUTPUT
%   D .......... The final dictionary or set of dictionaries.
%   state ...... Intermediate state. 
%
% AUTHOR: Ignacio Ramirez, ignacio.ramirez@gmail.com
%
function [D,state] = mdlsDictLearn(data,params,state)

  DEBUG_NONE = 0;
  DEBUG_MINIMAL = 1;
  DEBUG_MILD = 2;
  DEBUG_MEDIUM = 3;
  DEBUG_HEAVY = 4;

  %
  % PARSE INPUT
  %
  % no input means request default parameters
  if nargin == 0
      params = struct();
      params.max_iter     = 100;
      params.min_change   = 1e-4;
      params.batch_size   = 0; % no sampling of training set
      params.test_size    = 0; % no sampling of test set
      params.remember_factor = 0.8;
      params.xval_step = 10;
      params.mu0          = 0.0; % auto-incoherence
      params.xmu0         = 0.0; % cross incoherence
      params.mu_mode      = [-1 50]; % start when no further
                                     % improvement is seen for the
                                     % normal dictionary, and do a ramp
                                     % during 50 iterations.
      
      params.model = mdlsDefaultModelParams();
      params.cache_dir      = 'cache/dict';
      params.debug           = DEBUG_NONE;
      params.base_name       = 'global';
      params.discard_unused_atoms    = 0.001; % discard ones that are5
                                               % used 1% of the time
      params.discard_constant_patches = 0.001;
      params.D0 = [];
      D = params;
      state = [];
      return;
  end

  if ~isfield(data,'Xt')
      data.Xt = [];
  end

  if isempty(data.Xt)
      data.Xt = data.X;
      do_xval = false;
  else
      do_xval = true;
  end

  if ~iscell(data.X)
      data.X = {data.X};
  end

  if ~iscell(data.Xt)
      data.Xt = {data.Xt};
  end

  if params.test_size > 0
      data.Xt = subsample(data.Xt,params.test_size);
  end

  %
  % ------------------------------------
  % INITIALIZATION
  % ------------------------------------
  %
  if ~exist(params.cache_dir,'file')
      mkdir(params.cache_dir);
  end
  statname = sprintf('%s/%s-state.mat',params.cache_dir, ...
                     params.base_name);
  if ~exist('state','var')
      state = [];
  end

  if isempty(state)
      state = init_state(data,params);
  end
  
  if state.r > 0 % continue from previous run
      r0 =  state.r;
      params.max_iter = params.max_iter + r0;
  else
      r0 = 1;
  end

  %
  % ------------------------------------
  % MAIN LOOP
  % ------------------------------------
  %
  t0 = cputime();
  st0 = t0;
  xr = [];
  NC = state.NC;
  for r=r0:params.max_iter
      state.r = r;
      %
      % see how much time left
      %
      if (params.debug >= DEBUG_MILD) && (r == (r0+1))
          rem_time = (t1-t0)*(params.max_iter-r)/r;
          fprintf('** Batch %03d/%03d [%6.0f min to go]\n',...
                  r,params.max_iter,rem_time/60);
      end
      %
      % subsample
      %
      if params.batch_size > 0
          Xb = subsample(data.X,params.batch_size);
      else
          Xb = data.X;
      end
      %
      % sparse coding of each class
      %
      R = zeros(NC,1);
      T = zeros(NC,1);
      for c=1:NC
          rho = params.remember_factor;
          oldN = rho*state.N(c);
          newN = oldN + size(Xb{c},2);
          A = mdlsSparseCoding(Xb{c},state.D{c},params.model);
          state.XAt{c} = (oldN * state.XAt{c} + Xb{c}*A')*(1/newN);
          state.AAt{c} = (oldN * state.AAt{c} + A*A')*(1/newN);
          state.N(c) = newN;
      end
      %
      % dictionary update
      %
      mu_fac = comp_current_mu_fac(params,r);
      newD = state.D;
      for c=1:NC
          Kc = size(state.D{c},2);
          D2 = build_rival_dict(state,c);              
          K2 = size(D2,2);
          mu = params.mu0 / (Kc^2);
          xmu = params.xmu0 / (Kc*K2);
          dupar = mdlsDictUpdatePG();
          dupar.mu  = mu_fac * mu; 
          dupar.xmu = mu_fac * xmu;
          [Dc,stuck] = mdlsDictUpdatePG(state.D{c},...
                                        state.AAt{c},state.XAt{c},D2,...
                                        dupar);
          newD{c} = Dc;
      end

      %
      % compute current coherence
      %
      [maxco,meanco] = compute_coherence(state.D);
      state.maxco{end+1} = maxco;
      state.meanco{end+1} = meanco;

      %
      % cost function evaluation
      %
      is_xval_step = mod(r-1,params.xval_step)==0;
      if is_xval_step 
          xr = [xr r];
          for c=1:NC
              %
              % reconstruction energy
              %
              A = mdlsSparseCoding(data.Xt{c},state.D{c},params.model);
              R(c,1) = mean(mdlsModelEnergy(data.Xt{c},state.D{c},A,params.model));
          end
          mu = params.mu0 / (Kc^2);
          % plus mutual coherence
          T(c,1) = R(c,1) + mu*norm(state.D{c}'*state.D{c},'fro');
          %
          % plus cross coherence
          %
          for c2 = 1:(c-1)
              K2 = size(state.D{c2},2);
              xmu = params.xmu0 / (Kc*K2);
              T(c,1) = T(c,1) + xmu*norm(state.D{c}'*state.D{c2},'fro');
          end
          for c2 = (c+1):NC
              K2 = size(state.D{c2},2);
              xmu = params.xmu0 / (Kc*K2);
              T(c,1) = T(c,1) + xmu*norm(state.D{c}'*state.D{c2},'fro');
          end
          state.R = [state.R R];
          state.T = [state.T T];
          state.RT(end+1) = sum(state.R(end));
          state.TT(end+1) = sum(state.T(end));
          if r > 1
              % change in cost function
              dT = abs(state.TT(end)-state.TT(end-1))/(state.TT(end)+eps);
              % change in dict
              ndD =  0;
              nD  =  0;
              for c=1:NC
                  aux = (newD{c}-state.D{c}).^2;
                  ndD = ndD + sum(aux);
                  aux = (newD{c}).^2;
                  nD = nD + sum(aux);
              end
              dD = ndD / nD;
              if params.debug >= DEBUG_MINIMAL
                  fprintf('%4d/%4d:R=%8.2f\tT=%8.2g\tmax_co=%5.2f\tmean_co=%5.2f\tdcost=%6g\tdarg=%6g\n',...
                          state.r,...
                          params.max_iter,...
                          state.RT(end),...
                          state.TT(end),...
                          max(state.maxco{end}(:)),...
                          mean(state.meanco{end}(:)),...
                          dT,dD);
              end
              %
              % stopping condition
              %
              if (dT < params.min_change)
                  fprintf('Stopped (by dcost).\n');
                  break;
              elseif (dD < params.min_change)
                  fprintf('Stopped (by darg).\n');
                  break;
              end
          end
      end % if we are evaluating the cost function in this step
      state.D = newD;
      %
      % save state to cache every 10 min
      %
      if (cputime()-st0)>600 % 10 min
          state.rng_uni = rand('twister');
          stage.rng_normal = randn('state');
          state.r = state.r+1;
          save(statname,'state');
          state.r = state.r-1;
          st0 = cputime();
      end
      t1 = cputime();      
  end
  %
  % ------------------------------------
  % END
  % ------------------------------------
  %
  D = state.D;
end
      


function fac = comp_current_mu_fac(params,r)
    if params.mu_mode(1) >= 0 
        if length(params.mu_mode) == 1
            fac = (r> params.mu_mode);
        else
            fac = min(1,max(0,(r-params.mu_mode(1))/ ...
                            (params.mu_mode(2)- ...
                             params.mu_mode(1))));
        end
    else
        fac = 0;
    end
end

function idx = discard_atoms(params)
    if params.discard_unused_atoms > 0
        kt_usage = (state.Dusage(c,:) + 0.5)./( accumulated_N(state.c) + 1 );  % KT
                                                                               % estimate                
        if params.debug >= 2
            mdlsFigure(sprintf('Usage for %d',state.c),'dock',do_dock);
            plot(kt_usage,'x'); grid on;
        end
        thres = params.discard_unused_atoms;
        dead_atoms = find(kt_usage < thres); % used only 0.5%
                                             % of the time
        if length(dead_atoms) >= size(Dc,2)
            warning(['All atoms were discarded! Consider ' ...
                     'reducing threshold. I will reduce it ' ...
                     'anyway...']);
            params.discard_unused_atoms = params.discard_unused_atoms ...
                / 10;
        end
        aux = randperm(size(Xb,2));
        aux = aux( 1:min([ size(Dc,2) size(Xb,2) length(dead_atoms)]) );
        Dc(:,dead_atoms(1:length(aux))) = Xb(:,aux);
        state.Dusage(c,dead_atoms) = accumulated_N(state.c)/2; % reset stats
        if length(dead_atoms>0) && params.debug >= DEBUG_MEDIUM
            fprintf('Reset %d atoms:',length(dead_atoms));
            for i=dead_atoms
                fprintf('%d, ',i);
            end
            fprintf('\n');
        end
    end
end

function [maxco,meanco] = compute_coherence(D)
  NC=length(D);
  maxco=zeros(NC,NC);
  meanco=zeros(NC,NC);
  for i=1:NC
      for j=1:NC
          Gij = abs(D{i}'*D{j});
          if i ~= j
              maxco(i,j) = max(Gij(:));
              meanco(i,j) = mean(Gij(:));
          else
              Kc = size(D{i},2);
              Gij = Gij - diag(diag(Gij));
              maxco(i,j) = max(Gij(:));
              meanco(i,j) = sum(Gij(:))/(Kc*(Kc-1));
          end
      end
  end  
end

function state=init_state(data,params)
    state = struct();
    M = size(data.X{1},1);
    NC = length(data.X);
    K  = zeros(1,NC);    
    state.D = cell(1,NC);
    dipar = mdlsDictIni();
    if iscell(params.D0)
        for c=1:NC
            if length(params.D0) < c
                D0c = params.D0{1};
            else
                D0c = params.D0{c};
            end
            if isscalar(D0c)
                dipar.K = D0c;
                state.D{c} = mdlsDictIni(dipar);
            else
                state.D{c} = D0c;
            end
        end
    else % D0 is a scalar
        if isscalar(params.D0)
            dipar.K = params.D0;
            D0 = mdlsDictIni(dipar);
        else
            D0 = params.D0;
        end
        for c=1:NC
            state.D{c} = D0;
        end
    end
    

    KT = 0;
    for c=1:NC
        K(c) = size(state.D{c},2);
        KT = KT + K(c);
    end
    state.M = M;
    state.K = K;
    state.NC = NC;
    state.K = K;
    state.KT = KT;
    
    save(sprintf('%s/%s-params.mat',params.cache_dir,params.base_name),...
         'params');
    %
    % state of the algorithm
    %
    state.AAt = cell(1,NC);
    state.XAt = cell(1,NC);
    for c=1:NC
        state.AAt{c} = zeros(K(c),K(c));
        state.XAt{c} = zeros(M,K(c));
    end
    state.rho = zeros(NC,1);
    state.xrho = zeros(NC,1);
    %
    % initialize random number generator
    %
    rand('twister',M^2+M+K(1));
    randn('state',M^2+M+K(1));
    state.rng_uni = rand('state');    
    state.rng_normal = randn('state');
    
    state.r = 1;
    state.N = zeros(1,NC);
    
    state.R = [];
    state.T = [];
    state.RT = [];
    state.TT = [];
    state.maxco = {};
    state.meanco = [];
end

function Xb = subsample(X,batch_size)
    NC = length(X);
    Xb = cell(1,NC);
    for c=1:NC
        idx = ceil(rand(1,batch_size)*size(X{c},2));
        Xb{c} = X{c}(:,idx);
    end
end

function D2 = build_rival_dict(state,c)
    D2 = zeros(state.M,state.KT-state.K(c));
    K1 = 1;
    K2 = 0;
    for c2=1:c-1
        K2 = K2 + state.K(c2);
        D2(:,K1:K2) = state.D{c2};
        K1 = K2 + 1;
    end
    K2 = K2 + state.K(c);
    K1 = K2 + 1;
    for c2=c+1:state.NC
        K2 = K2 + state.K(c);
        D2(:,K1:K2) = state.D{c2};
        K1 = K2 + 1;
    end
end
