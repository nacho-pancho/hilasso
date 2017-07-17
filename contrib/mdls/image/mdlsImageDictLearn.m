%
% function D = mdlsLearnImageDict(params)
%
% Category: image processing
%
% Purpose: General purpose reconstructive dictionary learning for images images.
%
% Description:
%
% Input:
% training_dir ...... folder with training images
% params ............ algorithm parameters in a single struct.
%                     See help mdlsLearModel to see a description of them.
%                     This function sets some defaults of its own on these
%                     parameters, and adds a few more. These are:
%                       w (new) ...... width of the image patches. Defaults to 8.
%                       K ............ Size of dictionary. Defaults to 256.
%                       batch_size ... Set to 2560 by default.
%                       test_size .... Set to 5000 by default.                
%                       D0 ........... Empty by default. It is initialized  
%                                      as an average of 5 random patches plus
%                                      some slight noise.
%                                      corresponding to the data.
% Output:
%
% D .......... The final dictionary or set of dictionaries. If function
%              is called with no arguments, this is a set of default parameters.
%
% author Ignacio Ramirez <ignacio.ramirez@gmail.com>
%
function [D,state] = mdlsImageDictLearn(img_dir,params,state)
    %
    % SET DEFAULTS
    %
  DEBUG_NONE = 0;
  DEBUG_MINIMAL = 1;
  DEBUG_MILD = 2;
  DEBUG_MEDIUM = 3;
  DEBUG_HEAVY = 4;
    %
    % model parameters
    %
    %
    % PARSE INPUT
    %
    if nargin == 0
        % no input means request default parameters
        params.w            = 8;
        params.K            = params.w^2*4;
        params.max_iter     = 1000;
        params.min_change   = 1e-4;
        params.output_dir    = 'results/dictionaries';
        params.batch_size = 1000;
        params.test_size  = 10000;
        params.base_name  = 'images';
        params.xval_step  = 1;
        params.D0 = {};
        params.debug = 0;
        params.remember_factor = 0.9;
        params.xval_step = 10;
        params.mu0          = 0.5; % auto-incoherence
        params.mu_mode      = [200 300]; 
        params.discard_constant_patches = 1e-4;        
        params.model = mdlsDefaultModelParams();
        params.cache_dir      = 'cache/dict';
        params.debug           = DEBUG_NONE;
        D = params;
        return;
    end

    if ~isstruct(params)
        error('Input must be a struct. Call function with no arguments for a default set of parameters.');
    end

    if ~exist('state','var')
        state = [];
    end
    %
    % derived parameters
    %
    params.M  = params.w^2;    

    %
    % ------------------------------------
    % INITIALIZATION
    % ------------------------------------
    %
    %
    % prepare cache dir
    %
    if ~exist(params.cache_dir,'file')
        mkdir(params.cache_dir);
    end
    save(sprintf('%s/%s-params.mat',params.cache_dir,params.base_name),...
         'params');

    statname = sprintf('%s/%s-state.mat',...
                       params.cache_dir, ...
                       params.base_name);
    
    if isempty(state)
        state = init_state(img_dir,params);
    end
    
    if state.r > 0 % continue from previous run
        r0 =  state.r;
        params.max_iter = params.max_iter + r0;
    else
        r0 = 1;
    end
    
    %
    % Adapt dictionary
    %

    fprintf('Learning dictionary\n');
    %
    % ------------------------------------
    % MAIN LOOP
    % ------------------------------------
    %
    t0 = cputime();
    st0 = t0;
    xr = [];
    state.stuck_D = 0;
    state.stuck_T = 0;
    for r=r0:params.max_iter
        state.r = r;
        %
        % see how much time left
        %
        if (params.debug >= DEBUG_MILD) && (r == (r0+1))
            rem_time = (t1-t0)*(params.max_iter-r)/r;
            fprintf('** Estimated running time is %6.0f minutes **\n',...
                    rem_time/60);
        end
        %
        % sample training batch
        %
        Xb = sample_patches(img_dir,state.train_list, params, params.batch_size);
        %
        % sparse coding 
        %
        rho = params.remember_factor;
        oldN = rho*state.N;
        newN = oldN + size(Xb,2);
        A = mdlsSparseCoding(Xb,state.D,params.model);
        state.XAt = (oldN * state.XAt + Xb*A')*(1/newN);
        state.AAt = (oldN * state.AAt + A*A')*(1/newN);
        state.N = newN;
        %
        % coefficient statistics (for building models)
        %
        state.nnzA    = (oldN * state.nnzA + sum(A~=0,2)      ) * (1/newN);
        state.sumA    = (oldN * state.sumA + sum(A,2)         ) * (1/newN);
        state.sumAbsA = (oldN * state.sumAbsA + sum(abs(A),2) ) * (1/newN);
        state.sumA2   = (oldN * state.sumA2 + sum(A.^2,2)     ) * (1/newN);
        %
        % dictionary update
        %
        mu_fac = comp_current_mu_fac(params,r);
        mu = params.mu0 / (params.K^2);
        dupar = mdlsDictUpdatePG();
        dupar.mu  = mu_fac * mu; 
        [newD,stuck] = mdlsDictUpdatePG(state.D, state.AAt, ...
                                        state.XAt,[], dupar);
        %
        % compute current coherence
        %
        G = abs(state.D'*state.D);
        G = G - diag(diag(G));
        state.maxco(end+1) = max(G(:));
        state.meanco(end+1) = sum(G(:))/(params.K*(params.K-1));
        %
        % cost function evaluation
        %
        is_xval_step = (mod(r-1,params.xval_step)==0);
        if is_xval_step 
            xr = [xr r];
            %
            % reconstruction energy
            %
            A = mdlsSparseCoding(state.Xt,state.D,params.model);
            R = mean(mdlsModelEnergy(state.Xt,state.D,A,params.model));
            mu = params.mu0 / (params.K^2);
            %
            % plus mutual coherence
            %
            T = R + mu*norm(state.D'*state.D,'fro');
            state.R = [state.R R];
            state.T = [state.T T];
            state.nnz = [state.nnz nnz(A)/size(A,2)];
            if r > 1
                %
                % relative change in cost function
                %
                dT = (state.T(end-1)-state.T(end))/abs(state.T(end)+eps);
                %
                % relative change in argument
                %
                dD = norm(newD(:)-state.D(:)) / norm(newD(:)+eps);
                if params.debug >= DEBUG_MINIMAL
                    fprintf('%4d/%4d:NNZ=%4.1f\tR=%6g\tT=%6g\tmax_co=%4.2f\tmean_co=%4.2f\tdcost=%6g\tdarg=%6g\n',...
                            state.r,...
                            params.max_iter,...
                            state.nnz(end),...
                            state.R(end),...
                            state.T(end),...
                            state.maxco(end),...
                            state.meanco(end),...
                            dT,dD);
                end
                %
                % stopping condition
                %                
                if (dT < params.min_change)
                    state.stuck_T = state.stuck_T + 1;
                    if state.stuck_T >= 3
                        fprintf('Stopped (by dcost).\n');
                        break;
                    end
                else
                    state.stuck_T = 0;
                end
                if (dD < params.min_change)
                    state.stuck_D = state.stuck_D + 1;
                    if state.stuck_D >= 3
                        fprintf('Stopped (by darg).\n');
                    end
                    break;
                else
                    state.stuck_D = 0;
                end
            end
        end % if we are evaluating the cost function in this step
        %
        % update state to new dictionary
        %
        state.D = newD;
        if params.debug >= DEBUG_MEDIUM
            mdlsFigure('Dictionary');
            imagesc(mdlsDictDisplay(newD)); colormap gray;
        end
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
        %
        % end of main loop
        %
    end
    save(statname,'state');
    D = state.D;
end
    
function X=sample_patches(img_dir,img_list,params,num_patches)
    grabbed = 0;
    %
    % do not stop until we have num_patches
    %
    X = zeros(params.M,num_patches);
    w = params.w;
    
    while grabbed < num_patches
        i = ceil(rand()*length(img_list));
        I = imread([img_dir '/' img_list{i}]);
        I = double(I)*(1/255);
        [m,n] = size(I);
        % sample patches from random grid
        G = int32(w+rand(2,5*num_patches)*(m-w));  % take more than needed
        [Xb,DC,VAR]=mdlsGetPatches(I,w,G);
        %
        % keep only 'interesting' patches, that is, nonconstant
        %
        idx =  VAR > 1e-4;
        Xb = Xb(:,idx);
        Xb = Xb-repmat(DC(idx),params.M,1);
        Nb = min(num_patches-grabbed,size(Xb,2));
        X(:,(grabbed+1):(grabbed+Nb)) = Xb(:,1:Nb); 
        clear Xb;
        grabbed = grabbed + Nb;
    end
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


function state=init_state(img_dir,params)
    state = struct();
    %
    % initialize random number generator
    %
    rand('twister',458753); % monkey-typed
    randn('state',1293871);
    state.rng_uni = rand('state');    
    state.rng_normal = randn('state');
    %
    % read image dir
    %
    state.img_list = mdlsGetFileList(img_dir,'pgm');
    %
    % create training and cross-validation subsets
    %
    state.NI = length(state.img_list);
    if (state.NI == 0)
        error(sprintf('No images in training folder %s!',img_dir));
    end
    xval_idx = randperm(state.NI);
    state.train_idx = xval_idx(11:end);
    state.xval_idx = xval_idx(1:10);
    state.train_list = state.img_list(state.train_idx);
    state.xval_list  = state.img_list(state.xval_idx);
    %
    % cross validation samples
    %
    state.Xt = sample_patches(img_dir, state.xval_list, params, ...
                              params.test_size);
    %
    % initialize dictionary
    %
    if isempty(params.D0)
        fprintf('Initializing dictionary\n');
        % average 5 patches and add some noise
        state.D = zeros(params.M,params.K)
        for i =1:5
            state.D = state.D + sample_patches(img_dir, state.train_list, params, params.K);
        end        
        state.D = state.D * (1/5);
        state.D = mdlsDictNormalize(state.D + (5/255)*randn(size(state.D)));        
    else
        fprintf('Using a given initial dictionary\n');
        state.D = params.D0; 
        clear params.D0;
        if iscell(state.D)
            state.D = state.D{1};
        end
    end
    M = params.M; % convenience
    K = params.K;
    state.w = params.w;
    state.M = params.w^2;
    state.K = params.K;
    
    %
    % state of the algorithm
    %
    state.AAt = zeros(K,K);
    state.XAt = zeros(M,K);
    state.rho = 0;
    state.xrho = 0;
    
    state.r = 1;
    state.N = 0;
    
    state.R = [];
    state.T = [];
    state.nnz = [];
    state.maxco = [];
    state.meanco = [];
    %
    % statistics of coefficients
    %
    state.nnzA    = zeros(K,1);
    state.sumA    = zeros(K,1);
    state.sumAbsA = zeros(K,1);
    state.sumA2   = zeros(K,1);
end

