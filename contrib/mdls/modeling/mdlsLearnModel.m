%
% Category: core
%
% Purpose: Dictionary learning using MOD/MOCOD
%
% Description:
%
% Input:
% params ...... algorithm parameters in a single struct.
%               These are the fields:
%
%        params.D0 ............ initial dictionary(ies). Either a matrix
%                               with one dictionary, or a cell where each element is a
%                               dictionary for a class.
%                               Defaults to [0 1]. This is used to scale
%                               the problem to [0 1] since penalty parameters
%                               depend on this range. This is NOT
%                               equivalent to normalization.
%        params.model ......... Variables specifying the sparse model:
%        params.model.reg_mode .... See scLasso: params.mode
%        params.model.reg_type ...... 'l1','moe','joe','l2','rwl2'. 'l1'
%                                     is the default.
%        params.model.kappa ......... MOE shape parameter 
%        params.model.beta .......... MOE scale parameter

%        params.mu0 ........... base value for the coherence penalty:   
%                               mu*||DtD-I||_2^2, where mu is computed as
%                               mu0 * N / K^2, N = # of samples, K = # of atoms
%        params.xmu0 .......... base value for the cross-coherence penalty:   
%                               xmu*||DtD||_2^2, where mu is computed as
%                               xmu = xmu0 * N / C / K^2, N = # of
%                               samples, K = # of atoms, C = # of classes.
%        params.mu_mode ...... How will we start imposing
%                              incoherence. 
%                              If this is an integer rm,
%                              mu=mu0*heavyside(rm). This also affects xmu
%                              If this is two integers, it will produce a
%                              linear ramp from 0 to mu0 starting at the
%                              first number and ending at the second
%                              If it is two numbers and the first number
%                              is negative, it means that the starting time
%                              is automatically selected when the
%                              dictionaries are 'stuck'. Then the second
%                              number specifies the length of the ramp.
%        params.positive ..... Enforce positivity in the dictionary and
%                              the sparse coding.
%        params.batch_size ... Size of minibatch. Set to 0 to disable
%                              minibatch learning.
%        params.test_size .... How many samples per test image do we take for control.
%        params.max_iter ..... max number of iterations.
%        params.min_change ... minimum change in cost function that justifies continuation of algorithm.
%
%        params.resume ....... Set to true to resume previous
%                              learning. Defaults to false.
%        params.xval_step . Every how many interations do we perform the control.
%        params.output_dir ... Output directory
%        params.testing_format .... Specifies the type of data for
%                              training. Can be:
%            'cache' ......... patches from cache. In this case, the
%                              'testing_data' member is a struct with
%                              information about the cache, such as the
%                              one returned by mdlsCreatePatchBuffer.
%            'mat' ........... mat file containing a single variable 'X' 
%            'explicit' ...... params.?_samples is X
%        params.training_format .... Same.
%        params.testing_data ..... Directory where the PGM testing images are stored. 
%        params.training_data .... Directory where the PGM training images are stored.
%        params.testing_labels ... For classification, give labels
%                                  corresponding to the data.
%        params.training_labels .. For classification, give labels
%                                  corresponding to the data.
% Output:
% D .......... The final dictionary or set of dictionaries.
%
function D = mdlsLearnModel(params)

DEBUG_NONE = 0;
DEBUG_MINIMAL = 1;
DEBUG_MILD = 2;
DEBUG_MEDIUM = 3;
DEBUG_HEAVY = 4;
DEBUG_INSANE = 5;
%
% PARSE INPUT
%

    if nargin == 0
        % no input means request default parameters
        params = struct();
        params.model = mdlsDefaultModelParams();
        params.mu0          = 0.0; % auto-incoherence
        params.xmu0         = 0.0; % cross incoherence
        params.mu_mode      = [-1 50]; % start when no further
                                       % improvement is seen for the
                                       % normal dictionary, and do a ramp
                                       % during 50 iterations.
        params.positive     = false;
        params.max_iter     = 1000;
        params.min_change   = 1e-3;
        params.batch_size   = 0;
        params.test_size    = 0;
        params.resume       = false;
        params.do_control_class   = false;
        params.xval_step = 10;
        params.remember_factor = 0.8;
        params.output_dir      = 'results/dictionaries';
        params.training_format = 'explicit';
        params.testing_format  = 'explicit';
        params.training_data   = [];
        params.testing_data    = [];
        params.training_labels = [];
        params.testing_labels  = [];
        params.update_method   = 'pg';
        params.debug           = 0;
        params.base_name       = 'global';
        params.discard_unused_atoms    = 0.001; % discard ones that are5
                                               % used 1% of the time
        params.discard_constant_patches = 0.001;
        params.dict_update     = mdlsDictUpdatePG();
        params.D0 = [];
        D = params;
        return;
    elseif ~isstruct(params)
        error('Input must be a struct. Call function with no arguments for a default set of parameters.');
    end
    if isequal(class(params.D0),'cell')
        [M,K] = size(params.D0{1});
        NC = length(params.D0);
        D = params.D0;
    else
        [M,K] = size(params.D0);
        NC = 1;
        D = {params.D0};
    end
    ref = zeros(NC);
    if isempty(D) || (M == 0) || (K == 0)
        error('Invalid initial dictionary.');
    end
    if params.model.L == 0
        params.model.L = M - 1;
    end
    if params.debug >= DEBUG_MEDIUM
        params.dict_update.debug = true;
    end
    %
    %
    % INITIALIZATION
    %
    do_dock = false;
    if ~exist(params.output_dir,'file')
        mkdir(params.output_dir);
    end
    %
    % learning status
    %
    ei = 0;
    est_batch_time = 0;
    r0 = 1;
    colmap = colormap(jet); 
    ncolors = size(colmap,1);
    global labels;
    if isempty(params.training_labels)
        labels = 1;
    else
        labels = sort(unique(params.training_labels));
    end
    labels = double(labels);
    
    % for plotting dictionaries on a grid
    ng = ceil(sqrt(NC));
    mg = ceil(NC/ng);

    global ftrain_cache;
    global multiclass;
    dDn = zeros(1,NC);

    multiclass = false;
    params.M = M;
    params.K = K;
    prefix = mdlsGetModelPrefix(params);
    %
    % log file
    %
    if params.debug >= DEBUG_MINIMAL
        fprintf('** Prefix: %s\n',prefix);
        fprintf('** Parameters: %s\n',prefix);
        params
        params.model
    end
    params2 =params;
    params.testing_data = [];
    params.training_data = [];
    params.testing_labels = [];
    params.training_labels = [];
    params.D0 = [];
    %mdlsWriteConfiguration(params);
    save([prefix '-params.mat'],'params');
    params = params2;
    clear params2;
    %
    % save configuration
    %

    %
    % state of the algorithm
    %
    finished = zeros(1,NC);
    AAt = cell(1,NC);
    XAt = cell(1,NC);
    Dusage = zeros(NC,K);
    Dmean = zeros(NC,K);
    Dvar = zeros(NC,K);
    %
    % evolution of cost function
    %
    cost = zeros(params.max_iter,NC);
    l1_energy = zeros(params.max_iter,NC);
    mean_coherence = zeros(params.max_iter,NC);
    max_coherence = zeros(params.max_iter,NC);
    max_cross_coherence = zeros(params.max_iter,NC);
    mean_cross_coherence = zeros(params.max_iter,NC);
    %
    % cross validation cost functions
    %
    control_t = [];
    xval_cost = zeros(params.max_iter,NC);
    xval_fit = zeros(params.max_iter,NC);
    xval_l1_energy = zeros(params.max_iter,NC);
    avg_nnz = xval_cost;

    error_rate = zeros(params.max_iter,1);
    accumulated_N = zeros(1,NC);

    tmpname = sprintf('%s-status.mat',prefix);
    dicname = sprintf('%s-dict.mat',prefix);

    if params.resume
        %
        % reuse already trained dictionary
        %
        if exist(dicname,'file')
            %
            % dictionary exists, load and return
            %
            d = dir(dicname);
	    if params.debug >= DEBUG_MINIMAL
                fprintf('* Model already computed on %s\n',d.date);
	    end
            load(dicname);
            return;
        elseif exist(tmpname,'file') 
            %
            % dictionary partially learned, continue
            %
            d = dir(tmpname);
            %params.resume = mdlsInput(...
            %    sprintf(['There is an existing session from %s.\nDo you want to ' ...
            %             'resume?'],d.date),0,'n');
        else
            params.resume = false;
        end
    end
    %
    % initialization
    %
    if params.debug >= DEBUG_MINIMAL
        fprintf('*** Initializing\n');
    end
    if NC > 1
        multiclass = 1;
    end
    %
    % initialize random number generator
    %
    rand('twister',M^2+M+K);
    %
    % load training data
    %
    if isequal(params.training_format,'explicit')
        X = params.training_data;        
        total_training_patches = size(X,2);
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Working with an explicit training set of size %d\n',...
                    total_training_patches);
        end
    elseif isequal(params.training_format,'mat')
        load(params.training_data);
        total_training_patches = size(X,2);
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Working with a loaded training set of size %d\n',...
                    total_training_patches);
        end
    elseif isequal(params.training_format,'cache')
        train_cache_info = params.training_data;
        total_training_patches = train_cache_info.patch_count;
        fname = [train_cache_info.buffer_prefix '-patches.bin'];
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Training cache file: %s\n',fname);
        end
        ftrain_cache = fopen(fname,'r');
        if ftrain_cache < 0
            error('Cannot open cache file.');
        end
        fprintf('* Working with training cache file of %d samples\n',...
                total_training_patches);
    else
        error('Unknown format');
    end
    if total_training_patches < params.batch_size
        params.batch_size = 0;
    end
    %
    % load all data at once if not in minibatch mode
    %
    if params.batch_size == 0
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Using full training dataset\n');
        end
        if isequal(params.training_format,'mat') || ...
                isequal(params.training_format,'explicit')
            Xb = X; clear X;
        else
            Xb = mdlsGetPatchesFromBuffer(ftrain_cache, train_cache_info);
        end
        labelsb = params.training_labels;
    end
    
    %
    % load testing data
    %
    if isequal(params.testing_format,'explicit')
        Xt = params.testing_data;        
        total_testing_patches = size(Xt,2);
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Working with an explicit testing set of size %d\n',...
                    total_testing_patches);
        end
    elseif isequal(params.testing_format,'mat')
        tmp = load(params.testing_data);
        Xt = tmp.X ; clear tmp;
        total_testing_patches = size(Xt,2);
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Working with an explicit testing set of size %d\n',...
                    total_testing_patches);
        end
    elseif isequal(params.testing_format,'cache')
        test_cache_info = params.testing_data;
        total_testing_patches = test_cache_info.patch_count;
        fname = [test_cache_info.buffer_prefix '-patches.bin'];
        fprintf('* Testing cache file: %s\n',fname);
        ftest_cache = fopen(fname,'r');
        if ftest_cache < 0
            error('Cannot open test cache file.');
        end
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Working with testing cache file of %d samples\n',...
                    total_training_patches);
        end
    end

    if params.test_size == 0
        %
        % use full testing data
        %
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Using full testing dataset\n');
        end
        if isequal(params.testing_format,'cache')
            Xt = mdlsGetPatchesFromBuffer(ftest_cache, test_cache_info);
        end
    elseif total_testing_patches > params.test_size 
        %
        % subsample testing data, if larger than requested test size
        %
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Sampling %d patches from test set\n', ...
                    params.test_size);
        end
        indexes = ceil( rand(1, params.test_size) * total_testing_patches);
        indexes=sort(indexes);
        if isequal(params.testing_format,'mat') || isequal(params.testing_format,'explicit')
            Xt = Xt(:,indexes);        
        else 
            Xt = mdlsGetPatchesFromBuffer(ftest_cache, test_cache_info, ...
                                          indexes);
            fclose(ftest_cache);
        end            
        if multiclass
            params.testing_labels = ...
                params.testing_labels(indexes);
        end
    end
    %
    % if no testing data is given, we do not perform xval
    %
    do_cross_validation = ~isempty(Xt);

    %
    % data structures
    %
    est_batch_time = 0;
    r0 = 1;
    for c=1:NC
        AAt{c} = zeros(K,K);
        XAt{c} = zeros(M,K);
    end

    if params.debug >= DEBUG_MINIMAL
        fprintf('* Total patches for training: %d\n',total_training_patches);
        fprintf('* Total patches for testing: %d\n', ...
                total_testing_patches);
    end

    if params.resume
        %
        % resume previous execution
        %
        fprintf('*** Resuming session from %s\n',d.date);
        load(tmpname);
        rand('twister',rng_state);            
        randn('state',rngg_state);            
        r0=r;
    end

    %
    % file handles have to be re-created every time, so saved one is
    % invalid.
    %
    if isequal(params.training_format,'cache')
        fname = [train_cache_info.buffer_prefix '-patches.bin'];
        if params.debug >= DEBUG_MINIMAL
            fprintf('* Training cache file: %s\n',fname);
        end
        ftrain_cache = fopen(fname,'r');
        if ftrain_cache < 0
            error('Cannot open training cache file.');
        end
    end
    %
    %
    %
    if params.debug >= DEBUG_MINIMAL
        fprintf('*** Main loop\n');
    end
    %
    % main loop
    %       
    for r=r0:params.max_iter
	%
	% if not adding cross-coherence, each dictionary is updated
	% separately and therefore, if one doesn't change, we skip it.
	% 
        if isempty(find(finished < 3))
            %
            % when adding incoherence in  automatic mode
            % the first time that the dictionaries get stuck 
            % is the starting time for adding incoherence
            %
            if (params.mu_mode(1) < 0) && ((params.mu0 > 0) || (params.xmu0 > 0))
                fprintf('** Start adding incoherence now.');
                params.mu_mode(1) = r;
                % automatic with ramp: second value is duration of
                % ramp
                if length(params.mu_mode) == 2
                    params.mu_mode(2) = params.mu_mode(2) + ...
                        params.mu_mode(1);
                end
                finished(:) = 0;
            else
                %
                % otherwise this is the END
                %
                if params.debug >= DEBUG_MILD
                    fprintf('\n*** All dictionaries have finished.');
                else
                    fprintf('*\n');
                end
                break;
            end
        end                    
        rem_time = NC*est_batch_time*(params.max_iter-r);
        if params.debug >= DEBUG_MILD
            fprintf('** Batch %03d/%03d [%6.0f minremaining.]\n',...
                    r,params.max_iter,rem_time/60);
        else
            %fprintf('.');
        end
        %
        % subsampling with reposition
        %
        if (params.batch_size > 0) && (total_training_patches > ...
                                       params.batch_size)
            if params.debug >= DEBUG_HEAVY
                fprintf(['* Minibatch subsampling training dataset\' ...
                         'n']);
            end
            indexes = ceil(rand(1,params.batch_size)*total_training_patches);
            indexes = sort(indexes);
            if multiclass
                labelsb = params.training_labels(indexes);
            end
            if isequal(params.training_format,'mat') || ...
                    isequal(params.training_format,'explicit')
                Xb = X(:,indexes);        
            else 
                %
                % get some random samples
                %
                Xb = mdlsGetPatchesFromBuffer(ftrain_cache, ...
                                              train_cache_info, ...
                                              indexes);
            end                
            if params.discard_constant_patches  > 0 
	        thres = params.discard_constant_patches;
                Vb = var(Xb);
                %fprintf('Discarding %d patches\n', );
                if sum(Vb <= thres) > 0
                    fprintf('x');
                    good = Vb > thres;
                    Xb = Xb(:,good);
                    if multiclass
                        labelsb = labelsb(good);
                    end
                end
            end
            if params.debug >= DEBUG_INSANE
                mdlsFigure('Batch','dock',do_dock,'nomargin',1); 
                imagesc(mdlsDictDisplay(Xb,'orient',1)); 
                colormap gray;
            end
            if params.debug >= DEBUG_HEAVY
                mdlsFigure('Variance histogram','dock',do_dock); 
                hist(var(Xb),25);
            end
        end % minibatch sampling
        maxMeanL = 0;
        leg=cell(1,NC);
        leg_mc=cell(1,2*NC);
        for c=1:NC          
            if params.debug >= DEBUG_MEDIUM && multiclass
                fprintf('** Class %d/%d\n',labels(c),NC);
            end
            ll = sprintf('%d',labels(c));
            leg_mc{2*c-1}=ll;
            leg_mc{2*c}=ll;
            leg{c}=ll; clear ll;
            t0 = cputime();
            Dc = D{c};
            if multiclass
                %
                % assemble 'enemy' dictionary 
                % for cross-mutual coherence
                D2 = zeros(M,K*(NC-1));
                for c2=1:c-1
                    D2(:,(1+K*(c2-1)):K*(c2)) = D{c2};
                end
                for c2=c+1:NC
                    D2(:,(1+K*(c2-2)):K*(c2-1)) = D{c2};
                end 
            else
                D2 = [];
            end
            if (finished(c) < 3) || (params.xmu0 > 0) 
	        % this class didn't finish yet
                               %
                               % take samples from this class
                               %
                if multiclass 
                    Xbc = Xb(:,labelsb == labels(c));
                else
                    Xbc = Xb;
                end
                
                %
                % adjust data-size dependent parameters
                %
                Nc = size(Xbc,2);
                old_Nc = params.remember_factor*accumulated_N(c);
                if r > 1
                    l1_energy(r,c) = l1_energy(r-1,c) * old_Nc;
                else
                    l1_energy(r,c) = 0;
                end
                accumulated_N(c) = old_Nc + Nc;
                %
                % since the full energy increases linearly with the number
                % of samples in the L2 and L1 terms, but not in the incoherence
                % terms, we need to adjust these penalties so that he incoherence
                % has the same weight no matter how many samples we are
                % encoding.
                %
                % on the other hand, the incohernece penalty terms are
                % proportional to K^2 (for self-cohernece) and to NC*K^2
                % for cross-coherence, so we also scale down these terms
                % so that they have the same weight no matter the size of
                % the dictionaries.
                mu  = params.mu0  *  accumulated_N(c) / (K^2); 
                xmu  = params.xmu0  * accumulated_N(c) / ((NC-1)*K^2);

                %
                % if mu_mode(1) < 0 this means that the first value is
                % chosen automatically (see below after 'stopping
                % criterion').
                %
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
                mu = fac*mu;
                xmu = fac*xmu;
                %
                % sparse coding
                %
                if params.debug >= DEBUG_HEAVY
                    fprintf(['* Coding %d patches with lambda=%g L=%g, ' ...
                             'penalty mode\n'],size(Xbc,2),params.model.lambda, ...
                            params.model.L);
                    fprintf('* Accumulated samples for this class: %d\n',...
                            round(accumulated_N(c)));
                end
                A = mdlsSparseCoding(Xbc,D{c},params.model);
                new_energy = mdlsModelEnergy(Xbc,D{c},A,params.model);
                l1_energy(r,c) = ( l1_energy(r,c) + sum(new_energy) ) / accumulated_N(c);
                if params.debug >= DEBUG_MEDIUM
                    fprintf('* Updating statistics\n');
                end
                XAt{c} = params.remember_factor * XAt{c} + Xbc*A';
                AAt{c} = params.remember_factor * AAt{c} + A*A';
                batch_usage = sum(A ~= 0, 2)';
                Dusage(c,:) = params.remember_factor * Dusage(c,:) + batch_usage;
                Dmean(c,:)  = params.remember_factor * Dmean(c,:) + (sum(A, 2)')./batch_usage;
                Dvar(c,:)   = params.remember_factor * Dvar(c,:) + (sum(A.^2, 2)')./batch_usage;
                params.dict_update.mu = mu;
                params.dict_update.xmu = xmu;
                if params.debug  >= DEBUG_INSANE
                    fprintf('* Updating dictionary\n');
                    mdlsFigure(sprintf('AAt class %d',c),'nomargin',1,'dock',do_dock); 
                    imagesc(AAt{c}); colormap gray;
                    mdlsFigure(sprintf('XAt class %d',c),'nomargin',1,'dock',do_dock); 
                    imagesc(XAt{c}); colormap gray;
                    % count number of times that an atom was used
                    %mdlsFigure(sprintf('Dusage for %d',c),'nomargin',1,'dock',do_dock);
                    %                imagesc(mdlsDisplayAtomUsage(Dusage));
                    mdlsFigure(sprintf('diag(AAt) class %d',c),'dock',do_dock); 
                    semilogy(diag(AAt{c}),'x'); grid on; colormap gray;
                end
                %
                % dictionary update
                %
                params.dict_update.positive = params.positive;
                [Dc,stuck] = mdlsDictUpdatePG(D{c},AAt{c},XAt{c},D2, ...
                                              params.dict_update);
                if max(sum(Dc.^2)) > (1+sqrt(eps)) % unnormalized
                    warning('Unnormalized dictionary?? ');
                    keyboard
                end
                %
                % resetting of unused atoms
                %
                % frequency is computed using a 
                % Kirchevskii-Trofimov estimator: 
                % theta(k) =  ( n_k + 1/2 ) / (n + 1)
                %
                % if an atom is reset, its usage is reset
                % so that it starts over with p = 0.5
                if params.discard_unused_atoms > 0
                    kt_usage = (Dusage(c,:) + 0.5)./( accumulated_N(c) + 1 );  % KT
                                                                              % estimate                
                    if params.debug >= 2
                        mdlsFigure(sprintf('Usage for %d',c),'dock',do_dock);
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
                    Dc(:,dead_atoms(1:length(aux))) = mdlsDictNormalize(Xb(:,aux));
                    Dusage(c,dead_atoms) = accumulated_N(c)/2; % reset stats
                    if length(dead_atoms>0)% params.debug
                        fprintf('Reset %d atoms:',length(dead_atoms));
                        for i=dead_atoms
                            fprintf('%d, ',i);
                        end
                        fprintf('\n');
                    end
                end
                clear A;
                dD{c} = Dc - D{c}; dDn(c) = norm(dD{c}(:)); % change in dictionary
                D{c} = Dc;

                mc = abs(Dc'*Dc); mc = mc-diag(diag(mc));
                mc = mc(:);
                max_coherence(r,c) = max(mc);
                mean_coherence(r,c) = mean(mc);
                
            else % if finished < 3 
                %
                % copy from previous iterations
                %
                l1_energy(r,c) = l1_energy(r-1,c);
                cost(r,c)      = cost(r-1,c);
                max_coherence(r,c) = max_coherence(r-1,c);
                mean_coherence(r,c) = mean_coherence(r-1,c);
            end
            %
            % cost function evaluation
            %
            if multiclass && ~isempty(D2)
                xc = abs(D2'*Dc);
                max_cross_coherence(r,c) = max(xc(:));
                mean_cross_coherence(r,c) = mean(xc(:));                
            end
            cost(r,c)      = l1_energy(r,c);
            if params.mu0 > 0
                mu = params.mu0 /K^2;
                cost(r,c) = cost(r,c) + mu*sum(sum((Dc'*Dc).^2)) ;
            end
            if multiclass && params.xmu0 > 0
                xmu = params.xmu0 / ((NC-1)*K^2);
                cost(r,c) = cost(r,c) + xmu*sum(sum(xc.^2));
            end

            colidx = 1;
            if NC > 1
                colidx = round(1+(c-1)/(NC-1)*(ncolors-1));
            end
            col = colmap(colidx,:); 
            clear colidx;
            %
            % Perform cross-validation
            %
            is_xval_step =  do_cross_validation && (mod(r-1,params.xval_step)==0);
            if is_xval_step
                if c==1
                    ei = ei+1;
                end
                if params.debug >= DEBUG_MILD
                    fprintf('* Cross-validating Reconstruction\n');
                end
                if multiclass
                    Xtc = Xt(:,params.testing_labels == labels(c));
                else
                    Xtc = Xt;
                end
                %
                % Control: sparse encode the same image and see if the energy is reduced
                %
                Ntc = size(Xtc,2);                
                Ar = mdlsSparseCoding(Xtc,D{c},params.model);
                avg_nnz(ei,c) = nnz(Ar)/size(Ar,2);
                xval_l1_energy(ei,c) = mean(mdlsModelEnergy(Xtc,D{c},Ar,params.model));
                xval_fit(ei,c) = 0.5*mean(sum((Xtc-D{c}*Ar).^2));
                xval_cost(ei,c) = xval_l1_energy(ei,c);
                control_t(ei) = r;

                mu  = params.mu0  / (K^2); 
                xmu = params.xmu0 / (NC*K^2);
                if mu > 0
                    xval_cost(ei,c) = xval_cost(ei,c) + mu*sum(sum((Dc'*Dc-eye(K)).^2));
                end
                if multiclass && xmu > 0
                    xval_cost(ei,c) = xval_cost(ei,c) + xmu*norm(D2'*Dc,'fro')^2 ;
                end
                xval_cost(ei,c) = xval_cost(ei,c);
            end % cross validation

            %
            % Display evolution (debug mode)
            %
            if params.debug >= DEBUG_MEDIUM
                if params.debug >= DEBUG_HEAVY
                 mdlsFigure('Dictionaries','dock',do_dock);
                 subplot(mg,ng,c);
                 imagesc(mdlsDictDisplay(Dc,'orient',0)); axis off; 
                 colormap gray; title(sprintf('D(%d)',r));
                 mdlsFigure('dDict','dock',do_dock);
                 subplot(mg,ng,c);
                 imagesc(mdlsDictDisplay(dD{c},'orient',0)); axis off; 
                 colormap gray; title(sprintf('D(%d)',r));
                end                
                mdlsFigure('Evolution','dock',do_dock);
                subplot(1,3+multiclass,1); 
                if multiclass 
                    if is_xval_step
                        plot(control_t,xval_l1_energy(1:ei,c),'.-','Color', ...
                             col);
                    end
                else
                    plot(1:r,l1_energy(1:r,c),'.-','Color', col);
                end
                grid on;  hold on;
                title('L1 Energy');
                
                subplot(1,3+multiclass,2); 
                if multiclass
                    if is_xval_step
                        plot(control_t,xval_cost(1:ei,c),'.-','Color', ...
                             col);
                    end
                else
                    plot(1:r,cost(1:r,c),'.-','Color',col);
                end
                grid on;  hold on;
                title('Total cost function');
                
                subplot(1,3+multiclass,3); 
                plot(1:r, max_coherence(1:r,c),'.-','Color',col); 
                grid on; hold on;
                plot(1:r, mean_coherence(1:r,c),'.-.','Color',col);
                title('Mutual Coherence');
                if multiclass
                    subplot(1,3+multiclass,4);
                    plot(1:r, max_cross_coherence(1:r,c),'.-','Color',col); 
                    grid on; hold on;
                    plot(1:r, mean_cross_coherence(1:r,c),'.-.','Color',col);
                    title('Cross Coherence');
                end
            end 

            %
            % stopping criterion
            %
            if do_cross_validation % based on testing data
                if is_xval_step && (ei > 1)
                    cost_dif = (xval_cost(ei-1,c) - xval_cost(ei,c))/abs(xval_cost(1,c));
                    if params.debug >= DEBUG_MINIMAL
                        fprintf('%4d/%4d: fit=%6g\tfit+reg=%6g\ttotal=%6g\t<nnz>=%3.0f\tmax_co=%5.2f\tmean_co=%5.2f\tdcost=%6g\tdD=%6g\n',...
                                r,params.max_iter,...
                                xval_fit(ei,c),...
                                xval_l1_energy(ei,c),...
                                xval_cost(ei,c),...
                                avg_nnz(ei,c),...
                                max_coherence(r,c), ...
                                mean_coherence(r,c),cost_dif,dDn(c));
                    end
                    %
                    % requires energy to be stuck for three iterations in a row
                    %
                    %if abs(cost_dif) < params.min_change
                    % with sign, avoids overfitting...
                    if cost_dif < params.min_change
                        %fprintf('** Stuck.',labels(c));
                        finished(c) = finished(c) + 1;
                    else 
                        finished(c) = 0; % back to activity
                    end
                end
            elseif r > 1 % based on training data
                cost_dif = (cost(r-1,c) - cost(r,c))/cost(1,c);
                if params.debug >= DEBUG_MINIMAL
                    fprintf('%4d/%4d:lasso cost=%8.2f\ttotal cost=%8.2g\tmax_co=%5.2f\tmean_co=%5.2f\tdcost=%6g\n',...
                            r,params.max_iter,l1_energy(r,c),cost(r,c),max_coherence(r,c), ...
                            mean_coherence(r,c),cost_dif);
                end
                %
                % requires energy to be stuck for three iterations in a row
                %
                % with sign, avoids overfitting...
                if cost_dif < params.min_change
                    finished(c) = finished(c) + 1;
                else 
                    finished(c) = 0; % back to activity
                end
            end % 

        end % for each class

        %
        % show evolution
        %
        if is_xval_step
            %
            % drawing
            %            
            if params.debug >= DEBUG_MEDIUM
                mdlsFigure('Evolution');
                subplot(1,3+multiclass,1); hold off; %legend(leg);
                subplot(1,3+multiclass,2); hold off; %legend(leg);
                subplot(1,3+multiclass,3); ; hold off; % legend(leg);
                if multiclass
                    subplot(1,3+multiclass,4); legend(leg_mc); hold off;
                end
                drawnow;
            end
        end
        %
        % save state
        %
        rng_state = rand('twister');
        rngg_state = randn('state');
        r = r+1;
        save(tmpname,'D','AAt','XAt','accumulated_N','Dusage','xval_cost',...
             'xval_l1_energy','mean_coherence','max_coherence',...
             'mean_cross_coherence','max_cross_coherence',...
             'error_rate','finished','r','ei','rng_state','rngg_state',...
             'est_batch_time');
        r = r-1;
        est_batch_time = 0.5*est_batch_time + 0.5*(cputime()-t0);
        if params.debug
            %profile viewer
        end
    end % batches
    %
    % save model
    %
    if params.debug >= DEBUG_MILD
        fprintf('** Saving trained dictionaries\n');
    end
    save(sprintf('%s-dict.mat',prefix) ,'D','Dusage','Dmean','Dvar');
    %
    % cleanup
    %
    if ftrain_cache >= 0
        fclose(ftrain_cache);
    end
end
