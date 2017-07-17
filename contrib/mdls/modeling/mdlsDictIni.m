%
% function D = mdlsDictIni(params)
%
% input:
%
% params .. struct with the following members. Call with no arguments to
%           get a default value.
%   method ...... Initialization method. Can be 'dct','random' or
%                 'patches'. Defaults to 'dct'.
%   K ........... Number of atoms. Defaults to dimension of sample vectors.
%   X ........... (for method='patches') Sample patches to initialize
%                 dictionary. 
%   labels ...... If present, multiple dictionaries are returned in D (as
%                 a cell array). One for each distinct label.
%   layers ...... (for method='patches') Number of patches to be averaged to
%                 form an initial atom. Default = 1
%   sigma ....... (for method='patches') Noise std. dev. added to
%                 atoms. Defaults to 0.
%   seed ........ Random seed. Defaults to a fixed function of K and dim
%   scramble_dc . (for method='dct') Replace DC patch by a random
%                 Gaussian sample.
%   dim ......... Dimension of atoms.
%   file ........ Load dictionary from a file.
%
% output:
%
% D ....... initialized dictionary or dictionaries (if labels present).
%
function D = mdlsDictIni(params)
    if nargin==0
        params = struct();
        params.method = 'dct';
        params.K = 0;
        params.X = [];
        params.layers = 5;
        params.sigma = 0;
        params.seed = 0;
        params.scramble_dc = false;
        params.dim = 0;
        params.file = '';
        params.labels = [];
        params.variance_threshold = 0;
        D = params;
        return;
    end
    
    multiclass = ~isempty(params.labels);
    if multiclass
        classes = sort(unique(params.labels));
        NC = length(classes);
    else
        NC = 1;
    end
    if params.K == 0
        params.K = params.dim;
    end
    if multiclass && (length(params.K) == 1)
        params.K = params.K*ones(1,NC);
    end
    MD = cell(1,NC);
    for c=1:NC
        %
        % preprocessing
        %
        Kc = params.K(c);         

        if ~isempty(params.X)
            if multiclass
                Xc = params.X(:,params.labels==classes(c));
            else
                Xc = params.X;
            end
            size(Xc)
            M  = size(Xc,1);
            Nc = size(Xc,2);
            if params.variance_threshold > 0
                Xv = var(Xc);
                [vv,vi] = sort(Xv);
                good_ones = find(vv > params.variance_threshold);
                if length(good_ones) < Kc
                    vi = vi(1:Kc);
                else
                    vi = vi(good_ones);
                end
                Xc = Xc(:,vi);
                fprintf('%d samples survived variance threshold of %g\n',...
                        length(vi),params.variance_threshold);
            end
        end
        %
        % actual initialization
        %
        switch (params.method) 

          case 'dct'
            D = mdlsGenDCTDictionary(params.dim, Kc, params.scramble_dc);
            
          case 'patches'
            if isempty(params.X)
                error('Must specify training patches');
            end
            if Nc <= Kc
                D = Xc;
            else
                D = mdlsGenPatchesDictionary(Xc, Kc, params.sigma, ...
                                             params.layers, params.seed);
            end

          case 'variance'
            % non-uniform random sampling, with probability given by the
            % variance of the patches
            Xv = sqrt(var(Xc));
            pv = Xv/sum(Xv); % create a PDF out of this
            % cumulative distribution
            aux = repmat(1:Nc,Nc,1); 
            aux = (aux-aux')>=0; % upper triangular
            cv  = pv*aux;         
            % without reposition
            D = zeros(M,Kc);
            ki = 0;
            unsel = 1:Nc;
            while ki < Kc
                x = rand();
                k = find(cv<=x, 1, 'last'); 
                if unsel(k) 
                    D(:,ki+1) = Xc(:,k);
                    unsel(k) = 0;
                    ki = ki + 1;
                end
            end
            
          case 'greedy'
            %
            % greedy incoherent dictionary initialization
            % take first one at random, next one maximally incoherent
            % with previous ones. 
            %
            unsel = 1:Nc;
            aux = randperm(Nc);
            D = zeros(M,Kc);
            D(:,1)= Xc(:,aux(1));
            unsel(aux(1))=0;
            ki = 1;
            while ki < Kc
                %D = [D zeros(M,1)];
                minmc = 1e10;
                for kk = find(unsel)
                    %D(:,end) = Xc(:,kk);
                    mc = min(abs(Xc(:,kk)'*D(:,1:ki)));
                    if minmc > mc
                        minmc = mc;
                        k = kk;
                    end
                end
                D(:,ki+1) = Xc(:,k);
                unsel(k) = 0;
                ki = ki + 1;
            end

          % case 'trace'
          %   %
          %   % experimental mode: optimize initial dictionary
          %   % in terms of its 'perimeter'
          %   %
          %   %
          %   % solve SDP
          %   %
          %   %cvx_precision best;
          %   fprintf('Solving max trace(DDt)\n');
          %   cvx_begin
          %     variables z(Nc,1);
          %     minimize( -trace(Xc*sparse(diag(z))*Xc') );
          %     subject to
          %     sum(z)==K;
          %     0<=z;
          %     z<=1;
          %   cvx_end
          %   [zsorted,zi] = sort(z,'descend');
          %   D = Xc(:,zi(1:K));
          %   if params.sigma > 0
          %       D = D + params.sigma*randn(size(D));
          %   end

          % case 'logdet'
          %   if Nc <= Kc
          %       D = params.X;
          %   else
          %       %
          %       % this problem requires *a lot* of memory, so
          %       % we have either to split it and solve partial problems
          %       % or to do an initial refinement of the candidates, or both.
          %       %
          %       Nc = size(Xc,2)
          %       if Nc > 400
          %           %
          %           % split in two
          %           %
          %           fprintf('logdet: split %d samples in two\n',Nc);
          %           aux = randperm(Nc);
          %           params2 = params;
          %           params2.labels = [];
          %           params2.X = Xc(:,aux(1:floor(Nc/2)));
          %           fprintf('logdet: first subproblem');
          %           D1 = mdlsDictIni(params2);
          %           fprintf('logdet: second subproblem');
          %           params2.X = Xc(:,aux(floor(Nc/2)+1:end));
          %           D2 = mdlsDictIni(params2);
          %           Xc = [D1 D2];
          %       end
          %       Nc = size(Xc,2);
          %       %
          %       % solve SDP
          %       %
          %       %cvx_precision best;
          %       fprintf('Solving max log det(DDt) with %d patches\n',Nc);
          %       Im = 0.1*speye(M); % to make sure we have strict positive
          %                           % definite matrices...
          %       cvx_begin
          %       variables z(Nc,1);
          %       minimize( -log_det(Xc*diag(z)*Xc' + Im) );
          %       subject to
          %       sum(z)==Kc;
          %       0<=z;
          %       z<=1;
          %       cvx_end
          %       [zsorted,zi] = sort(z,'descend');
          %       D = Xc(:,zi(1:Kc));
          %       if params.sigma > 0
          %           D = D + params.sigma*randn(size(D));
          %       end
          %   end

          case 'random'
            D = mdlsGenRandomDictionary( params.dim,Kc,params.seed );

          case 'file'
            D = load(params.method);
            if isempty(D) || ...
                    (~isequal(class(D),'double') && ...
                     ~isequal(class(D),'single')) || ...
                    (size(D,1) ~= M)
                error(sprintf('File %s did not contain a valid dictionary\n',params.method));
            end                
          otherwise
            error('Unknown initialization method or file not found.');
        end
        %
        % multiple dictionaries are stored in a cell
        %
        if multiclass
            MD{c} = mdlsDictNormalize(D);
        end
    end
    %
    % in this case D is a cell array of NC dictionaries
    %
    if multiclass
        D = MD;
    end
end
