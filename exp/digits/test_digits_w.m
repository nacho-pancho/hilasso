

outdir0 = 'results/hilasso/digits/';
if ~exist(outdir0,'file')
    system(['mkdir -p ' outdir0]);
end

%
% noise level
%
if ~exist('sigma','var')
    sigma = 0; % no noise by default
end

%
% parameters to try with the different algorithms
%
lambda1L = [0.125 0.25 0.5 1.0 1.5];
lambda1 = [0.5 0.25 .125];
lambda2 = [4.0 3.0 2.0 1.0];

errorL2 = zeros(length(lambda1L),1);
errorH = zeros(length(lambda1),1);
errorG = zeros(length(lambda1),1);
errorC = zeros(length(lambda1),length(lambda2));

aSetL2 = zeros(length(lambda1L));
aSetG = zeros(length(lambda1));
aSetH = zeros(length(lambda1));
aSetC = zeros(length(lambda1),length(lambda2));

%
% initialize RNGs for repeatability
%
rand('twister',123987234);
randn('state',123987234);

%
% which algorithms shall we run
%
do_chilasso = true;
do_hilasso  = false;
do_lasso    = false;
do_glasso   = false;
do_cglasso  = false; 

%
% sort which sources will be active in each run
%
k = 1;
for i=1:10
    ii = randperm(10);
active{i} = ii(1:2)-1;
end
%
% load the data
%
load exp/hilasso/digits/dataUSPS.mat

N = 100;
K  = size(D{1},2);
NC = length(D);
aux = repmat(1:NC,K,1);
%
% Construct combined dictionary
%
Do = [];
groups = [];
for i=1:NC
    ng(i) = size(D{i},2);
    Do = [Do D{i}];
    groups = [groups i*ones(1,size(D{i},2))];
end

fh = fopen(sprintf('%s/digits-sigma%g.txt',outdir0,sigma),'w');
%try
%
% main batch
%
for a=1:length(active)
    %
    % Create data
    %
    fprintf('Creating data\n');
    if sigma == 0
        outdir = [outdir0 'digits'];
    else
        outdir = [outdir0 sprintf('digits-sigma=%g',sigma)];
    end
    [Y,X] = createDataDigits(data,N,active{a},k,sigma);
    Ao = zeros(NC*K,N);
    for aa = 1:length(active{a})
        Ao((K*active{a}(aa)+1):(K+1)*active{a}(aa),:) = 1;
        outdir = [outdir '-' num2str(active{a}(aa))];
    end
    gAo = group_act_set(Ao,K,1e-4);
    fprintf('outdir:%s\n',outdir);
    if ~exist(outdir,'file')
        mkdir(outdir);
    end
    %
    % LASSO
    %
    fprintf(fh,'GROUND TRUTH: %s\n', show_group_activity(mean(gAo')));
    if do_lasso
        for h = 1:length(lambda1L)
            % lasso method
            fres = sprintf('%s/digits-lambda1=%g-lasso.mat',...
                           outdir,lambda1L(h));
            fprintf(fh,'Lasso: lambda1=%g\t',lambda1L(h));
            if ~exist(fres,'file')
                [Xr,A] = lassoMethod(Y,D,lambda1L(h));
                save(fres,'A');
            else        
                load(fres);            
            end
            [Yols,Aols] = compute_ols(Y,Do,A); % experim
            clear Yols;
            Xo = cell(1,NC);
            for i=1:NC
                idx = find(groups == i);
                Xo{i} = D{i}*Aols(idx,:);
                v(i) = sum(sum((Aols(idx,:).^2)));
                %    Xr{i} = compute_ols(Z,Do,A)
            end

            se = separationError(X,Xo); % experim
            gA = group_act_set(A,K,1e-4);
            mean(gA');
            eA = group_energy(A,K);
            mean(eA');
            hamm = mdlsHammingDistance(gAo,gA);
            fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            %    errorL2(h) = separationError(X,Xr);
            aSetL(h) = mean(sum(A~=0));
            errorL(h) = se;
            % errorL2
        end
    end

% HILASSO AND HILASSO COLLABORATIVE
%
for h = 1:length(lambda1)
    
% 
%     % Group lasso method
%     [Xr,A] = groupLassoMethod(Ys,D,lambda1(h));
%     errorG(h) = separationError(Xs,Xr);
%     aSetG(h) = mean(sum(A~=0));
    prevAH = [];
    prevAC = [];
    for f=1:length(lambda2)
        %
        % Hilasso
        %
        if do_hilasso
            fres = sprintf('%s/digits-lambda1=%g-lambda2=%g-hilasso.mat',...
                           outdir,lambda1(h),lambda2(f));
            fprintf(fh,'HiLasso: lambda1=%g\tlambda2=%g\t',lambda1(h), ...
                    lambda2(f));
            if ~exist(fres,'file')
                [Xr,A,v] = HiLassoMethod2(Y,D,prevAH,lambda1(h),lambda2(f)/sqrt(K));
                save(fres,'A');
            else        
                clear A;
                load(fres);            
            end
            [Yo,Aols] = compute_ols(Y,Do,A); % experim
            clear Yo;
            Xo = cell(1,NC);
            for i=1:NC
                idx = find(groups == i);
                Xo{i} = D{i}*Aols(idx,:);
                v(i) = sum(sum((Aols(idx,:).^2)));
                %    Xr{i} = compute_ols(Z,Do,A)
            end
            se = separationError(X,Xo);
            gA = group_act_set(Aols,K);
            eA = group_energy(Aols,K);
            hamm = mdlsHammingDistance(gAo,gA);
            fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            errorH(h) = se;
            aSetH(h) = mean(sum(A~=0));
            prevAH = A; clear A;
        end % HiLasso
        %
        % C-HiLasso
        %
        clear A Aols Xr Xo Yo;
        if do_chilasso
            fres = sprintf('%s/digits-lambda1=%g-lambda2=%g-chilasso.mat',...
                           outdir,lambda1(h),lambda2(f));
            fprintf(fh,'C-HiLasso: lambda1=%g\tlambda2=%g\t',lambda1(h),lambda2(f));
            if ~exist(fres,'file')
                % collaborative lasso method
                lambdaL = 0; % no second pass with Lasso
                tol = 0.001;
                max_iter = 200;
                H = [];
                c = 10;
                [Xr,A,v] = HiLassoColMethodW(Y,D,prevAC,...
                                             lambda1(h),lambda2(f)/sqrt(N*K),...
                                             lambdaL,tol,[],max_iter,c);
                %se = separationError(X,Xr);
                save(fres,'A');
            else
                clear A;
                load(fres);
            end
            sum(sum(A))
            [Yo,Aols] = compute_ols(Y,Do,A); % experim
            clear Yo;
            Xo = cell(1,NC);
            for i=1:NC
                idx = find(groups == i);
                Xo{i} = D{i}*Aols(idx,:);
                v(i) = sum(sum((Aols(idx,:).^2)));
                %    Xr{i} = compute_ols(Z,Do,A)
            end
            se = separationError(X,Xo);
            gA = group_act_set(Aols,K);
            hamm = mdlsHammingDistance(gAo,gA);
            eA = group_energy(Aols,K);
            fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            errorC(h,f) = se;
            detectedGroups{h}{f} = v;
            aSetC(h,f) = mean(sum(A~=0));
            prevAC = A; 
            clear A Aols Xr Xo Yo;
        end % C-Hilasso
    end %lambda2 
end % lambda1
end % each sample combination
    %  fclose(fh);
    %catch 
    %  fclose(fh);
    %end
