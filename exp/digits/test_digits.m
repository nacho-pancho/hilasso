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
% Lasso:
lambda1L = [1.0];
% GL and CGL
lambda2G = [10];
% HL and CHL
lambda1 = [0.5];
lambda2 = [1.5]; 

% initialize RNGs for repeatability
%
rand('twister',123987234);
randn('state',123987234);

%
% which algorithms shall we run
%
do_chilasso = false;
do_hilasso  = true;
do_lasso    = false;
do_glasso   = false;
do_cglasso  = false; 

%
% sort which sources will be active in each run
%
active = cell(1,8+10);
k = 1;
for i=1:8
   ii = randperm(10);
   active{i} = ii(1:2)-1;
end
for i=9:18
   active{i}=i-9;
end
%
% load the data
%
load data/dataUSPS.mat

N = 10;
K  = size(D{1},2);
NC = length(D);
aux = repmat(1:NC,K,1);
reg = 1e-2;
ga_thres = 1e-5; % threshold to declare a group is active
%olsfun = @mdlsOLS;
olsfun = @compute_ols2;
%olsfun = @compute_ols;
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

NA = length(active);
errorL = zeros(NA,length(lambda1L));
errorGL = zeros(NA,length(lambda2G));
errorCGL = zeros(NA,length(lambda2G));
errorHL = zeros(NA,length(lambda1),length(lambda2));
errorCHL = zeros(NA,length(lambda1),length(lambda2));

hammL = zeros(NA,length(lambda1L));
hammGL = zeros(NA,length(lambda2G));
hammCGL = zeros(NA,length(lambda2G));
hammHL = zeros(NA,length(lambda1),length(lambda2));
hammCHL = zeros(NA,length(lambda1),length(lambda2));

%try
%
% main batch
%
for a=1:NA
    %
    % Create data
    %
    %fprintf('Creating data\n');
    if sigma == 0
        outdir = [outdir0 'digits'];
    else
        outdir = [outdir0 sprintf('digits-sigma=%g',sigma)];
    end
    [Y,X] = createDataDigits(data,N,active{a},k,sigma);
    Ao = zeros(NC*K,N);
    for aa = 1:length(active{a})
        Ao((K*active{a}(aa)+1):K*(active{a}(aa)+1),:) = 1;
        outdir = [outdir '-' num2str(active{a}(aa))];
    end
    gAo = group_act_set(Ao,K,ga_thres);
    fprintf('outdir:%s\n',outdir);
    if ~exist(outdir,'file')
        mkdir(outdir);
    end
    %
    % LASSO
    %
    fprintf(fh,'GROUND TRUTH: %s\n', show_group_activity(mean(gAo')));
    fprintf('GROUND TRUTH: %s\n', show_group_activity(mean(gAo')));
    if do_lasso
        for h = 1:length(lambda1L)
            % lasso method
            fres = sprintf('%s/digits-lambda1=%g-lasso.mat',...
                           outdir,lambda1L(h));
            fprintf(fh,'Lasso: lambda1=%g\t',lambda1L(h));
            fprintf('Lasso: lambda1=%g\t',lambda1L(h));
            if ~exist(fres,'file')
                [Xr,A] = lassoMethod(Y,D,lambda1L(h));
                save(fres,'A');
            else        
                load(fres);            
            end
            [Yols,Aols] = olsfun(Y,Do,A,reg); % experim
            clear Yols;
            Xo = cell(1,NC);
            for i=1:NC
                idx = find(groups == i);
                Xo{i} = D{i}*Aols(idx,:);
                v(i) = sum(sum((Aols(idx,:).^2)));
                %    Xr{i} = olsfun(Z,Do,A)
            end

            se = separationError(X,Xo); % experim
            gA = group_act_set(A,K,ga_thres);
            mean(gA');
            eA = group_energy(A,K);
            mean(eA');
            hamm = mdlsHammingDistance(gAo,gA);
            fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            fprintf('se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            %    errorL2(h) = separationError(X,Xr);
            errorL(a,h) = se;
            hammL(a,h) = hamm;
            % errorL2
        end
    end

%
% GLASSO AND C-GLASSO
%
for h = 1:length(lambda2G)
    %
    % GLASSO
    %
    prevAG = [];
    if do_glasso
        fres = sprintf('%s/digits-lambda2=%g-glasso.mat',...
                       outdir,lambda2G(h)-3);
        fprintf(fh,'GLasso: lambda2=%g\t',lambda2G(h)-3);        
        fprintf('GLasso: lambda2=%g\t',lambda2G(h)-3);        
        if ~exist(fres,'file')            
            [Xr,A] = groupLassoMethod(Y,D,prevAG,(lambda2G(h)-3)/sqrt(K));
            save(fres,'A');
        else        
            load(fres);            
        end
        prevAG = A;
        %[Yols,Aols] = olsfun(Y,Do,A,reg); % experim        
        %clear Yols;
        Xr = cell(1,NC);
        for i=1:NC
            idx = find(groups == i);
            Xr{i} = D{i}*A(idx,:);
        end        
        se = separationError(X,Xr); % experim
        gA = group_act_set(A,K,ga_thres);
        mean(gA');
        eA = group_energy(A,K);
        mean(eA');
        hamm = mdlsHammingDistance(gAo,gA);
        fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                show_group_activity(mean(eA')));
        fprintf('se=%g\thamm=%g\tact=%s\n',se,hamm,...
                show_group_activity(mean(eA')));
        %    errorL2(h) = separationError(X,Xr);
        errorGL(a,h) = se;
        hammGL(a,h) = hamm;
        % errorL2
    end
    prevAC = [];
    %
    % C-GLASSO
    %
    if do_glasso
        fres = sprintf('%s/digits-lambda2=%g-cglasso.mat',...
                       outdir,lambda2G(h));
        fprintf(fh,'C-GLasso: lambda2=%g\t',lambda2G(h));
        fprintf('C-GLasso: lambda2=%g\t',lambda2G(h));
        if ~exist(fres,'file')
            [Xr,A] = GLassoColMethod(Y,D,prevAC,lambda2G(h)/sqrt(K)/sqrt(size(Y,2)));
            save(fres,'A');
        else        
            load(fres);            
        end
        prevAC = A;
        Xr = cell(1,NC);
        for i=1:NC
            idx = find(groups == i);
            Xr{i} = D{i}*A(idx,:);
        end        
        se = separationError(X,Xr); % experim
        gA = group_act_set(A,K,ga_thres);
        eA = group_energy(A,K);
        hamm = mdlsHammingDistance(gAo,gA);
        fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                show_group_activity(mean(eA')));
        fprintf('se=%g\thamm=%g\tact=%s\n',se,hamm,...
                show_group_activity(mean(eA')));
        %    errorL2(h) = separationError(X,Xr);
        errorCGL(a,h) = se;
        hammCGL(a,h) = hamm;
        % errorL2
    end
end

%
% HILASSO AND HILASSO COLLABORATIVE
%
for h = 1:length(lambda1)
    
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
            fprintf('HiLasso: lambda1=%g\tlambda2=%g\t',lambda1(h), ...
                    lambda2(f));
            if ~exist(fres,'file')
                [Xr,A,v] = HiLassoMethod2(Y,D,prevAH,lambda1(h),lambda2(f)/sqrt(K));
                save(fres,'A');
            else        
                clear A;
                load(fres);            
            end
            [Yo,Aols] = olsfun(Y,Do,A,reg); % experim
            clear Yo;
            Xo = cell(1,NC);
            for i=1:NC
                idx = find(groups == i);
                Xo{i} = D{i}*Aols(idx,:);
                v(i) = sum(sum((Aols(idx,:).^2)));
                %    Xr{i} = olsfun(Z,Do,A)
            end
            se = separationError(X,Xo);
            gA = group_act_set(Aols,K);
            eA = group_energy(Aols,K);
            hamm = mdlsHammingDistance(gAo,gA);
            fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            fprintf('se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            errorHL(a,h,f) = se;
            hammHL(a,h,f) = hamm;
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
            fprintf('C-HiLasso: lambda1=%g\tlambda2=%g\t',lambda1(h),lambda2(f));
            if ~exist(fres,'file')
                % collaborative lasso method
                lambdaL = 0; % no second pass with Lasso
                tol = 0.0001;
                H =[];
                max_iter = 100;
                c = lambda2(f)*50;
                [Xr,A,v] = HiLassoColMethod2(Y,D,prevAC,...
                                             lambda1(h),lambda2(f)/sqrt(N*K),...
                                             lambdaL,tol,H,max_iter,c);
                %se = separationError(X,Xr);
                save(fres,'A');
            else
                clear A;
                load(fres);
            end
            [Yo,Aols] = olsfun(Y,Do,A,reg); % experim
            clear Yo;
            Xo = cell(1,NC);
            for i=1:NC
                idx = find(groups == i);
                Xo{i} = D{i}*Aols(idx,:);
                v(i) = sum(sum((Aols(idx,:).^2)));
                %    Xr{i} = olsfun(Z,Do,A)
            end
            se = separationError(X,Xo);
            gA = group_act_set(Aols,K);
            hamm = mdlsHammingDistance(gAo,gA);
            eA = group_energy(Aols,K);
            fprintf(fh,'se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            fprintf('se=%g\thamm=%g\tact=%s\n',se,hamm,...
                    show_group_activity(mean(eA')));
            errorCHL(a,h,f) = se;
            hammCHL(a,h,f) = hamm;
            prevAC = A; 
            clear A Aols Xr Xo Yo;
        end % C-Hilasso
    end %lambda2 
end % lambda1
end % each sample combination
fclose(fh);

%
% summary
%  
sigma
l=[mean(hammL(:,:)) mean(errorL(:,:));...
 mean(hammL(1:8,:)) mean(errorL(1:8,:));...
 mean(hammL(9:end,:)) mean(errorL(9:end,:)) ]

gl=[mean(hammGL(:,:)) mean(errorGL(:,:)); ...
mean(hammGL(1:8,:)) mean(errorGL(1:8,:));...
mean(hammGL(9:end,:)) mean(errorGL(9:end,:))]

hl=[mean(hammHL(:,:)) mean(errorHL(:,:));...
mean(hammHL(1:8,:)) mean(errorHL(1:8,:));...
mean(hammHL(9:end,:)) mean(errorHL(9:end,:))]

cgl=[mean(hammCGL(:,:)) mean(errorCGL(:,:));...
mean(hammCGL(1:8,:)) mean(errorCGL(1:8,:));...
mean(hammCGL(9:end,:)) mean(errorCGL(9:end,:))]

chl=[mean(hammCHL(:,:)) mean(errorCHL(:,:));...
mean(hammCHL(1:8,:)) mean(errorCHL(1:8,:));...
mean(hammCHL(9:end,:)) mean(errorCHL(9:end,:))]

l1 = l(3,:);
l2 = l(2,:);

gl1 = gl(3,:);
gl2 = gl(2,:);

hl1 = hl(3,:);
hl2 = hl(2,:);

cgl1 = cgl(3,:);
cgl2 = cgl(2,:);

chl1 = chl(3,:);
chl2 = chl(2,:);

fprintf('1dig & %4.2f/%4.2f & %4.2f/%4.2f & %4.2f/%4.2f & %4.2f/%4.2f & %4.2f/%4.2f \\\\\n',...
        l1(2),l1(1),gl1(2),gl1(1),hl1(2),hl1(1),cgl1(2),cgl1(1),chl1(2),chl1(1));
fprintf('2dig & %4.2f/%4.2f & %4.2f/%4.2f & %4.2f/%4.2f & %4.2f/%4.2f & %4.2f/%4.2f \\\\\n',...
        l2(2),l2(1),gl2(2),gl2(1),hl2(2),hl2(1),cgl2(2),cgl2(1),chl2(2),chl2(1));
%%%

  %catch 
  %  fclose(fh);
  %end
