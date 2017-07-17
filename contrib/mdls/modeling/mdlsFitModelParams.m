%
% Given a dictionary D and a list of images, this function 
% computes the best fitting Laplacian, MOE and JOE parameters,
% globally and per-row, so that we can use them later in applications.
%
% INPUT
%
% D ................. dictionary 
% list_of_images .... cell with one filename per cell element. Optionally can
%                     contain path, but it must contain extension.
%
% OUTPUT
%
% stats ............. struct with all the estimated parameters
%
function stats=mdlsFitModelParams(prefix,D,list_of_images,img_dir)
    doMOL = true;
    doLASSO = true;
    debug = false;
    kappa0_row = 2.5;
    kappa0 = -1; 
    if iscell(D)
        D = D{1};
    end
    [M,K] = size(D);
    l2err = M*1.15*(1/255)^2;
    width = sqrt(M);
    overlap = width-1; 
    robust_par = 1;
    ignore_dc = 1;
    %
    % gather statistics
    %
    if ~exist('img_dir','var')
        img_dir = 'data/pascal06/train';
    end
    if ~exist('list_of_images','var')
        list_of_images = mdlsGetFileList(img_dir, 'pgm');
        aux = randperm(length(list_of_images));
        NI = min(100,length(list_of_images));
        list_of_images = list_of_images(aux(1:NI));
    end
    mu1_row = zeros(K,1);
    mu2_row = zeros(K,1);
    % per image averages: used for more robust MOLJ and JOE estimations
    NI = length(list_of_images);
    mu1_row_img = zeros(K,NI);
    mu2_row_img = zeros(K,NI);
    n1_row = zeros(K,1);
    n1_row_img = zeros(K,NI);
    cA = -1.5:(1/200):1.5;
    nPoints=length(cA);
    H_row = zeros(K,length(cA));
    maxAbsA = 0;
    N = 0;
    fprintf('Processing ');
    datafile = sprintf('%s-fitting.mat',prefix);
    img_time = 0;
    ii0 = 1;
    if exist(datafile,'file')
        load(datafile);
        ii0=ii+1;
    end
    err = zeros(M,1);
    %
    %==========================================================
    % Gather statistics
    %==========================================================
    %
    for ii = ii0:NI
        img = list_of_images{ii};
        fprintf('%4d/%4d:theta=',ii,NI);
        if ischar(img)
            I      = double(imread([img_dir '/' img]))*(1/255);
        else
            I = img; % if we want to train with an explicit image
        end
        if ignore_dc
            [X,DC] = mdlsDeconstruct(I,width,overlap);
        else
            X = double(mdlsDeconstructFast(I,width,overlap));
        end
        aux = randperm(size(X,2));
        Naux = min(1000,length(aux));
        X = X(:,aux(1:Naux));
        N = N + size(X,2);
        A      = scLasso(X,D,M-1,l2err,1); % min ||A||_1 s.t. ||X-DA|| <=
                                           % (1/255)^2
        E = X - D*A;
        err = err + mean(E.^2,2); 
        clear E;
        %
        % error parameters
        %
        
        for i = 1:K
            Ai  = A(i,:);
            nzi = find(Ai ~= 0);
            Anz = full(Ai(nzi));
            n1_row(i) = n1_row(i) + length(nzi);
            %
            % moments
            %
            mu1_row_img(i,ii) = sum(abs(Anz));
            n1_row_img(i,ii) = length(nzi);
            mu1_row(i) = mu1_row(i) + mu1_row_img(i,ii);
            mu2_row_img(i,ii) = sum(Anz.^2);
            if robust_par <= 0
                mu2_row(i) = mu2_row(i) + mu2_row_img(i,ii);
            else % saturate at robust_par
                mu2_row(i) = mu2_row(i) + sum(min(Anz.^2,robust_par^2));
            end
            maxAbsA = max(maxAbsA,max(abs(Anz)));
            %
            % histomgram
            %
            h = hist(Anz(:),cA);
            H_row(i,:) = H_row(i,:) + h;
        end
        fprintf('%6.3f\n',sum(n1_row)/sum(mu1_row));
        if mod(ii,50)==0 || (ii==NI)
            save(datafile,'ii','mu1_row','mu1_row_img','mu2_row','n1_row_img','n1_row','N','H_row','err');
        end
    end   
    fprintf(' [done]\n');
    %
    %==========================================================
    % parameter estimation
    %==========================================================
    %
    sigma2 = err/NI; 
    fprintf('Gaussian error sigma2=%f\tstdev=%f\n',mean(sigma2),std(sigma2));

    n1 = sum(n1_row);
    mu1 = sum(mu1_row)/n1;
    mu2 = sum(mu2_row)/n1;
    H = sum(H_row)*(1/n1); % global histogram
    for i = 1:K
        H_row(i,:) = H_row(i,:)/n1_row(i);
    end
    %
    % Laplacian
    %
    theta_row = n1_row./mu1_row; % per row
    theta  = 1/mu1; % global

    %
    % MOL with parameters given by ML
    %
    mu1_row = mu1_row./n1_row;
    mu2_row = mu2_row./n1_row;
    %    keyboard
    if kappa0 < 0
        kappa = 2*(mu2-mu1.^2)./(mu2-2*mu1.^2);
        if kappa < 2
            kappa = 2;
        end
    else
        kappa = kappa0;
    end
    beta = (kappa-1)*mu1;
    if kappa0_row < 0
        kappa_row = 2*(mu2_row-mu1_row.^2)./(mu2_row-2*mu1_row.^2);
    else
        kappa_row = repmat(kappa0_row, K, 1);
    end
    beta_row = (kappa_row-1).*mu1_row;
    %
    % Conditional Jeffreys: pick three random samples
    molj_kappa = max(3,min(3,NI));
    molj_kappa_row = molj_kappa*ones(K,1);
    
    theta_row_img =  n1_row_img./(mu1_row_img+eps); % avoid NaN
    for k=1:K
        aux = randperm(NI);        
        % we form MOLJ estimates taking the average of three 'samples'
        % but here the samples are themselves averages over a single
        % image each. Let's see if this gives more robust estimations
        molj_beta_row(k) = sum(1./theta_row_img(k,aux(1:min(3,NI))));
    end
    theta_row_img2 = theta_row_img(:);
    aux = randperm(length(theta_row_img2));
    molj_beta = sum( 1./theta_row_img2(aux(1:min(3,NI))) ); 
    clear aux;
    %
    % Constrained Jeffreys: pick empirical range of theta
    % also 'robust'
    theta_min = min(theta_row_img2); 
    theta_max = max(theta_row_img2);
    theta_min_row = min(theta_row_img')'; 
    theta_max_row = max(theta_row_img')';

    %
    %==========================================================
    % Plot the global fitting
    %==========================================================
    %
    hfg1 = mdlsFigure('Global fitting'); clf;
    semilogy(cA,H,'.','Color',[0,0.5,0]); 
    axis([min(cA) max(cA) 1e-8 1+1e-16])
    hold on
    indL = 1;
    leg = {};
    leg{indL} = 'Empirical distribution'; indL = indL + 1;

    %
    % Plot the laplacian
    %
    fprintf('Global Laplacian: theta = %0.2f\n',theta);
    f_lap = mdlsLapDisc(cA,theta);
    semilogy(cA,f_lap,'b')
    KL_laplace = mdlsKLDiv(H,f_lap);
    leg{indL} = sprintf('Laplacian (KL: %.3f)',KL_laplace); indL = indL + 1;

    %
    % MOL
    %
    f_mol = mdlsMOLDisc(cA,kappa,beta);
    fprintf('Global MOL: kappa = %3.1f beta = %5.2f\n',kappa,beta);
    semilogy(cA,f_mol,'r');
    KL_mol = mdlsKLDiv(H,f_mol);
    leg{indL} = sprintf('MOL (KL: %.3f)',KL_mol); indL = indL + 1;
    %
    % MOL-J : using conditional Jeffreys formulation
    %
    f_molj = mdlsMOLDisc(cA, molj_kappa, molj_beta);
    fprintf('Conditional Jeffreys: kappa = %3.1f beta = %5.2f\n',...
                 molj_kappa, molj_beta);
    semilogy(cA,f_molj,'m'); % magenta
    KL_molj = mdlsKLDiv(H,f_molj);
    leg{indL} = sprintf('MOLJ (KL: %.3f)',KL_molj); indL = indL + 1;

    %
    % Constrained Jeffreys (JOE)
    %
    symmet = true;
    f_joe = mdlsJoeDisc(cA,theta_min,theta_max,symmet);
    fprintf('Constrained Jeffreys: a = %3.1f b = %5.2f\n',...
                 theta_min, ...
                 theta_max);
    semilogy(cA,f_joe,'c'); % cyan
    KL_joe = mdlsKLDiv(H,f_joe);
    leg{indL} = sprintf('JOE (KL: %.3f)',KL_joe); indL = indL + 1;
    legend(leg);
    grid on;
    %
    %==========================================================
    % return results
    %==========================================================
    %
    stats = struct();
    %
    % Histogram
    %
    stats.hist_x = cA;
    stats.hist = H;
    stats.hist_row = H_row;
    %
    % variance of reconstruction error per row of X
    %
    stats.sigma2_row = sigma2;
    stats.sigma2 = mean(sigma2);
    %
    % Bernoulli 
    %
    stats.rho = n1*(1/K/N);
    stats.rho_row = n1_row*(1/N);
    %
    % Laplacian
    %
    stats.theta = theta;
    stats.theta_row = theta_row;
    %
    % MOE
    %
    stats.kappa = kappa;
    stats.beta = beta;
    stats.kappa_row = kappa_row;
    stats.beta_row = beta_row;
    %
    % JOE
    %
    stats.theta_max = theta_max;
    stats.theta_min = theta_min;
    stats.theta_max_row = theta_max_row;
    stats.theta_min_row = theta_min_row;
    %
    % MOL-J
    % 
    stats.molj_kappa = molj_kappa;
    stats.molj_beta = molj_beta;
    stats.molj_kappa_row = molj_kappa_row;
    stats.molj_beta_row = molj_beta_row;
end % function
