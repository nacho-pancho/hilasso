% 
% function [S,D] = mdlsInitializeClusters(X, Do, k, T,lambda,Dsize,normalized, seed)
%
% Input:
% X ............... data samples
% Do .............. initial dictionaries
% params .......... struct with configuration stuff:
%     k ............... if an initialization mode other than 0 is used, this
%                       is the number of desired classes
%     T ............... maximum number of non-zero coefficients in s.c.
%     lambda .......... L1 penalty term
%     mode ............ 1
%     Dsize ........... 
%
% Output:
% clustering ...... cluster assignement for each sample in X
% D ............... cell of initial dictionaries for each cluster.
% S ............... tree structure for hierarchical clustering.
%
% @author Pablo Sprechmann (psprechmann@gmail.com)
%
function [S,D] = mdlsInitializeClusters(X,Do,k,T,lambda,Dsize,normalized, ...
                                        seed,use_angle)
if nargin == 0
  params.k = 0;
  params.lambda = 0.1;
  params.Dsize = 100;
  params.normalized = true;
  params.seed = -1;
  params.mode = 1;
  params.T = 100;
  params.exp = false;
  params.local_sigma = false;
  S=params;
  D=[];
  return;
end
if exist('seed','var') && (seed > 0)
  rand('twister',seed);
  randn('state',seed);
end
if ~exist('use_angle','var')
  use_angle = false;
end

sD = size(Do,2);
[m,n] = size(X);


Xo = X;

% Initialization
MAXP=4000;
if size(X,2) > MAXP
  aux=randperm(size(X,2));
  X = X(:,aux(1:MAXP));
end
fprintf('Initial sparse coding\n');
A = scLasso(X, Do ,T,lambda,1);
if normalized
   for jj=1:size(A,2)
      A(:,jj) = A(:,jj)*(1 / norm(A(:,jj)));
   end   
end
if ~use_angle
  K = abs(A)'*abs(A);%.*As;%.*(X'*X);
else
  K = abs(A'*A);%.*(X'*X);
end
%     K2 = X'*X;
%     K = K2.*K;
fprintf('Spectral clustering\n');
% EXPERIMENTAL: WEN-YEN CHEN
[clustering,centers,eigvecs,eigvals] = GD_SpectralClustering(K,2,1);
%[clustering,dum1,dum2,dum3] = sc(K,-1,2);
for ii = 1:2
    idx{ii} = find(clustering == ii);

end

if length(idx{2}) > length(idx{1})
    idxA = idx{1};
    idx{1} = idx{2};
    idx{2} = idx{1};
end

S(1) = createDic(X(:,idx{1}), n, sD, T, lambda, Do, Dsize, k==2, [], [], use_angle);
S(2) = createDic(X(:,idx{2}), n, sD, T, lambda, Do, Dsize, k==2, [], [], use_angle);


if k == 2
    D{1} = S(1).D;
    D{2} = S(2).D;
    return
end

for i = 1:k-2
    
    Vsave = zeros(1,length(S));
    for j = 1:length(S)
        
        
        if isempty(S(j).splitD)
            Saux = createDic(S(j).X, n, sD, T, lambda, Do, Dsize, 0, S(j).D, ...
                             S(j).costD, use_angle);
            S(j) = Saux;
            clear Saux
        end
        Vsave(j) = S(j).save;
    end
    
    disp(Vsave)
    [m,id] = max(Vsave);
    
    Saux = S(id);
    
    S(id).X = Saux.splitX{1}; 
    S(id).D = Saux.splitD{1};
    S(id).splitD = [];
    S(id).costD = Saux.costsplitD{1};
    
    jj = length(S) + 1;
    S(jj).X = Saux.splitX{2}; 
    S(jj).D = Saux.splitD{2};
    S(jj).splitD = [];
    S(jj).costD = Saux.costsplitD{2};
    
end

for i = 1:k
    D{i} = S(i).D;
end

%==========================================================================

function S = createDic(X,n,sD,T,lambda,Do,Dsize,twoClasses,D,costD,use_angle)
fprintf('d');
L =  size(X,2);
% Dsize = 150;

% Create set dictionary
if ~exist('D','var') || isempty(D)
%     Dini = mdlsGenPatchesDictionary(X,round(sD*L/n),0,1);
    Dini = mdlsGenPatchesDictionary(X,Dsize,0,1);
    %[D,a] = scDictLearn(X,Dini,60,lambda,2);
    [D,a] = learnDict(X,Dini,Dsize,lambda);
    costD = cost(X,D,a,lambda);
end

if not(twoClasses)

    % Split
    A = scLasso(X, Do ,T,lambda,1);
    if ~use_angle
        K = abs(A)'*abs(A);
    else
        K = abs(A'*A);
    end

    fprintf('Spectral clustering\n');
    [clustering,centers,eigvecs,eigvals] = GD_SpectralClustering(K,2,1);
    %[clustering,dum1,dum2,dum3] = sc(K,-1,2);
% fprintf('Spectral clustering\n');
% [clustering,centers,eigvecs,eigvals] = GD_SpectralClustering(K,2,1);
% %[clustering,dum1,dum2,dum3] = sc(K,-1,2);

    for ii = 1:2
        idx{ii} = find(clustering == ii);

    end

    X1 = X(:,idx{1});
    if length(idx{1}) > Dsize

        % Compute first dictionary
        % Dini = mdlsGenPatchesDictionary(X1,round(sD*length(idx{1})/n),0,1);
        Dini = mdlsGenPatchesDictionary(X1,Dsize,0,1);
        % [D1,a1] = scDictLearn(X1,Dini,60,0.1,2);
        [D1,a1] = learnDict(X1,Dini,Dsize,0.1);

        costD1 = cost(X1,D1,a1,lambda);

    else
        costD1 = 2*costD;
        D1 = zeros(size(X,1),Dsize);
    end

    X2 = X(:,idx{2});
    if length(idx{2}) > Dsize
        % Compute second dictionary
        % Dini = mdlsGenPatchesDictionary(X2,round(sD*length(idx{2})/n),0,1);
        Dini = mdlsGenPatchesDictionary(X2,Dsize,0,1);
        %[D2,a2] = scDictLearn(X2,Dini,60,0.1,2);
        [D2,a2] = learnDict(X2,Dini,Dsize,0.1);

        costD2 = cost(X2,D2,a2,lambda);

    else
        costD2 = 2*costD;
        D2 = zeros(size(X,1),Dsize);
    end


    saveSplit = costD - costD1 - costD2;

else
    X1 = [];
    X2 = [];
    D1 = [];
    D2 = [];
    costD1 = [];
    costD2 = [];
    costD = [];
    saveSplit = 0;
end
    
S = struct('X',X,...
           'D',D,...
           'splitX',{{X1,X2}},...
           'splitD',{{D1,D2}},...
           'save',saveSplit,...
           'costD',costD,...
           'costsplitD',{{costD1,costD2}});

%==========================================================================

function c = cost(X,D,A,lambda)
fprintf('c');


g = X - D*A;
for j=1:size(A,2)
    normasError(j) = norm( g(:,j) )^2;
end

c = sum( normasError+lambda*sum(abs(A),1) );

function [D,A]=learnDict(X,Dini,J,lambda)
  dl_params = mdlsLearnModel();
  dl_params.training_data = X;
  dl_params.D0 = Dini;
  dl_params.lambda0 = lambda*sqrt(size(X,1));
  dl_params.xmu0 = 0;
  dl_params.mu0 = 0;
  dl_params.batch_size = 0;
  dl_params.test_size = 0;
  dl_params.resume = false;
  D = mdlsLearnModel(dl_params);
  D = D{1};
  L = ceil(size(X,1)/2);
  A = scLasso(X,D,L,lambda,2);
