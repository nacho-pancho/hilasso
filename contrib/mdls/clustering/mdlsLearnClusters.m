% 
% function [clustering2 D] = mdlsLearnClusters(X,Do,k,T,lambda,mode,iter,l2)
%
% Input:
% X ............... data samples
% Do .............. initial dictionaries
% lambda .......... L1 penalty term
% xmu ............. Can be 0,1,2. Use 0 means use Do as initiali
%                   dictionary.
% iter ............ total iterations of the algorithm
%
% Output:
% clustering ...... cluster assignement for each sample in X
% D ............... cell of dictionaries for each cluster.
%
% @author Pablo Sprechmann (psprechmann@gmail.com)
%
function [clustering D] = mdlsLearnClusters(X,Do,lambda,xmu,iter,seed,true_labels,figfile)
if exist('seed','var')
  rand('twister',seed);
  randn('state',seed);
end
if ~exist('true_labels','var')
    true_labels = [];
end
if ~exist('xmu','var')
    xmu = 0; 
end
if ~exist('iter','var')
    iter = 10; 
end
if ~exist('figfile','var')
    n=round(10000*rand());
    figfile = sprintf('anonymous_clustering_vs_iters_%d.eps',n);
end

clustering = zeros(size(X,2),1);
D = Do; clear Do;
C = length(D);
if ~isempty(true_labels)
    N = length(true_labels);
    class_labels  = sort(unique(true_labels));
    class_labels2 = 1:length(class_labels);
    true_labels2 = zeros(1,length(true_labels));
    for j=1:N
        true_labels2(j) = find(true_labels(j)==class_labels);
    end
end

%
% default set of parameters for dictionary initialization
%
di_params  = mdlsDictIni();
di_params.sigma = 0.01;
di_params.layers = 1;
di_params.method = 'patches';
%di_params.method = 'trace'; % more expensive but better
%di_params.method = 'greedy'; % maximally incoherent greedy
di_params.K = zeros(1,C);
di_params.X = X;
for c=1:C
    di_params.K(c) = size(D{c},2);
end
%
% default set of parameters for dictionary learning
%
dl_params = mdlsLearnModel();
dl_params.mode = 2;
dl_params.max_iter = 100;
dl_params.batch_size = 0; % non-batch
dl_params.test_size = 0;
dl_params.mu_mode = [1 20]; % slowly add incoherence
dlmode = 2;
dl_params.xmu0 = xmu;
dl_params.training_data = X;
dl_params.resume = 0;
dl_params.discard_unused_atoms = 0.001;
dl_params.discard_constant_patches = 0;
dl_params.model.lambda = lambda;
 dl_params.batch_size = 4000;
%
% main loop
%
dmode = 'lasso';
err=zeros(1,iter);
for i=1:iter
    Dold = D;
    clu_old = clustering;
    %
    % Assign classes using Lasso energy
    %
    F = mdlsDiscriminant(D,X,dmode,lambda,0);
    [minValue,clustering] = min(F);
    if ~isempty(true_labels)
        cm   = mdlsConfusionMatrix(clustering, true_labels2);
        clu2 = mdlsFixUnsupervisedLabels(class_labels2, cm, clustering);
        cm   = mdlsConfusionMatrix(clu2, true_labels2);
        mdlsPrintConfusionMatrix(cm, class_labels2);
        err(i) = 100 - cm(end,end);
        mdlsFigure('Clustering error');
        plot(err(1:i),'*-'); grid on;
    end
    %dl_params.D0 = D;
    %D = mdlsLearnModel(dl_params);    
    %
    % re-learn dictionarioes from updated clusters
    %
    di_params.labels = clustering;
    dl_params.D0 = mdlsDictIni(di_params);

    dl_params.training_labels = clustering;
    D = mdlsLearnModel(dl_params);
    %
    % see if there is any change...
    %
    change = 0;
    if length(D) < C
        warning('Clusters collapsed into less classes!');
        D = Dold;
        clustering = ceil(C*randn(size(clustering))); % random assignement
        return;
    end
    for c=1:C        
        change = change + abs(norm(D{c}-Dold{c},'fro') / norm(Dold{c},'fro'));
    end
    fprintf('|D-Dold|/|Dold|=%f\n',change);
    if abs(change) < 1e-3
        break;
    end
end
mdlsFigure('Clustering error');
plot(err(1:i),'*-'); hold on; grid on;
h=xlabel('iterations');
p=get(gcf());
set(p.CurrentAxes,'Fontsize',16);
set(h,'FontSize',16);
h=ylabel('error');
set(h,'FontSize',16);
title('Clustering error vs. iterations');
n=round(10000*rand());
saveas(gcf(),figfile,'epsc');
hold off; close;