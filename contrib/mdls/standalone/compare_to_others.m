D = load('-ascii','ksvd-dict.ascii');
D= mdlsDictNormalize(D);
addpath([WORKSPACE '/contrib/lars2']); % taken from MATLAB central
w = 8;
ov= 4;
I = double(imread('lena.pgm'));
[X,DC] = mdlsDeconstruct(I,w,ov);
%
% same indexes as chosen randomly in test_cpp
%
idx = [1506 612 2856 2479 3242] + 1; 
Xs = X(:,idx);
%mdlsFigure('Samples'); mdlsDictDisplay(Xs,'mag',2);
for i=2 %1:length(idx)
    fprintf('Sample #%d\n',idx(i));
    %
    % load solution path obtained with CSC
    %
    Acsc = load('-ascii',sprintf('A_%06d.ascii',i-1));
    state_csc = load('-ascii',sprintf('state_%06d.ascii',i-1));    
    J=size(Acsc,1); % number of iterations
    %
    % compute path using LARS in its three modes
    %
    usegram=1;
    Xj    = Xs(:,i);
    % the function expects normalized Xj
    Alars = lars(D,Xj/var(Xj),'LASSO');
    Alars = Alars*var(Xj); % un-scale
    %
    % compute path using OMP
    %
    Aomp = zeros(size(Acsc));
    for j=1:J
        Aj = scOMP(Xj,D,j,0.0);
        Aomp(j,:)= Aj';
    end
    figure(i);
    state_csc(:,1) = state_csc(:,1)/size(X,1);
    subplot(221); plot(1:J,Acsc(1:J,:)); axis([1 J 0 600]); 
    subplot(222); h=plot(1:J,state_csc(1:J,[1 3 7])); axis([1 J 0 1000]);
    legend(h,'MSE','Lr','L');
    subplot(223); plot(1:J,Alars(1:J,:)); axis([1 J 0 600]);
    subplot(224); plot(1:J,Aomp(1:J,:)); axis([1 J 0 600]);
end
tic
A=scOMP(X,D,15,0);
toc 