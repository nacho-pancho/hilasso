


Ngroups = [8];
sizeD = 32;
Nactive = 2;
sp = 3;
rep = 1;
N = 1;
m = 128;

sigma = 0.5;
lambda1 = 2.^(-3:-2);
lambda2 = [0.25 0.5 0.75];
randn('state',12345897);
rand('twister',12345897);

[Y,D,X,R,Ao] = createData(m,N,sp,Ngroups,sizeD,sigma,Nactive,0);
XT=X{1};
for g=2:length(X)
  XT = XT + X{g};
end
size(Y)    
close all
if ~exist('results/sim','file')
    mkdir('results/sim');
end
for h = 1:length(lambda1)
    fprintf('lambda1=%g lambda2=%g\n',lambda1(h),0);
    fname = sprintf('results/sim/simple-lambda1=%g-lambda2=%g-',lambda1(h),0);
    % lasso method
    [Xr,A] = lassoMethod(Y,D,lambda1(h));
    %    errorL(h,i,j) = separationError(X,Xr);
    %    mseL(h,i,j) = norm(full(Ao-AL),'fro')/N;
    %    aSetL(h,i,j) = mean(sum(AL~=0));
    [group_energy(Ao,sizeD)';...
    group_energy(A,sizeD)']

    figure(100*h); 
    stem([Ao A]);; legend('A0','A');
    title(sprintf('lambda1=%g lambda2=%g',lambda1(h),0));
    plot2svg([fname '-stem.svg'],gcf());
    %    keyboard
    for f=1:length(lambda2)
        fprintf('lambda1=%g lambda2=%g\n',lambda1(h),lambda2(f));
        fname = sprintf('results/sim/simple-lambda1=%g-lambda2=%g-', ...
                        lambda1(h),lambda2(f));        
        tic
            [Xr,A] = HIlassoMethod(Y,D,lambda1(h),lambda2(f),1e-6);
        toc
        %        errorH(h,i,j) = separationError(Xs,Xr);
        %        mseH(h,i,j) = norm(full(Ao-A),'fro')/Ns;
        %        aSetH(h,i,j) = mean(sum(A~=0));
        [group_energy(Ao,sizeD)';...
        group_energy(A,sizeD)']
        
        figure(100*h+f); 
        fh=stem([Ao A],'MarkerSize',10




,'LineWidth',1); 
        set(fh(1),'Color','red');
        set(fh(1),'Marker','x');
        set(fh(2),'Color','blue');
        legend('A0','A','Location','South');
        title(sprintf('lambda1=%g lambda2=%g',lambda1(h),lambda2(f)));
        axis([44 164 -1 1])
        plot2svg([fname 'stem.svg'],gcf());
    end
end        
