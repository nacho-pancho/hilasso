load exp/psprech/digits/dataUSPS.mat

Ngroups = [4,8,12];
% Ngroups = [4,8,12,20];%:0.1:0.4;
rep = 1;

lambda1 = [0.025,0.05,0.1,0.15];
% lambda2 =[0.2:0.1:0.7 1 2];
lambda2 =[0.25 0.5 0.75 1 1.25];
% sp = [4,8,12,16,32];
sp = [8,12,16];

errorL = zeros(length(lambda1),length(sp),length(Ngroups));
errorG = zeros(length(lambda1),length(sp),length(Ngroups));
errorC = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));
errorH = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups)); 

aSetL = zeros(length(lambda1),length(sp),length(Ngroups));
aSetG = zeros(length(lambda1),length(sp),length(Ngroups));
aSetC = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));
aSetH = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));

mseL = zeros(length(lambda1),length(sp),length(Ngroups));
mseG = zeros(length(lambda1),length(sp),length(Ngroups));
mseC = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));
mseH = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));

hammL = zeros(length(lambda1),length(sp),length(Ngroups));
hammG = zeros(length(lambda1),length(sp),length(Ngroups));
hammC = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));
hammH = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));


disp('echi')

V{1} = lambda1;
V{2} = lambda2;
V{3} = sp;
V{4} = Ngroups;


N = 200;
Ns = 100;
sizeD = 64;
m = 64;
n=m;
sigma = 0;
Nactive = 2;
clear error aux aux1
system('mkdir -p results/hilasso/test_direct');
for j=1:length(Ngroups)

    for i=1:length(sp)
        [Y,D,X,R,Ao] = createData(m,N,sp(i),Ngroups(j),sizeD,sigma,Nactive,0);
        Aos = Ao(:,1:Ns);
        Ys = Y(:,1:Ns);
        for ii=1:length(X)
            Xs{ii} = X{ii}(:,1:Ns);
        end

        
        for h = 1:length(lambda1)
            
            resfile = sprintf('results/hilasso/test_direct/LHG-ng%d-sp%d-lambda1%g.mat',...
                              Ngroups(j),sp(i),lambda1(h));
            if ~exist(resfile,'file')
                % lasso method
                [Xr,A] = lassoMethod(Y,D,lambda1(h));
                errorLt = separationError(X,Xr);
                mseLt = norm(full(Ao-A),'fro')^2/N/n;
                aSetLt = mean(sum(A~=0));
                hammLt = mdlsHammingDistance(Ao,A);
                
                [Xr,A] = HIlassoMethod(Ys,D,lambda1(h),lambda1(h));
                errorHt = separationError(Xs,Xr);
                mseHt = norm(full(Aos-A),'fro')^2/Ns/n;
                aSetHt = mean(sum(A~=0));
                hammHt = mdlsHammingDistance(Aos,A);
                
                [Xr,A] = groupLassoMethod(Ys,D,lambda1(h));
                errorGt = separationError(Xs,Xr);
                mseGt = norm(full(Aos-A),'fro')^2/Ns/n;
                aSetGt = mean(sum(A~=0));
                hammGt = mdlsHammingDistance(Aos,A);

                save(resfile,'errorLt','mseLt','aSetLt','hammLt',...
                     'errorHt','mseHt','aSetHt','hammHt',...
                     'errorGt','mseGt','aSetGt','hammGt');
            else
                load(resfile);
            end
            %
            % collaborative
            %
            for f=1:length(lambda2)
                resfile = sprintf('results/hilasso/test_direct/C-ng%d-sp%d-lambda1%g-lambda2%g.mat',...
                                  Ngroups(j),sp(i),lambda1(h),lambda2(f));
                if ~exist(resfile,'file')
                    % collaborative lasso method
                    [Xr,AC,v] = HIlassoColMethod(Y,D,lambda1(h),lambda2(f),0.00001);
                    detectedGroups{h}{f}{i}{j} = v;
                    
                    errorCt= separationError(X,Xr);
                    mseCt  = norm(full(Ao-AC),'fro')^2/N/n;
                    aSetCt = mean(sum(AC~=0));
                    hammCt = mdlsHammingDistance(Ao,AC);
                    save(resfile,'errorCt','mseCt','aSetCt','hammCt');
                else
                    load(resfile);
                end
                errorC(h,f,i,j)= errorCt;
                mseC(h,f,i,j)  = mseCt;
                aSetC(h,f,i,j) = aSetCt;
                hammC(h,f,i,j) = hammCt;
            end
            
            errorL(h,i,j)= errorLt;
            mseL(h,i,j)  = mseLt;
            aSetL(h,i,j) = aSetLt;
            hammL(h,i,j) = hammLt;
            
            errorH(h,i,j)= errorHt;
            mseH(h,i,j)  = mseHt;
            aSetH(h,i,j) = aSetHt;
            hammH(h,i,j) = hammHt;
            
            errorG(h,i,j)= errorGt;
            mseG(h,i,j)  = mseGt;
            aSetG(h,i,j) = aSetGt;
            hammG(h,i,j) = hammGt;                       
        end
        
        disp(['NG: ' num2str(Ngroups(j)) ' SP: ' num2str(sp(i))])
        fprintf('Lasso/Separation:');
        errorL(:,i,j)'
        fprintf('Lasso/Hamming:');
        hammL(:,i,j)'

        fprintf('GrLasso/Separation:');
        errorG(:,i,j)'
        fprintf('GrLasso/Hamming:');
        hammG(:,i,j)'

        fprintf('HiLasso/Separation:');
        errorH(:,i,j)'
        fprintf('HiLasso/Hamming:');
        hammH(:,i,j)'

        fprintf('C-HiLasso/Separation:\n');
        errorC(:,:,i,j)
        fprintf('C-HiLasso/Hamming:\n');
        hammC(:,:,i,j)
        
    end
    
    save exp/psprech/source/gridHilassoAux errorL aSet* error* mse* V detectedGroups lambda1 lambda2


end


save exp/psprech/source/gridHilasso64_B errorL aSet* error* mse* V