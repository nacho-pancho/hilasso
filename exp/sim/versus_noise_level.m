


Ngroups = 8;%[4,8,12,20];%:0.1:0.4;
rep = 1;

lambda1 = [0.1,,0.2,0.3 0.4 0.5];
lambda2 =[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

sp = [8];%[4,8,12,16,32];
sigmas = [0.05,0.1,0.2,0.4];

errorL = zeros(length(lambda1),length(sp),length(sigmas));
errorG = zeros(length(lambda1),length(sp),length(sigmas));
errorC = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas));
errorH = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas)); 

aSetL = zeros(length(lambda1),length(sp),length(sigmas));
aSetG = zeros(length(lambda1),length(sp),length(sigmas));
aSetC = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas));
aSetH = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas));

mseL = zeros(length(lambda1),length(sp),length(sigmas));
mseG = zeros(length(lambda1),length(sp),length(sigmas));
mseC = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas));
mseH = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas));

hammL = zeros(length(lambda1),length(sp),length(sigmas));
hammG = zeros(length(lambda1),length(sp),length(sigmas));
hammC = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas));
hammH = zeros(length(lambda1),length(lambda2),length(sp),length(sigmas));


disp('echi')

V{1} = lambda1;
V{2} = lambda2;
V{3} = sp;
V{4} = Ngroups;


N = 100;
Ns = 30;
sizeD = 64;
m = 64;
n=m;
Nactive = 2;
clear error aux aux1

for j=1:length(sigmas)

    for i=1:length(sp)
        [Y,D,X,R,Ao] = createData(m,N,sp(i),Ngroups,sizeD,sigmas(j),Nactive,0);
        Aos = Ao(:,1:Ns);
        Ys = Y(:,1:Ns);
        for ii=1:length(X)
            Xs{ii} = X{ii}(:,1:Ns);
        end
        % % direct method createData(n,N,k,numD,sizeD,sigma,Nactive,ort)
%         [Xr,A] = directMethod(Y,D,sigma(j));
%         errorD(i,j) = separationError(X,Xr);
%         aSetD(i,j) = mean(sum(A~=0));

        
        for h = 1:length(lambda1)
            
            % lasso method
            [Xr,A] = lassoMethod(Y,D,lambda1(h));
            errorL(h,i,j) = separationError(X,Xr);
            mseL(h,i,j) = norm(full(Ao-A),'fro')^2/N/n;
            aSetL(h,i,j) = mean(sum(A~=0));
            hammL(h,i,j) = mdlsHammingDistance(Ao,A);
            
            [Xr,A] = HIlassoMethod(Ys,D,lambda1(h),lambda1(h));
            errorH(h,i,j) = separationError(Xs,Xr);
            mseH(h,i,j) = norm(full(Aos-A),'fro')^2/Ns/n;
            aSetH(h,i,j) = mean(sum(A~=0));
            hammH(h,i,j) = mdlsHammingDistance(Aos,A);
            
            [Xr,A] = groupLassoMethod(Ys,D,lambda1(h));
            errorG(h,i,j) = separationError(Xs,Xr);
            mseG(h,i,j) = norm(full(Aos-A),'fro')^2/Ns/n;
            aSetG(h,i,j) = mean(sum(A~=0));
            hammG(h,i,j) = mdlsHammingDistance(Aos,A);
            
            for f=1:length(lambda2)

                % collaborative lasso method
                [Xr,AC] = HIlassoColMethod(Y,D,lambda1(h),lambda2(f),0.001);
                errorC(h,f,i,j) = separationError(X,Xr);
                mseC(h,f,i,j) = norm(full(Ao-AC),'fro')^2/N/n;
                aSetC(h,f,i,j) = mean(sum(AC~=0));
                hammC(h,f,i,j) = mdlsHammingDistance(Ao,AC);
                
            end

        end
        
        disp(['sigma: ' num2str(sigmas(j)) ' SP: ' num2str(sp(i))])
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
    
%     save exp/psprech/source/gridHilassoAux errorL aSet* error* mse* V


end


save exp/psprech/source/gridHilasso64_noise errorL aSet* error* mse* V