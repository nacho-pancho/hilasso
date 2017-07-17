


Ngroups = 8;%[4,8,12,20];%:0.1:0.4;
rep = 1;

lambda1 = [0.01,0.05,0.1,0.15,0.2,0.3];
lambda2 =[0.2:0.1:0.7 1 2];
sp = [5,8];%[4,8,12,16,32];

errorL = zeros(length(lambda1),length(sp),length(Ngroups));
errorC = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));
errorH = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups)); 

aSetL = zeros(length(lambda1),length(sp),length(Ngroups));
aSetC = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));
aSetH = zeros(length(lambda1),length(lambda2),length(sp),length(Ngroups));


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
sigma = [0.05,0.1,0.2];
Nactive = 2;
clear error aux aux1

for j=1:length(sigma)

    for i=1:length(sp)
        [Y,D,X,R,Ao] = createData(m,N,sp(i),Ngroups,sizeD,sigma(j),Nactive,0);
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
            [Xr,AL] = lassoMethod(Y,D,lambda1(h));
            errorL(h,i,j) = separationError(X,Xr);
            mseL(h,i,j) = norm(full(Ao-AL),'fro')^2/N/n;
            aSetL(h,i,j) = mean(sum(AL~=0));
            
            [Xr,A] = HIlassoMethod(Ys,D,lambda1(h),lambda1(h));
            errorH(h,i,j) = separationError(Xs,Xr);
            mseH(h,i,j) = norm(full(Aos-A),'fro')^2/Ns/n;
            aSetH(h,i,j) = mean(sum(A~=0));
            
            [Xr,A] = groupLassoMethod(Ys,D,lambda1(h));
            errorG(h,i,j) = separationError(Xs,Xr);
            mseG(h,i,j) = norm(full(Aos-A),'fro')^2/Ns/n;
            aSetG(h,i,j) = mean(sum(A~=0));

            
            for f=1:length(lambda2)

                % collaborative lasso method
                [Xr,AC] = HIlassoColMethod(Y,D,lambda1(h),lambda2(f),0.001);
                errorC(h,f,i,j) = separationError(X,Xr);
                mseC(h,f,i,j) = norm(full(Ao-AC),'fro')^2/N/n;
                aSetC(h,f,i,j) = mean(sum(AC~=0));
                
            end

        end
        
        disp(['NG: ' num2str(sigma(j)) ' SP: ' num2str(sp(i))])
        errorL(:,i,j)
        errorG(:,i,j)
        errorH(:,i,j)
        errorC(:,:,i,j)
        
    end
    
%     save exp/psprech/source/gridHilassoAux errorL aSet* error* mse* V


end


save exp/psprech/source/gridHilasso64_noise errorL aSet* error* mse* V