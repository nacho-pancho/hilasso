
clear detectedGroups

sigma = [100];%:0.1:0.4;
rep = 3;


lambda1 = [0.5];
lambda2 = [1.6]/sqrt(230)/sqrt(300);


errorL2 = zeros(length(lambda1),1);
errorH = zeros(length(lambda1),1);
errorG = zeros(length(lambda1),1);
errorC = zeros(length(lambda1),length(lambda2));


aSetL2 = zeros(length(lambda1));
aSetG = zeros(length(lambda1));
aSetH = zeros(length(lambda1));
aSetC = zeros(length(lambda1),length(lambda2));


k = 2;
active = [3 5]; 
int = 7;

N = 150;
Ns = 30;

Nc = 1

     [Y,X] = createDataDigits(data,N,active,k,0);
for a=1:Nc
    

     [Y2,X2] = createDataDigits(data,30,[3 5 int],2,0);
     
     Y = [Y,Y2];
     
     for i=1:length(X)
         X{i} = [X{i} X2{i}];
     end
     

for h = 1:length(lambda1)

    % lasso method
%     [Xr,AL] = lassoMethod(Y,D,lambda1(h));
%     errorL2(h) = separationError(X,Xr);
%     aSetL2(h) = mean(sum(AL~=0));

    % Hilasso method
%     [Xr,A] = HIlassoMethod(Ys,D,lambda1(h),lambda1(h));
%     errorH(h) = separationError(Xs,Xr);
%     aSetH(h) = mean(sum(A~=0));
% 
%     % Group lasso method
%     [Xr,A] = groupLassoMethod(Ys,D,lambda1(h));
%     errorG(h) = separationError(Xs,Xr);
%     aSetG(h) = mean(sum(A~=0));

    for f=1:length(lambda2)


        % collaborative lasso method
        [Xr,AC,v] = HIlassoColMethod(Y,D,lambda1(h),lambda2(f));
        errorC(h,f) = separationError(X,Xr);
        detectedGroups{h}{f} = v;
        A{h}{f} = AC;
        aSetC(h,f) = mean(sum(AC~=0));
        
        errorC
        v

    end
    
    save('aux','errorC','aSetC','lambda1','lambda2','detectedGroups','Y','X','A')

end

file = ['Mix3digits' num2str(a)];
save(file,'errorC','aSetC','lambda1','lambda2','detectedGroups','Y','X','A')

end