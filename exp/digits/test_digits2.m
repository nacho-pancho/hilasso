


sigma = [100];%:0.1:0.4;
rep = 3;


lambda1 = [0.1:0.4];
lambda2 = [0.6:1];


errorL2 = zeros(length(lambda1),1);
errorH = zeros(length(lambda1),1);
errorG = zeros(length(lambda1),1);
errorC = zeros(length(lambda1),length(lambda2));


aSetL2 = zeros(length(lambda1));
aSetG = zeros(length(lambda1));
aSetH = zeros(length(lambda1));
aSetC = zeros(length(lambda1),length(lambda2));



active = [4 5];


N = 100;
Ns = 30;


[Y,X] = createDataDigits(data,N,active);

Ys = Y(:,1:Ns);
for ii=1:length(X)
    Xs{ii} = X{ii}(:,1:Ns);
end




for h = 1:length(lambda1)

    % lasso method
    [Xr,AL] = lassoMethod(Y,D,lambda1(h));
    errorL2(h) = separationError(X,Xr);
    aSetL2(h) = mean(sum(AL~=0));

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
        [Xr,AC] = HIlassoColMethod(Y,D,lambda1(h),lambda2(f));
        errorC(h,f) = separationError(X,Xr);
        aSetC(h,f) = mean(sum(AC~=0));
        
        errorC
        errorL
        errorL2

    end

end


save test45 errorC errorL2 aSetL2 aSetC lambda1 lambda2

