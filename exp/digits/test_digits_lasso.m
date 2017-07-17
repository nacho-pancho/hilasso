

Ns = 20;

Ys = Y(:,101:(100+Ns));
for ii=1:length(X)
    Xs{ii} = X{ii}(:,1:Ns);
end

for h=1:length(lambda)
    
%     [Xr,AL,v] = lassoMethod(Y,D,lambda(h));
%     errorL(h) = separationError(X,Xr);
%     detectedGroupsL{h} = v;
%     aSetL(h) = mean(sum(AL~=0));
    
    
%     Hilasso method
 %       [Xr,A] = groupLassoMethod(Ys,D,lambda(h));
%    errorG2(h) = separationError(Xs,Xr)
for j=1:length(lambda2)


    [Xr,A] = HIlassoMethod(Ys,D,lambda(h),lambda2(j));
    errorH(h,j) = separationError(Xs,Xr);
    aSetH(h,j) = mean(sum(A~=0));
    errorH
end
% 
%     % Group lasso method
%     detectedGroupsG{h} = v;
%     aSetG(h) = mean(sum(A~=0));
        
    
end