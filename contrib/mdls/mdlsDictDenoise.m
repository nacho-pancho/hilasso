
function [imd psnr A] = mdlsDictDenoise(D,X,Dc,lambda,imo)

im_size = size(imo);
Y =  zeros(size(X));

for i=1:length(lambda)

    if length(D)>1
    n = mdlsDiscriminant(D,X,'lasso',10,1);
    [dum m] = min(n);
    else
        m = ones(1,size(X,2));
    end
    
    lambda(i)
    for ii=1:length(D)
        idx = find(m==ii);
%         A = scLasso(X, D{ii}, 64, lambda(i), 0);
        A=scOMP(X(:,idx),D{ii},64,lambda(i));
        Y(:,idx) = D{ii}*A;
    end


    imD{i} = reconstructImage(Y,7,im_size,Dc);


    mse = (imo-imD{i}).^2;
    mse = sum(mse(:))/im_size(1)/(im_size(2));


    psnr(i) = 10*log10(255^2/mse);
    
    psnr(i)

end

if length(lambda)==1
    imd = imD{1};
else
    imd = imD;
end