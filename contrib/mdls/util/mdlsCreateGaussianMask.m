function W = mdlsCreateGaussianMask(width,sigma)
    W = zeros(width,width);
    if mod(width,2) == 1
        offset = ceil(width/2);
    else
        offset = width/2 + 0.5;
    end
    x = repmat((1:width)-offset,width,1);
    y = x';
    W = exp(-(x.^2+y.^2)*(1/(2*sigma^2)));
    W = W*(1/sum(W(:)));
end