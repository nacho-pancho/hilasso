
function [I,A,error] = separateImages(Y1,Y2,Do,sp,lambda1,factor)


% extract patches
% Y = mdlsDeconstructFast(im{2}+im,sp,sp-1);
Y = Y1+Y2;

% normalize - [mejor no hacerlo!]
% for i=1:size(Y,2)
%     Y(:,i) = Y(:,i)/norm(Y(:,i));
% end

% perform collaborativ Hilasso
[Xr,A] = HIlassoColMethod(Y,Do,lambda1,factor);

% reconstructed images
for i=1:length(Xr)
    I{i} = reconstructImage(Xr{i},sp-1,size(im),zeros(size(Xr{1},2)));
end


error = separationError({Y1,Y2},Xr);



