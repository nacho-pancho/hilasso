function D = getDict(file,lambda,sp,ov,sD)


im = imread(file);
im = double(im);
im = imresize(im, 0.5);

[N1,N2] = size(im);
im = im(1:round(N1/2),1:round(N2/2));


% [X,Dc] = mdlsPatchify(im,sp,ov);
X = mdlsDeconstructFast(im,sp,ov);

% for i=1:size(X,2)
%     X(:,i) = X(:,i)/norm(X(:,i));
% end


Dini = mdlsGenPatchesDictionary(X,sD,0,1);
D =learnDict(X,[],{Dini},300,lambda,0,4000);