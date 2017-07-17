

function [Di,D,psnr] = evaluate(file,sigma)

% read image
im = imread( file );
im = double(im);
[N1,N2] = size(im);

% add noise
im = im + sigma * randn(N1, N2);

% extract patches
[X,Dc] = mdlsPatchify(im,8,7);


[S,Di2] = mdlsInitializeClusters(X,D{2}{1},4,sizeD,lambdaL1,sizeD,true);

[clustering Di2s101 Irr] = Texcluster(im,X,Di2,sizeD,7,lambdaL1,1);

[imd psnrS1 A] = mdlsDictDenoise(Di2s1,X,Dc,1.2*64*100,imo);

[clustering Di2s102 Irr] = Texcluster(im,X,Di2s101,sizeD,7,lambdaL1,1);

[imd psnrS2 A] = mdlsDictDenoise(Di2s82,X,Dc,1.2*64*100,imo);