function [Y,I]=getTestImage(file,N,M,sp,ov)

im = imread(file);
im = double(im);
im = imresize(im, 0.25);

[N1,N2] = size(im);
im = im(1:round(N1/2),1:round(N2/2));
I = im(end-N+1:end,end-M+1:end);


Y = mdlsDeconstructFast(I,sp,ov);