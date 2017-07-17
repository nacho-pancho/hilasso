%
% function I = mdlsReconstruct(X,m,n,overlap)
%
% Given a set of patches X corresponding to an m x n image
% and the overlap that was used to obtain the patches from it,
% this function reconstructs the image (or a modified version of it).
%
% Input:
% 1 X ............ patches (n x N)
% 2 m ............ height of output image
% 3 n ............ width of output image
% 4 overlap ...... overlap of patches
% 5 DC ........... optional, DC of the patches
%
% output:
% 1 I ............ Output image.
%
