%
% function X=mdlsDeconstructFast(I,width,overlap)
%
% purpose:
% Decompose an image as a series of patches. 
%
% input:
% I ......... an image to be decomposed
% width ..... width of the patches
% overlap ... overlap of the patches (both vert and horiz). Defaults to width-1
%
% output:
%
% X ......... patches given as a d x n matrix, where
%             d=width*width is the dimension of the patches
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
