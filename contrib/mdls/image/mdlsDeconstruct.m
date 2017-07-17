%
% function [X,DC]=mdlsDeconstruct(I,width,overlap)
%
% purpose:
% decompose an image as a series of patches. In this version,
% the DC of each patch is removed. If you don't want this to happen
% and care about speed, use mdlsDeconstructFast instead.
%
% input:
% I ......... an image to be decomposed
% width ..... width of the patches
% overlap ... overlap of the patches (both vert and horiz)
%
% output:
%
% X ......... patches given as a d x n matrix, where
%             d=width*width is the dimension of the patches
% DC ........ DC of each patch given as a 1 x n vector.
%
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
