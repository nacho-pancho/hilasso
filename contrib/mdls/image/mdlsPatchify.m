%
% function [X,DC]=patchify(I,width,overlap)
%
% purpose:
% decompose an image as a series of patches, with an optional
% DC removal of each patch.
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
function [X,DC]=mdlsPatchify(I,width,overlap)
% 
% image domain dimensions are in capital case
%
[M,N]=size(I);
if nargin < 3 
    overlap=width/2;
end
remove_dc=0;
if nargout > 1
    remove_dc=1;
end

dw=width-overlap;
%
% sparse domain dimensions are in lower case
%
d=width*width;
n=ceil((M-width)/dw)*ceil((N-width)/dw);
% 
% create X from overlappedpatches
%
X=zeros(d,n);
DC=[];
if remove_dc
    DC=zeros(1,n);
end
step=width-overlap;
ii=1;
for i=1:dw:M-width+1
  for j=1:dw:N-width+1
    X(:,ii)=reshape(I(i:i+width-1,j:j+width-1),d,1);
    if (remove_dc)
        DC(ii)=sum(X(:,ii))/d;
        X(:,ii)=X(:,ii)-DC(ii);
    end
    ii=ii+1;
  end
end

