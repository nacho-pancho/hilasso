%
% function [D,alpha,dc]=mdlsEncode(I,width,overlap, K, sparsity,lambda,maxerr,metric)
%
% purpose:
%    take an image and encode it using a sparseland model.
% The user must specify the width of the patches, the overlap 
% the size of the dictionary (to be trained specifically for this
% image), the desired sparsity and a penalty term (lambda).
%
% input:
%
% I ......... image to be encoded
% width ..... width of the patches
% overlap ... desired overlap
% K ......... size of the dictionary
% sparsity .. desired sparsity in reconstruction coefficients
% lambda .... penalty term
% maxerr .... only useful in OMP based reconstruction, specify
%             maximum desired error variance.
%
% output:
% D ......... dictionary created for this image
% alpha ..... reconstruction coefficients
% dc ........ DC components of each patch
%
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [D,alpha,dc]=mdlsEncode(I,width,overlap, K, sparsity,lambda,maxerr,metric)

if ~exist('metric','var')
  metric=2;
end
if ~exist('maxerr','var')
  maxerr=0;
end

if nargin < 6
  error('Must specify at least first 6 parameters');
end

[x,dc]=patchify(I,width,overlap);
D=dicIni(x,K,10,0); % add some gaussian noise with sigma = 10
if metric == 0
  D=mexDL(x,D,10,lambda,2);
  alpha=mexOMP(x,D,sparsity,maxerr);
elseif metric == 1
  [D,alpha]=dicLearnL1(x,D0,10,lambda,1e-4);
  alpha=sparseCodingL1(x,D,sparsity,lambda);
else
  D=mexDL(x,D,10,lambda,2);
  alpha=mexLasso(x,D,sparsity,lambda,2);
end

