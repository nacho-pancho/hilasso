%
% function X=reconstruct(D,alpha,sparsity,threshold)
%
% reconstruct a set of patches, either exactly or approximately,
% by optionally specifying sparsity and/or threshold constrains
% on the reconstruction coefficients.
%
% The DC is specified separately.
%
% input:
% D ........ d x k, the codebook or dictionary consisting of k
%            d-dimensional atoms.
% alpha .... k x n reconstruction coefficients to express n
%            patches in terms of the k atoms of the dictionary.
% sparsity . if greater than 0, take only the 'sparsity' absolutly 
%            largest coefficients in the reconstruction of each
%            patch.
% threshold  Ignore the contribution of coefficients of absolute
%            magnitude below this value. The default value of 0
%            indicates no thresholding will be performed.
%
% output:
% Y ........ the reconstructed patches.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%

function Y=mdlsReconstruct(D,alpha,sparsity,threshold)
if ~exist('sparsity','var')
  sparsity=0;
end
if ~exist('threshold','var')
  threshold=0;
end
%
% not implemented yet
%if ~exist('maxerr','var')
%  maxerr=Inf;
%end
alpha=full(alpha);

[k,n]=size(alpha);
[d,k]=size(D);
%
% thresholding
%
if threshold > 0
  alpha=alpha.*(abs(alpha)>threshold);
end
%
% sparsity
%
if sparsity <= 0
  Y=zeros(d,n);
  return;
else
  if sparsity < k
    alpha=mdlsSparsify(alpha,sparsity);
  end
  Y=D*alpha;
end

