%
% purpose:
% compute the two-part code length of an image given as a set
% of non-overlapping patches, given the dictionary, the
% reconstruction coefficients and  the hiperparameters lambda 
% (penalty) and s (sparsity level).
% 
% The code is formed as L(y|D,alpha,lambda,sparsity) +
% L(alpha|lambda,sparsity) + L(lambda) + L(sparsity) + L(D) 
%
% the first term actually includes a sepparate term for describing
% the DC of each patch. This is recovered exactly since it can
% take at most d*A values, with A=256 and d=w*w the dim of the
% patches. This takes log(d*A) bits per patch.
%
% This is the base case implicit in most sparsity problems with
% an L2 prior on the reconstruction error, which in turn implies
% a m-dimensional symmetric Gaussian prior on the error
% (covariance matrix is diagonal with sigma^2=||E||_F2) and zero
% mean.
%
% input:
%
% X ....... patches of the image
% D ....... reconstruction dictionary
% alpha ... rec. coefficients
% lambda .. rec. penalty (L1 prior parameter)
% s ....... desired sparsity
% precision_D 
% ......... desired precision of the dictionary coefficients, given
%           in bits. 0 would mean some MDL-optimal value, but this is not explored yet. 
% precision_alpha 
% ......... desired precision of the alpha coefficients, given
%           in bits. 0 would mean some MDL-optimal value, but
%           this is not explored yet.
% precision_error
% ......... desired precision of the error coefficients, 
%           in bits. 0 would mean some MDL-optimal value, but
%           this is not explored yet.
%
% output:
%
% LDC .... Description length of the DC component of the patches
% LE ..... Desc. length of the error term, namely p(X|D,alpha,lambda) for X without the DC
% LA ..... Desc. length associated with the reconstruction coefficients alpha
% LD ..... Desc. length of the Dictionary.
%
% if only one output value are requested, it will be the sum of the above four,
% i.e., the total length.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [LDC,LE,LA,LD]=mdlsTwoPartCode(X,XDC,D,alpha,lambda,sparsity, ...
  precision_D, precision_alpha, precision_error, ... 
  prior_E, prior_alpha, prior_D)

% 
% priors for each coding element
%
if ~exist('prior_E','var')
  prior_E='Gaussian';
end
if ~exist('prior_D','var')
  prior_D='Ball';
end
if ~exist('prior_alpha','var')
  prior_alpha='Laplacian';
end

if (nargout ~=  1) && (nargout ~= 4)
  error('Must have either one or four output arguments');
end
%
% alphabet size
%
asize=256;
[d,n]=size(X);

fprintf('lambda=%f sparsity=%d precision_alpha=%d precision_D=%d\n',lambda,sparsity,precision_alpha,precision_D);
%
% L(DC)
%
LDC=mdlsDCDL(XDC,1/d,'Uniform',asize);
%
% approximate signal with quantized versions of D and alpha
%
D=round(D.*2^precision_D)/(2^precision_D);
alpha=round(alpha.*2^precision_alpha)/(2^precision_alpha);
Y=mdlsReconstruct(D,alpha,sparsity);
%
% L(y|D,alpha,lambda,sparsity) = L(E)
%
E=round((X-Y)*2^precision_error)/2^precision_error;
fprintf('max E=%f\n',max(max(abs(E))));
LE=mdlsErrorDL(E,precision_error,prior_E,0);
%
% L(alpha)
%
LA=mdlsReconstructionDL(alpha,asize,precision_alpha,sparsity,prior_alpha,0);
%
%
% L(D)
%
LD=mdlsDictionaryDL(D,precision_D,prior_D); 
%
%
% finally, convert cost to bits
%
LDC=LDC/log(2);
LA=LA/log(2);
LE=LE/log(2);
fprintf('n=%6d d=%2d LDC=%10d LE=%10d LA=%10d LD=%10f \n',n,d,ceil(LDC),ceil(LE),ceil(LA),ceil(LD));
L=LDC+LE+LA+LD;
if nargout == 1
  LDC = L;
end