% 
% psnr = mdlsReconstructionError(X,Xr,varargin)
%
% Name: mdlsReconstructionError
%
% Category:
%
% Description: Given a signal X and an approximation Xr, compute the error between the two
%              and various statistics.
%
% Input:
% X ........ original data
% Xr ....... approximation
%
% Output:
% psnr ..... 
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%

function errorStats = mdlsReconstructionError(X,Xr,varargin)
showError = false;
for arg = 1:2:length(varargin)
  if isequal(varargin{arg}, 'ShowError')
    showError = varargin{arg + 1};
  end
end

[n,N] = size(X);
E = X - Xr; 
MSE = sum(sum(E.^2))*(1/n/N);
psnr = -10*log10(MSE);

errorStats.mse = MSE;
errorStats.psnr = psnr;
errorStats.mean = mean(E(:));
errorStats.median = median(E(:));
errorStats.max = max(E(:));
errorStats.min = min(E(:));
errorStats.maxabs = max(abs(E(:)));

if showError 
  imagesc(E)
  %   axis image
  axis off
  colorbar
end

end % function
