% 
% draw samples from a GAMMA of LAPLACIAN distribution.
% See mdlsGamLapSample for a description of this distribution.
% 
% inputs: 
% N ........ number of samples to draw
% k  ....... shape parameter of Gamma prior
% b  .......  beta parameter of Gamma prior
%
% outputs:
% x ........ a vector of N samples from this distribution.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function x=mdlsGamLapSample(N,k,b)
%
% use the X=F^-1(A) technique with A chosen uniformly between 0
% and 1
%
a=rand(1,N);
x=b*( (1-a).^(-1/k) - 1 );
