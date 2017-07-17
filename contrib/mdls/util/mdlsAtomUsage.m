%
% function F=mdlsAtomUsage(A,mode)
%
% Category: utility function
%
% Description: returns the relative frequency with which each atom 
%              is used for a given data. This is given by the
%              reconstruction coefficients A.
%              The algorithm has a 'soft' and a 'hard' mode.
%              In the first one, the absolute values of the coefficients
%              related to each atoms are added up. In the hard mode,
%              every time a coefficient is nonzero, the 
%              usage count for the corresponding atom is incremented in 1.
%
% Input:
% A ........ coefficients (K x N)
% mode ..... 'hard' or 'soft'. Defaults to 'soft'
% threshold . only meaningful for mode='hard', counts occurence if absolute value
%             of coefficient is greater than this value. Defaults to 1e-5.
% Output:
% F ........ atom usage as a Kx1 column vector
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Idn$
%
function F=mdlsAtomUsage(A,mode,threshold)
  if ~exist('mode','var')
    mode = 'soft';
  end
  if ~exist('threshold','var')
    threshold = 1e-5;
  end

  [K,N]=size(A);

  if isequal(mode,'soft')
     F = sum(abs(A),2);
     F = F * (1/max(F));
  else
     F = sum((abs(A) > threshold) , 2) * (1/N);
  end 
end
